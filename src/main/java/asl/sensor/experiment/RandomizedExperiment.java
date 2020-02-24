package asl.sensor.experiment;

import static asl.utils.FFTResult.DEFAULT_TAPER_WIDTH;
import static asl.utils.NumericUtils.TAU;
import static asl.utils.NumericUtils.atanc;
import static asl.utils.NumericUtils.complexRealsFirstSorter;
import static asl.utils.NumericUtils.getComplexSDev;
import static asl.utils.NumericUtils.multipointMovingAverage;
import static asl.utils.NumericUtils.rewrapAngleDegrees;
import static asl.utils.NumericUtils.unwrap;
import static asl.utils.NumericUtils.unwrapArray;
import static asl.utils.ReportingUtils.complexListToString;
import static asl.utils.ReportingUtils.complexListToStringWithErrorTerms;
import static asl.utils.ResponseUnits.SensorType.TR120;
import static asl.utils.ResponseUnits.SensorType.TR240;
import static asl.utils.ResponseUnits.SensorType.TR360;
import static asl.utils.ResponseUnits.getFilenameFromComponents;
import static asl.utils.TimeSeriesUtils.ONE_HZ_INTERVAL;
import static asl.utils.input.InstrumentResponse.loadEmbeddedResponse;

import asl.sensor.input.DataStore;
import asl.utils.FFTResult;
import asl.utils.ResponseUnits.ResolutionType;
import asl.utils.ResponseUnits.SensorType;
import asl.utils.input.DataBlock;
import asl.utils.input.InstrumentResponse;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresBuilder;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresOptimizer;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresProblem;
import org.apache.commons.math3.fitting.leastsquares.LevenbergMarquardtOptimizer;
import org.apache.commons.math3.fitting.leastsquares.MultivariateJacobianFunction;
import org.apache.commons.math3.fitting.leastsquares.ParameterValidator;
import org.apache.commons.math3.linear.DiagonalMatrix;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

/**
 * This experiment takes in a randomized calibration signal and the corresponding output from a
 * seismic sensor. It calculates the implicit response by deconvolving the calibration signal from
 * the sensor output, and then finds the best-fit poles (lowest 2 for low-frequency calibrations,
 * all remaining poles for high-frequency calibrations) to match the magnitude and rotation angle of
 * the calculated response curve produced from the deconvolution. Plottable data includes the
 * sensor's response curve (Bode plot), the calculated response from deconvolution, and the plot of
 * the response from the best-fit parameters. These plots are returned in a list of two datasets:
 * the first holds the absolute values per-frequency of those curves, and the second holds the angle
 * of each such point in complex space. For more details on the algorithm, see Ringler, Hutt, et
 * al., "Estimating Pole-Zero Errors in GSN-IRIS/USGS Network Calibration Metadata", Bulletin of the
 * Seismological Society of America, Vol 102 (Apr. 2012).
 *
 * @author akearns - KBRWyle
 */
public class RandomizedExperiment extends Experiment {

  /**
   * Maximum possible frequency bound as a multiple of Nyquist rate of input (90%). Also used to set
   * the maximum value of plottable data for high-frequency calibrations.
   */
  public static final double PEAK_MULTIPLIER = 0.9;

  /**
   * Default Nyquist rate limit for high-frequency calibrations, which usually preserves enough
   * corner freq. information without becoming too susceptible to noise. Determined by
   * experimentation.
   */
  public static final double DEFAULT_NYQUIST_PERCENT_LIMIT = 0.5;

  /**
   * Minimum nyquist rate percent for high-frequency calibrations.
   */
  public static final double MIN_MULTIPLIER = 0.3;

  /**
   * Sets the default normalization point for low-frequency calibration data (0.02 Hz)
   */
  private static final double LOW_FREQ_ZERO_TARGET = 0.02;

  /**
   * Sets the default normalization point for high-frequency calibration data (1.0 Hz) This value is
   * subject to change based on verification of cal fitting performance.
   */
  private static final double HIGH_FREQ_ZERO_TARGET = 1.0;

  /**
   * Choice of correction factor to use
   */
  public static final SensorType[] VALID_CORRECTIONS = {null, TR120, TR240, TR360};

  private InstrumentResponse correctionToApply;
  private double initialResidual, fitResidual;
  private List<Complex> initialPoles;
  private List<Complex> fitPoles;
  private Map<Complex, Complex> poleErrors;
  private List<Complex> initialZeros;
  private List<Complex> fitZeros;
  private Map<Complex, Complex> zeroErrors;
  /**
   * True if calibration is low frequency. This affects which poles are fitted, either low or high
   * frequencies.
   */
  private boolean isLowFrequencyCalibration;
  /**
   * Used to ensure that scaling is done correctly on calibrations that use a capacitive setting.
   * Default value will be false.
   */
  private boolean isCapacitive;
  private InstrumentResponse fitResponse;
  private double[] freqs;
  private boolean plotUsingHz;
  private double maxMagWeight, maxArgWeight; // max values of magnitude, phase
  private double nyquistMultiplier; // region up to nyquist to take for data
  private int numIterations; // how much the solver ran

  public RandomizedExperiment() {
    super();
    isCapacitive = false;
    isLowFrequencyCalibration = false;
    numIterations = 0;
    plotUsingHz = true;
    nyquistMultiplier = DEFAULT_NYQUIST_PERCENT_LIMIT; // defaults to 0.8
  }

  /**
   * Get the frequency at which the data will be normalized.
   *
   * @return Correct frequency to normalize data at for the given type of calibration curve
   * @see #LOW_FREQ_ZERO_TARGET
   * @see #HIGH_FREQ_ZERO_TARGET
   */
  static double getFrequencyForNormalization(boolean isLowFrequencyCalibration) {
    return isLowFrequencyCalibration ? LOW_FREQ_ZERO_TARGET : HIGH_FREQ_ZERO_TARGET;
  }

  private static Complex[] getResponseCorrection(double[] frequencies, InstrumentResponse resp) {
    if (resp == null) {
      Complex[] returnValue = new Complex[frequencies.length];
      Arrays.fill(returnValue, Complex.ZERO);
      return returnValue;
    }

    return resp.applyResponseToInputUnscaled(frequencies);
  }

  /**
   * Backend function to set instrument response according to current test variables (for best-fit
   * calculation / backward difference) and produce a response from that result. The passed response
   * is copied on start and is not modified directly. Which values (poles) are modified depends on
   * high or low frequency calibration setting.
   *
   * @param variables Values to set the response's poles to
   * @param freqs Set of frequencies to get the response curve over
   * @param numZeros How many (paired) variables represent zeros (to determine first pole index)
   * @param fitResponse Response to apply these variables to
   * @param isLowFreq True if the calibration being fit to is low-frequency
   * @param correctionCurve Response correction used for some Trillium responses
   * @return Doubles representing new response curve evaluation
   */
  private static double[] evaluateResponse(double[] variables, double[] freqs, int numZeros,
      InstrumentResponse fitResponse, boolean isLowFreq, Complex[] correctionCurve) {

    InstrumentResponse testResp = new InstrumentResponse(fitResponse);

    // prevent terrible case where, say, only high-freq poles above nyquist rate
    if (variables.length > 0) {
      testResp = fitResponse.buildResponseFromFitVector(
          variables, isLowFreq, numZeros);
    } else {
      System.out.println("NO VARIABLES TO SET. THIS IS AN ERROR.");
    }

    Complex[] appliedCurve = testResp.applyResponseToInputUnscaled(freqs);
    double[] curValue = new double[freqs.length * 2];

    if (correctionCurve != null) {
      for (int i = 0; i < appliedCurve.length; ++i) {
        appliedCurve[i] = appliedCurve[i].subtract(correctionCurve[i]);
      }
    }

    for (int i = 0; i < freqs.length; ++i) {
      int argIdx = freqs.length + i;
      Complex c = appliedCurve[i];
      curValue[i] = c.abs();
      curValue[argIdx] = atanc(c);
    }

    scaleValues(curValue, freqs, isLowFreq);

    return curValue;
  }

  /**
   * Function to run evaluation and backward difference for Jacobian approximation given a set of
   * points to set as response. Mainly a wrapper for the evaluateResponse function.
   *
   * @param variables Values to set the response's poles to
   * @param freqs Set of frequencies to get the response curve over
   * @param numZeros How many (paired) variables represent zeros (to determine first pole index)
   * @param fitResponse Response to apply these variables to
   * @param isLowFreq True if the calibration being fit to is low-frequency
   * @param correctionCurve Response correction used for some Trillium responses
   * @return RealVector with evaluation at current response value and RealMatrix with forward
   * difference approximation of that response's Jacobian
   */
  static Pair<RealVector, RealMatrix> jacobian(RealVector variables, double[] freqs,
      int numZeros, InstrumentResponse fitResponse, boolean isLowFreq,
      Complex[] correctionCurve) {
    int numVars = variables.getDimension();

    double[] currentVars = new double[numVars];

    for (int i = 0; i < numVars; ++i) {
      currentVars[i] = variables.getEntry(i);
    }

    double[] mag =
        evaluateResponse(currentVars, freqs, numZeros, fitResponse, isLowFreq, correctionCurve);

    double[][] jacobian = new double[mag.length][numVars];
    // now take the backward difference of each value
    for (int i = 0; i < numVars; ++i) {

      if (i % 2 == 1 && currentVars[i] == 0.) {
        // imaginary value already zero, don't change this
        // we assume that if an imaginary value is NOT zero, it's close enough
        // to its correct value that it won't get turned down to zero
        for (int j = 0; j < mag.length; ++j) {
          jacobian[j][i] = 0.;
        }
        continue;
      }

      double[] changedVars = Arrays.copyOf(currentVars, currentVars.length);

      // forward difference approximation -- ulp here gives us the error between this variable
      // and the smallest double larger than it. We multiply this by 100 to get the forward diff.
      // so that the difference is small relative to the value of the variable but also able to
      // give us a measurable change in the actual response curve function generated by it
      // and unlike having a fixed decimal step above the variable this is able to function on
      // floating-point numbers of arbitrary magnitude (useful for very high-freq poles in STS-6)
      double diffX = 100 * Math.ulp(changedVars[i]);
      changedVars[i] = changedVars[i] + diffX;

      double[] diffY =
          evaluateResponse(changedVars, freqs, numZeros, fitResponse, isLowFreq, correctionCurve);

      for (int j = 0; j < diffY.length; ++j) {
        jacobian[j][i] = diffY[j] - mag[j];
        jacobian[j][i] /= diffX;
      }

    }

    RealVector result = MatrixUtils.createRealVector(mag);
    RealMatrix jacobianMatrix = MatrixUtils.createRealMatrix(jacobian);

    return new Pair<>(result, jacobianMatrix);
  }

  /**
   * Given a candidate value for error terms for a given variable in the best-fit response, evaluate
   * the response given that value and estimate the jacobian by forward-difference.
   *
   * @param variables A vector with 2 entries representing a pole's real and complex value
   * respectively
   * @param freqs Frequencies to compute the response over
   * @param varIndex The entry in the fitResponse to be replaced with this given value
   * @param fitResponse The best-fit response being modified
   * @param pole True if the given variable being fit is a pole
   * @param isLowFreq True if the given calibration is low-frequency
   * @param freqToScaleAt Frequency that the source calculated response from cal is normalized on
   * @param correctionCurve Response correction used for some Trillium responses
   * @return Pair object holding the evaluation and the estimated jacobian at the given point
   * @see #jacobian(RealVector, double[], int, InstrumentResponse, boolean, Complex[])
   */
  static Pair<RealVector, RealMatrix> errorJacobian(RealVector variables, double[] freqs,
      int varIndex, InstrumentResponse fitResponse, boolean pole, boolean isLowFreq,
      double freqToScaleAt, Complex[] correctionCurve) {

    // variables should always be size 2 (fitting one pole or zero value at a time)
    Complex currentVar = new Complex(variables.getEntry(0), variables.getEntry(1));

    double[] mag = evaluateError(currentVar, freqs, varIndex, fitResponse, pole, isLowFreq,
        freqToScaleAt, correctionCurve);

    double[][] jacobian = new double[mag.length][2];

    // same procedure as forward difference in main jacobian function
    double diff = 100 * Math.ulp(currentVar.getReal());
    double next = diff + currentVar.getReal();
    Complex diffX = new Complex(next, currentVar.getImaginary());
    double[] diffY = evaluateError(diffX, freqs, varIndex, fitResponse, pole, isLowFreq,
        freqToScaleAt, correctionCurve);

    for (int j = 0; j < diffY.length; ++j) {
      jacobian[j][0] = diffY[j] - mag[j];
      jacobian[j][0] /= diff;
    }

    if (currentVar.getImaginary() != 0) {
      diff = 100 * Math.ulp(currentVar.getImaginary());
      next = diff + currentVar.getImaginary();
      diffX = new Complex(currentVar.getReal(), next);
      diffY = evaluateError(diffX, freqs, varIndex, fitResponse, pole, isLowFreq, freqToScaleAt,
          correctionCurve);

      for (int j = 0; j < diffY.length; ++j) {
        jacobian[j][1] = diffY[j] - mag[j];
        jacobian[j][1] /= diff;
      }
    } else {
      for (int j = 0; j < diffY.length; ++j) {
        jacobian[j][1] = 0.;
      }
    }

    RealVector result = MatrixUtils.createRealVector(mag);
    RealMatrix jacobianMatrix = MatrixUtils.createRealMatrix(jacobian);

    return new Pair<>(result, jacobianMatrix);
  }

  /**
   * Calculate the response of modified single pole-zero values in order to get an estimation of the
   * error for each term. Currently this is only performed on low-frequency calibrations for reasons
   * of performance.
   *
   * @param currentVar P/Z value to fit
   * @param freqs Range of frequencies to calculate the given response of
   * @param varIndex Index noting where in the best-fit response this value belongs
   * @param fitResponse Best-fit response returned by original solver
   * @param pole True if the value being fit is from the response's poles
   * @param isLowFreq True if the calibration being checked against is low-frequency
   * @param freqToScaleAt Frequency that the full calibration was normalized on
   * @param correctionCurve Response correction used for some Trillium responses
   * @return Response curve with the modified pole/zero value replaced in the first response
   */
  static double[] evaluateError(Complex currentVar, double[] freqs, int varIndex,
      InstrumentResponse fitResponse, boolean pole, boolean isLowFreq, double freqToScaleAt,
      Complex[] correctionCurve) {

    InstrumentResponse testResp = new InstrumentResponse(fitResponse);
    if (pole) {
      testResp.replaceFitPole(currentVar, varIndex, isLowFreq);
    } else {
      testResp.replaceFitZero(currentVar, varIndex, isLowFreq);
    }


    Complex[] appliedCurve = testResp.applyResponseToInputUnscaled(freqs);
    double[] curValue = new double[2 * freqs.length];

    if (correctionCurve != null) {
      for (int i = 0; i < appliedCurve.length; ++i) {
        appliedCurve[i] = appliedCurve[i].subtract(correctionCurve[i]);
      }
    }

    Complex[] normalizationData =
        testResp.applyResponseToInputUnscaled(new double[]{freqToScaleAt});
    assert (normalizationData.length == 1);
    Complex respAtNormalization = normalizationData[0];
    double normalizeAmp = 20 * Math.log10(respAtNormalization.abs());
    double normalizePhase = atanc(respAtNormalization);
    double phiPrev = atanc(appliedCurve[3 * appliedCurve.length / 4]);

    for (int i = 0; i < freqs.length; ++i) {
      int argIdx = i + freqs.length;
      Complex c = appliedCurve[i];
      curValue[i] = (20 * Math.log10(c.abs())) - normalizeAmp;
      double phi = atanc(c) - normalizePhase;
      phi = unwrap(phi, phiPrev);
      phiPrev = phi;
      curValue[argIdx] = Math.toDegrees(phi);
    }

    return curValue;
  }

  /**
   * Subtract a constant value from every point in the resp curve components such that the value at
   * a fixed given frequency is zero. This value is derived from the calibration type. If it is a
   * low-frequency cal, it is {@link #LOW_FREQ_ZERO_TARGET}, and if it is a high-frequency cal, then
   * it is {@link #HIGH_FREQ_ZERO_TARGET}. In the process the response amplitude is set to units of
   * dB, and the phase is set to units of degrees (instead of radians).
   *
   * @param unrot Unrotated response curve. First half is amplitude, second is phase
   * @param freqs Frequencies associated with the given response curve points
   * @param isLowFreq True if the calibration is low-frequency
   */
  static void scaleValues(double[] unrot, double[] freqs, boolean isLowFreq) {

    int normalIdx = FFTResult.getIndexOfFrequency(freqs,
        getFrequencyForNormalization(isLowFreq));
    // note that cal curve first half is the amplitude, second is phase
    int argStart = unrot.length / 2;
    // scale the first half of the data (amplitude)
    scaleMagnitude(unrot, freqs, isLowFreq);

    double unrotScaleArg = unrot[argStart + normalIdx];
    double phiPrev = 0;
    if (isLowFreq) {
      phiPrev = unrot[3 * unrot.length / 4];
    }
    for (int i = argStart; i < unrot.length; ++i) {
      double phi = unrot[i] - unrotScaleArg;
      phi = unwrap(phi, phiPrev);
      phiPrev = phi;
      unrot[i] = Math.toDegrees(phi);
    }
  }

  /**
   * Scale amplitude curve such that it is zero at the normalization frequency for the cal type.
   * This also puts the associated magnitude data into a dB-scale (20 * log10(amplitude))
   *
   * @param unscaled Unscaled data, first half of which is magnitude data
   * @param freqs Frequencies associated with each unmatched point
   * @param isLowFreq True if the calibration is low-frequency
   */
  static void scaleMagnitude(double[] unscaled, double[] freqs, boolean isLowFreq) {
    int normalIdx = FFTResult.getIndexOfFrequency(freqs,
        getFrequencyForNormalization(isLowFreq));
    double unrotScaleAmp = 20 * Math.log10(unscaled[normalIdx]);
    for (int i = 0; i < freqs.length; ++i) {
      double db = 20 * Math.log10(unscaled[i]);
      unscaled[i] = db - unrotScaleAmp;
    }
  }

  @Override
  public String[] getDataStrings() {

    List<Complex> fitPoles = getFitPoles();
    List<Complex> initialPoles = getInitialPoles();
    List<Complex> fitZeros = getFitZeros();
    List<Complex> initialZeros = getInitialZeros();

    if (fitPoles == null) {
      return new String[]{""};
    }

    double initialResidual = getInitResidual();
    double fitResidual = getFitResidual();

    StringBuilder sbInitialPoles = new StringBuilder();
    StringBuilder sbFitPoles = new StringBuilder();
    // add poles, initial then fit (single loop, append the two builders)
    sbInitialPoles.append("Initial poles: \n");
    sbFitPoles.append("Fit poles: \n");

    sbInitialPoles.append(complexListToString(initialPoles));
    if (poleErrors.size() > 0) {
      sbFitPoles.append(complexListToStringWithErrorTerms(fitPoles, poleErrors));
    } else {
      sbFitPoles.append(complexListToString(fitPoles));
    }

    sbInitialPoles.append("\n");
    sbFitPoles.append("\n");

    StringBuilder sbInitZ = new StringBuilder();
    StringBuilder sbFitZ = new StringBuilder();

    if (fitZeros.size() > 0) {
      sbInitZ.append("Initial zeros: \n");
      sbFitZ.append("Fit zeros: \n");
    }

    sbInitZ.append(complexListToString(initialZeros));
    if (zeroErrors.size() > 0) {
      sbFitZ.append(complexListToStringWithErrorTerms(fitZeros, zeroErrors));
    } else {
      sbFitZ.append(complexListToString(fitZeros));
    }

    sbFitPoles.append("\n");
    sbInitialPoles.append("\n");
    sbInitZ.append("\n");
    sbFitZ.append("\n");

    sbInitialPoles.append(sbFitPoles);

    sbInitZ.append(sbFitZ);

    String sbR = "Residuals:\n"
        + "Initial (nom. resp curve): "
        + DECIMAL_FORMAT.get().format(initialResidual)
        + "\nBest fit: "
        + DECIMAL_FORMAT.get().format(fitResidual);

    return new String[]{sbInitialPoles.toString(), sbInitZ.toString(), sbR};
  }

  public void setCorrectionResponse(SensorType typeToCorrect) throws IOException {
    if (typeToCorrect == null) {
      correctionToApply = null;
      return;
    }

    int index =
        Arrays.binarySearch(VALID_CORRECTIONS,1, VALID_CORRECTIONS.length, typeToCorrect);
    if (index <= 0) {
      correctionToApply = null;
      return;
    }
    String filename =
        getFilenameFromComponents(typeToCorrect, ResolutionType.STANDARD) + ".calresp";
    correctionToApply = loadEmbeddedResponse(filename);
  }

  /*
   * (non-Javadoc)
   * BACKEND FUNCTION BEGINS HERE
   * @see asl.sensor.Experiment#backend(asl.sensor.input.DataStore)
   */
  @Override
  protected void backend(DataStore dataStore) {
    numIterations = 0;

    DataBlock calib = dataStore.getBlock(0);
    DataBlock sensorOut = dataStore.getBlock(1);
    fitResponse = new InstrumentResponse(dataStore.getResponse(1));

    dataNames.add(calib.getName());
    String name = sensorOut.getName();
    dataNames.add(name);
    dataNames.add(fitResponse.getName());
    XYSeries calcMag = new XYSeries("Calc. resp. (" + name + ") magnitude");
    XYSeries calcArg = new XYSeries("Calc. resp. (" + name + ") phase");

    InstrumentResponse initResponse = new InstrumentResponse(fitResponse);
    initialPoles = new ArrayList<>(fitResponse.getPoles());
    initialZeros = new ArrayList<>(fitResponse.getZeros());

    // get the plots of the calculated response from deconvolution
    // PSD(out, in) / PSD(in, in) gives us PSD(out) / PSD(in) while removing
    // imaginary terms from the denominator due to multiplication with the
    // complex conjugate
    // PSD(out) / PSD(in) is the response curve (i.e., deconvolution)
    // also, use those frequencies to get the applied response to input
    fireStateChange("Getting PSDs of data...");
    FFTResult numeratorPSD, denominatorPSD, crossPSD;
    long interval = sensorOut.getInterval();
    {
      double taperWidth = DEFAULT_TAPER_WIDTH;
      // we should already have matching sample rates for data on experiment pre-processing steps
      double[] calData = calib.getData();
      double[] sensorData = sensorOut.getData();
      if (isLowFrequencyCalibration) {
        // we use a single-side taper of 25% because LF cals had consistent significant artifacts
        // around the corner that were consistent across instrument types and stations
        taperWidth = 0.25;
        // other attempts to fix those artifacts included downsampling and smoothing
        // but both had significant problems in replicating the shape of the curves
      }

      // now get the PSD data (already declared outside of this scope, so will persist)
      numeratorPSD = FFTResult.spectralCalc(sensorData, sensorData, interval, taperWidth);
      denominatorPSD = FFTResult.spectralCalc(calData, calData, interval, taperWidth);
      crossPSD = FFTResult.spectralCalc(sensorData, calData, interval, taperWidth);

    }
    double[] freqsUntrimmed = numeratorPSD.getFreqs(); // should be same for both results

    // store nyquist rate of data because freqs will be trimmed down later based on that value
    double nyquist = ONE_HZ_INTERVAL / (interval * 2.);

    // trim frequency window in order to restrict range of response fits
    int startIndex, endIndex, maxPlotIndex;
    // extfreq is how far out to extend data past range of fit
    // maxFreq is the largest frequency that we use in the solver
    // low frequency cal fits over a different range of data
    if (isLowFrequencyCalibration) {
      double minFreq = 0.001; // 1000s period
      double maxFreq = 0.05; // 20s period
      //double maxPlotFreq = maxFreq;
      startIndex = FFTResult.getIndexOfFrequency(freqsUntrimmed, minFreq) + 1;
      endIndex = FFTResult.getIndexOfFrequency(freqsUntrimmed, maxFreq);
      maxPlotIndex = endIndex;
    } else {
      double minFreq = .2; // lower bound of .2 Hz (5s period) due to noise
      // get factor of nyquist rate, again due to noise
      double maxFreq = nyquistMultiplier * nyquist;
      double maxPlotFreq = PEAK_MULTIPLIER * nyquist; // always plot up to 90% of nyquist
      startIndex = FFTResult.getIndexOfFrequency(freqsUntrimmed, minFreq);
      endIndex = FFTResult.getIndexOfFrequency(freqsUntrimmed, maxFreq);
      maxPlotIndex = FFTResult.getIndexOfFrequency(freqsUntrimmed, maxPlotFreq);
      // maxFreq = extFreq;
    }

    fireStateChange("Finding and trimming data to relevant frequency range");

    // now trim frequencies to in range
    // start and end are the region of the fit area, extIdx is the full plotted region
    double[] plottingFreqs = Arrays.copyOfRange(freqsUntrimmed, startIndex, maxPlotIndex);
    freqs = Arrays.copyOfRange(freqsUntrimmed, startIndex, endIndex);

    // index to fix the curve values to 0; both amplitude and phase set to 0 at that frequency
    int normalIdx = FFTResult.getIndexOfFrequency(freqs,
        getFrequencyForNormalization(isLowFrequencyCalibration));

    // trim the PSDs to the data in the trimmed frequency range
    Complex[] numeratorPSDVals = numeratorPSD.getFFT();
    Complex[] denominatorPSDVals = denominatorPSD.getFFT();
    Complex[] crossPSDVals = crossPSD.getFFT();
    double[] untrimmedAmplitude = new double[freqsUntrimmed.length];
    double[] untrimmedPhase = new double[freqsUntrimmed.length];

    // calculated response from deconvolving calibration from signal
    // (this will be in displacement and need to be integrated)
    for (int i = 0; i < freqsUntrimmed.length; ++i) {
      Complex ampNumer = numeratorPSDVals[i];
      Complex phaseNumer = crossPSDVals[i];
      double denom = denominatorPSDVals[i].abs(); // phase is 0
      // the actual complex value which we'll immediately convert to doubles for use in plots/fits
      Complex ampValue = ampNumer.divide(denom);
      Complex scaleFactor = getSignalScalingFactor(freqsUntrimmed[i]);
      // convert from displacement to velocity
      ampValue = ampValue.multiply(scaleFactor.pow(2));
      untrimmedAmplitude[i] = 10 * Math.log10(ampValue.abs());
      Complex phaseValue = phaseNumer.divide(denom);
      untrimmedPhase[i] = atanc(phaseValue.multiply(scaleFactor));
    }

    int offset = 0;

    if (!isLowFrequencyCalibration) {
      fireStateChange("Smoothing calculated resp data...");
      // smoothingPoints should be an even number to prevent biasing on the half-sample
      // since the averaging isn't done from the center of a sample but from either the left or right
      offset = 3; // how far to shift the data after the smoothing has been done
      int smoothingPoints = 2 * offset;
      // now smooth the data
      // scan starting at the high-frequency range for low-freq data (will be trimmed) & vice-versa
      untrimmedAmplitude = multipointMovingAverage(untrimmedAmplitude,
          smoothingPoints, !isLowFrequencyCalibration);
      // phase smoothing also includes an unwrapping step
      untrimmedPhase = unwrapArray(untrimmedPhase);
      untrimmedPhase = multipointMovingAverage(untrimmedPhase, smoothingPoints,
          !isLowFrequencyCalibration);
    }

    fireStateChange("Trimming calculated resp data down to range of interest...");
    double[] plottedAmp = Arrays.copyOfRange(untrimmedAmplitude, startIndex, maxPlotIndex);
    double[] plottedPhs = Arrays.copyOfRange(untrimmedPhase, startIndex, maxPlotIndex);

    fireStateChange("Scaling & weighting calculated resp data in preparation for solving");
    // get the data at the normalized index, use this to scale the data
    double freqToScaleAt = freqs[normalIdx];
    double ampScale = plottedAmp[normalIdx];
    double phsScale = plottedPhs[normalIdx];
    double[] observedResult = new double[2 * freqs.length];
    double[] weights = new double[observedResult.length];
    maxArgWeight = Double.MIN_VALUE;
    maxMagWeight = Double.MIN_VALUE;
    for (int i = 0; i < freqs.length; ++i) {
      double xAxis = freqs[i];
      if (!plotUsingHz) {
        xAxis = 1. / freqs[i];
      }

      // scale
      double amp = plottedAmp[i];
      double phs = plottedPhs[i];
      amp -= ampScale;
      phs -= phsScale;
      phs = Math.toDegrees(phs);
      // get weight
      maxMagWeight = Math.max(maxMagWeight, Math.abs(amp));
      maxArgWeight = Math.max(maxArgWeight, Math.abs(phs));

      // add to plot (re-wrap phase), add to fitting array
      int argIdx = i + freqs.length;
      observedResult[i] = amp;
      observedResult[argIdx] = phs;
      calcMag.add(xAxis, amp);
      phs = rewrapAngleDegrees(phs);
      calcArg.add(xAxis, phs);
    }

    maxMagWeight = 1. / maxMagWeight;
    if (maxArgWeight != 0) {
      maxArgWeight = 1. / maxArgWeight;
    }

    // apply weights
    Arrays.fill(weights, 0, freqs.length, maxMagWeight);
    Arrays.fill(weights, freqs.length, weights.length, maxArgWeight);

    // get the rest of the plotted data squared away
    for (int i = freqs.length; i < plottingFreqs.length; ++i) {
      double xAxis = plottingFreqs[i];
      if (!plotUsingHz) {
        xAxis = 1. / plottingFreqs[i];
      }
      // scale
      double amp = plottedAmp[i];
      double phs = plottedPhs[i];
      amp -= ampScale;
      phs -= phsScale;
      phs = Math.toDegrees(phs);
      // add to plot (re-wrap phase)
      calcMag.add(xAxis, amp);
      phs = rewrapAngleDegrees(phs);
      calcArg.add(xAxis, phs);
    }

    DiagonalMatrix weightMat = new DiagonalMatrix(weights);

    fireStateChange("Getting estimate and setting up solver...");

    // now to set up a solver for the params -- first, get the input variables
    // complex values are technically two variables, each a double
    // so, let's split up their real and im components and make each one a
    // variable. (we also need to ignore conjugate values, for constraints)
    RealVector initialGuess, initialPoleGuess, initialZeroGuess;

    // The peak value here used to be set to nyquistMultiplier * nyquist, so ONLY points
    // that are within the frequency band being fit could be modified.
    // This has since been changed specifically because some responses like the STS-2.5
    // have curves that are very dependent on poles which are out-of-band, and so we will
    // not present an upper bound on the frequency range for HF cals in order to better
    // produce fits for cases of such data.
    initialPoleGuess = fitResponse
        .polesToVector(isLowFrequencyCalibration, Double.MAX_VALUE);
    initialZeroGuess = fitResponse
        .zerosToVector(isLowFrequencyCalibration, Double.MAX_VALUE);
    int numZeros = initialZeroGuess.getDimension();
    initialGuess = initialZeroGuess.append(initialPoleGuess);

    // now, solve for the response that gets us the best-fit response curve
    // RealVector initialGuess = MatrixUtils.createRealVector(responseVariables);
    RealVector obsResVector = MatrixUtils.createRealVector(observedResult);

    // we don't need to get the response correction here because the normalization is done
    // when the full response curve is separated into its amplitude and phase
    Complex[] responseCorrection = getResponseCorrection(freqs, correctionToApply);
    MultivariateJacobianFunction jacobian = new MultivariateJacobianFunction() {

      final double[] freqsSet = freqs;
      final int numZerosSet = numZeros;
      final boolean isLowFrequency = isLowFrequencyCalibration;
      final InstrumentResponse fitSet = fitResponse;

      @Override
      public Pair<RealVector, RealMatrix> value(final RealVector point) {
        ++numIterations;
        fireStateChange("Fitting, iteration count " + numIterations);
        return jacobian(point, freqsSet, numZerosSet, fitSet, isLowFrequency, responseCorrection);
      }

    };

    // We define these variable separately in the event they may need to be tweaked for odd cases,
    // though evidence suggests that in both high AND low frequency cals that 1E-10 is a suitable
    // value to apply to for both parameters and produces good fit curves.
    // These are relative tolerance precisions that are already lower in magnitude than the last
    // digits of precision of most resp files' poles and zeros, at least.
    final double costTolerance = 1.0E-10;
    final double paramTolerance = 1.0E-10;

    LeastSquaresOptimizer optimizer = new LevenbergMarquardtOptimizer().
        withCostRelativeTolerance(costTolerance).
        withOrthoTolerance(1E-25). // arbitrarily low value, which should ideally never be relevant
        withParameterRelativeTolerance(paramTolerance);

    // set up structures that will hold the initial and final response plots
    name = fitResponse.getName();
    XYSeries initMag = new XYSeries("Initial param (" + name + ") magnitude");
    XYSeries initArg = new XYSeries("Initial param (" + name + ") phase");

    XYSeries fitMag = new XYSeries("Fit resp. magnitude");
    XYSeries fitArg = new XYSeries("Fit resp. phase");

    LeastSquaresProblem lsp = new LeastSquaresBuilder().
        start(initialGuess).
        target(obsResVector).
        model(jacobian).
        weight(weightMat).
        parameterValidator(new PoleValidator(numZeros)).
        lazyEvaluation(false).
        maxEvaluations(Integer.MAX_VALUE).
        maxIterations(Integer.MAX_VALUE).
        build();

    fireStateChange("Built least-squares problem; evaluating intial guess...");

    // residuals used to determine quality of solution convergence

    LeastSquaresProblem.Evaluation initEval = lsp.evaluate(initialGuess);
    initialResidual = initEval.getCost();

    fireStateChange("Got initial evaluation; running solver...");

    RealVector finalResultVector;

    LeastSquaresOptimizer.Optimum optimum = optimizer.optimize(lsp);
    finalResultVector = optimum.getPoint();
    numIterations = optimum.getIterations();

    LeastSquaresProblem.Evaluation evaluation = lsp.evaluate(finalResultVector);
    fitResidual = evaluation.getCost();
    double[] fitParams = evaluation.getPoint().toArray();
    // get results from evaluating the function at the two points

    XYSeries initResidMag = new XYSeries("Percent error of init. amplitude");
    XYSeries initResidPhase = new XYSeries("Percent error of with init. phase");
    XYSeries fitResidMag = new XYSeries("Percent error of fit amplitude");
    XYSeries fitResidPhase = new XYSeries("Percent error of with fit phase");

    fitResponse = fitResponse.buildResponseFromFitVector(
        fitParams, isLowFrequencyCalibration, numZeros);
    fitPoles = fitResponse.getPoles();
    fitZeros = fitResponse.getZeros();

    // response error term calculation here (3-sigma bounds)
    poleErrors = new HashMap<>();
    zeroErrors = new HashMap<>();

    if (isLowFrequencyCalibration) {
      constructErrorTerms(observedResult, numZeros, fitParams, freqToScaleAt);
    }

    fireStateChange("Getting extended resp curves for high-freq plots...");
    // we use the apply response method here to get the full range of plotted data, not just fit
    Complex[] correction = getResponseCorrection(plottingFreqs, correctionToApply);
    Complex[] init = initResponse.applyResponseToInputUnscaled(plottingFreqs);
    Complex[] fit = fitResponse.applyResponseToInputUnscaled(plottingFreqs);
    double[] initialValues = new double[plottingFreqs.length * 2];
    double[] fitValues = new double[plottingFreqs.length * 2];
    for (int i = 0; i < plottingFreqs.length; ++i) {
      Complex initEntry = init[i].subtract(correction[i]);
      Complex fitEntry = fit[i].subtract(correction[i]);
      int argIdx = plottingFreqs.length + i;
      initialValues[i] = initEntry.abs();
      initialValues[argIdx] = atanc(initEntry);
      fitValues[i] = fitEntry.abs();
      fitValues[argIdx] = atanc(fitEntry);
    }
    fireStateChange("Scaling extended resps...");
    scaleValues(initialValues, plottingFreqs, isLowFrequencyCalibration);
    scaleValues(fitValues, plottingFreqs, isLowFrequencyCalibration);

    fireStateChange("Compiling data into plots...");

    for (int i = 0; i < plottingFreqs.length; ++i) {
      double xValue;
      if (plotUsingHz) {
        xValue = plottingFreqs[i];
      } else {
        xValue = 1. / plottingFreqs[i];
      }

      int argIdx = initialValues.length / 2 + i;
      initMag.add(xValue, initialValues[i]);
      initArg.add(xValue, rewrapAngleDegrees(initialValues[argIdx]));
      fitMag.add(xValue, fitValues[i]);
      fitArg.add(xValue, rewrapAngleDegrees(fitValues[argIdx]));

      if (i < freqs.length) {
        int obsArgIdx = i + freqs.length; // observedResult cuts off before freqsFull does
        double initAmpNumer = Math.pow(10, initialValues[i] / 20);
        double fitAmpNumer = Math.pow(10, fitValues[i] / 20);

        double obsAmpDbl = observedResult[i];
        if (obsAmpDbl == 0.) {
          obsAmpDbl = Double.MIN_VALUE;
        }

        obsAmpDbl = Math.pow(10, obsAmpDbl / 20);

        double errInitMag = 100. * (initAmpNumer - obsAmpDbl) / obsAmpDbl;
        double errFitMag = 100. * (fitAmpNumer - obsAmpDbl) / obsAmpDbl;
        if (!Double.isInfinite(errInitMag)) {
          initResidMag.add(xValue, Math.abs(errInitMag));

        }
        if (!Double.isInfinite(errFitMag)) {
          fitResidMag.add(xValue, Math.abs(errFitMag));
        }

        double obsPhase = observedResult[obsArgIdx];
        if (obsPhase != 0.) {
          double errInitPhase = Math.abs(100 * (initialValues[argIdx] - obsPhase) / obsPhase);
          double errFitPhase = Math.abs(100 * (fitValues[argIdx] - obsPhase) / obsPhase);
          initResidPhase.add(xValue, errInitPhase);
          fitResidPhase.add(xValue, errFitPhase);
        }
      }
    }

    XYSeriesCollection xysc = new XYSeriesCollection();
    xysc.addSeries(initMag);
    xysc.addSeries(calcMag);
    xysc.addSeries(fitMag);
    xySeriesData.add(xysc);

    xysc = new XYSeriesCollection();
    xysc.addSeries(initArg);
    xysc.addSeries(calcArg);
    xysc.addSeries(fitArg);
    xySeriesData.add(xysc);

    xysc = new XYSeriesCollection();
    xysc.addSeries(initResidMag);
    xysc.addSeries(fitResidMag);
    xySeriesData.add(xysc);

    xysc = new XYSeriesCollection();
    xysc.addSeries(initResidPhase);
    xysc.addSeries(fitResidPhase);
    xySeriesData.add(xysc);
  }

  private void constructErrorTerms(double[] observedResult, int numZeros, double[] fitParams,
      double freqToScaleAt) {
    int currentZeroIndex = 0; // where zero under analysis lies in the response
    int currentPoleIndex = 0; // as above for pole

    LeastSquaresOptimizer optimizer = new LevenbergMarquardtOptimizer().
        withCostRelativeTolerance(1E-5).
        withOrthoTolerance(1E-25).
        withParameterRelativeTolerance(1E-5);
    for (int i = 0; i < fitParams.length; i += 2) {
      boolean pole = i >= numZeros;
      Complex fitTerm = new Complex(fitParams[i], fitParams[i + 1]);
      double corner = fitTerm.abs() / TAU;
      // get the frequency range over the octave centered on the p/z corner frequency
      int lowIndex = FFTResult.getIndexOfFrequency(freqs, corner / Math.sqrt(2.));
      int highIndex = FFTResult.getIndexOfFrequency(freqs, Math.sqrt(2.) * corner);

      if (lowIndex == highIndex) {
        // in this case the pole/zero is either well beyond the minimum plotted frequency value
        // or exists well inside the flat band of the sensor, so we will not fit it
        continue;
      }

      double[] errorTermFreqsFull = Arrays.copyOfRange(freqs, lowIndex, highIndex);
      double[] trimmedMagAndPhaseCurve = new double[errorTermFreqsFull.length * 2];
      double[] weights = new double[trimmedMagAndPhaseCurve.length - 2];
      Arrays.fill(weights, 0, weights.length / 2, maxMagWeight);
      Arrays.fill(weights, weights.length / 2, weights.length, maxArgWeight);
      // copy the relevant portion of the magnitude curve
      System.arraycopy(observedResult, lowIndex, trimmedMagAndPhaseCurve, 0,
          highIndex - lowIndex);
      // and copy the relevant portion of the phase curve (second half is where phase starts)
      System.arraycopy(observedResult, observedResult.length / 2 + lowIndex,
          trimmedMagAndPhaseCurve, highIndex - lowIndex, highIndex - lowIndex);

      final int index;
      // we keep track of count so that we can have a 1:1 mapping between error terms and
      // listed p/z values in the

      if (pole) {
        index = currentPoleIndex++; // return before incrementing
        // increment again to skip over complex conjugate for nonzero imaginary terms
        if (fitTerm.getImaginary() != 0) {
          ++currentPoleIndex;
        }
      } else {
        index = currentZeroIndex++; // again, return before incrementing
        if (fitTerm.getImaginary() != 0) {
          ++currentZeroIndex;
        }
      }

      List<Complex> bestFits = new ArrayList<>();

      for (int j = 0; j < errorTermFreqsFull.length; ++j) {
        String message = "Estimating error for variable" + (i / 2 + 1) + " of " +
            fitParams.length / 2 + " using frequency range " + (j + 1) + " of " +
            errorTermFreqsFull.length;
        fireStateChange(message);

        // get all but one frequency (and corresponding magnitude) term
        final double[] errorTermFreqs = new double[errorTermFreqsFull.length - 1];
        System.arraycopy(errorTermFreqsFull, 0, errorTermFreqs, 0, j);
        if (j + 1 < errorTermFreqsFull.length) {
          System.arraycopy(errorTermFreqsFull, j + 1,
              errorTermFreqs, j, errorTermFreqs.length - j);
        }
        // subtract 2 to reflect removal of specific pole and specific zero
        final double[] observedCurve = new double[trimmedMagAndPhaseCurve.length - 2];
        int trimOffset = trimmedMagAndPhaseCurve.length / 2; // where the phase component starts
        // copy the first j points to this array (from 0 to j-1)
        System.arraycopy(trimmedMagAndPhaseCurve, 0, observedCurve, 0, j);
        // now copy from from j+1 to where the phase component starts
        System.arraycopy(trimmedMagAndPhaseCurve, j + 1, observedCurve, j,
            errorTermFreqs.length - j);
        // now the first j phase components (note offset for destination due to missing phase term)
        System.arraycopy(trimmedMagAndPhaseCurve, trimOffset, observedCurve,
            trimOffset - 1, j);
        System.arraycopy(trimmedMagAndPhaseCurve, trimOffset + j + 1, observedCurve,
            trimOffset - 1 + j, errorTermFreqs.length - j);

        Complex[] correctedCurve = getResponseCorrection(errorTermFreqs, correctionToApply);
        // we DO need to do normalization here because the normalization requires calculating
        // the curve at a specific value that otherwise does not have this correction applied
        // and this will (hopefully) save us having to do scaling each iteration
        Complex normalizationForCorrection =
            getResponseCorrection(new double[]{freqToScaleAt}, correctionToApply)[0];
        for (int w = 0; w < correctedCurve.length; ++w) {
          correctedCurve[w] = correctedCurve[w].subtract(normalizationForCorrection);
        }

        MultivariateJacobianFunction errorJacobian = new MultivariateJacobianFunction() {
          final double[] freqsSet = errorTermFreqs;
          final int variableIndex = index;
          final boolean isLowFrequency = isLowFrequencyCalibration;
          final InstrumentResponse fitSet = fitResponse;

          @Override
          public Pair<RealVector, RealMatrix> value(final RealVector point) {
            ++numIterations;
            fireStateChange("Fitting, iteration count " + numIterations);
            return errorJacobian(point, freqsSet, variableIndex, fitSet, isLowFrequency, pole,
                freqToScaleAt, correctedCurve);
          }
        };

        RealVector initialError = MatrixUtils.createRealVector(
            new double[]{fitParams[i], fitParams[i + 1]});
        RealVector observed = MatrixUtils.createRealVector(observedCurve);

        LeastSquaresProblem errorLsq = new LeastSquaresBuilder().
            start(initialError).
            target(observed).
            model(errorJacobian).
            parameterValidator(new PoleValidator(numZeros)).
            weight(new DiagonalMatrix(weights)).
            lazyEvaluation(false).
            maxEvaluations(Integer.MAX_VALUE).
            maxIterations(Integer.MAX_VALUE).
            build();

        LeastSquaresOptimizer.Optimum optimum = optimizer.optimize(errorLsq);
        RealVector errorVector = optimum.getPoint();

        Complex c = new Complex(errorVector.getEntry(0), errorVector.getEntry(1));
        bestFits.add(c);
      } // end loop over frequency range (error term estimation for a given point)

      // now that we have a list of deviations for this pole or zero,
      // it's time to do statistics to it to get our variation
      Complex threeSigma =
          getComplexSDev(bestFits.toArray(new Complex[]{})).multiply(3);

      if (pole) {
        poleErrors.put(fitTerm, threeSigma);
        if (fitTerm.getImaginary() != 0) {
          poleErrors.put(fitTerm.conjugate(), threeSigma);
        }
      } else {
        zeroErrors.put(fitTerm, threeSigma);
        if (fitTerm.getImaginary() != 0) {
          zeroErrors.put(fitTerm.conjugate(), threeSigma);
        }
      }
    } // end loop over the fit poles and zeros
  }

  @Override
  public int blocksNeeded() {
    return 2;
  }

  /**
   * Calculates the scaling/rotation factor to be done on the calibration signal. If the calibration
   * is capacitive this value is always 1.
   *
   * @param freq Given frequency to apply the calculation to
   * @return Either 1 or 2*pi*i*freq depending on calibration type
   */
  private Complex getSignalScalingFactor(double freq) {
    if (isCapacitive) {
      return new Complex(1.);
    }
    return new Complex(0., TAU * freq);
  }

  /**
   * Get the poles that the solver has found to best-fit the est. response
   *
   * @return new poles that should improve fit over inputted response, as a list
   */
  public List<Complex> getFitPoles() {
    List<Complex> polesOut = new ArrayList<>();
    Set<Complex> retain = new HashSet<>(fitPoles);
    retain.removeAll(initialPoles);
    for (Complex c : fitPoles) {
      if (retain.contains(c)) {
        polesOut.add(c);
      }
    }
    complexRealsFirstSorter(polesOut);
    return polesOut;
  }

  /**
   * Get the residual value from the solved response parameters
   *
   * @return the residual of the solved-for poles (best-fit response)
   */
  public double getFitResidual() {
    return fitResidual;
  }

  /**
   * Provides a handle to the fit response object for parameter output
   *
   * @return the best-fit response
   */
  public InstrumentResponse getFitResponse() {
    return fitResponse;
  }

  /**
   * Gets the correction being applied to the data, which can be null if not one listed in
   * {@link #VALID_CORRECTIONS}
   * @return Response defining correction to zeros for some Trillium sensors
   */
  InstrumentResponse getCorrectionToApply() {
    return correctionToApply;
  }

  /**
   * Get the zeros fitted from the experiment
   *
   * @return List of zeros (complex numbers) that are used in best-fit curve
   */
  public List<Complex> getFitZeros() {
    List<Complex> zerosOut = new ArrayList<>();
    Set<Complex> retain = new HashSet<>(fitZeros);
    retain.removeAll(initialZeros);
    for (Complex c : fitZeros) {
      if (retain.contains(c)) {
        zerosOut.add(c);
      }
    }
    complexRealsFirstSorter(zerosOut);
    return zerosOut;
  }


  /**
   * Get the error terms of all the fitted poles. A map is used so as to allow efficient access for
   * an error term by lookup from a value gotten in getFitPoles, which will be sorted into an order
   * that may not match the original response.
   *
   * @return Map of complex pole values to the respective real/imaginary error terms, each as a
   * single complex value
   */
  public Map<Complex, Complex> getPoleErrors() {
    return poleErrors;
  }

  /**
   * Get the error terms of all the fitted zeros. A map is used so as to allow efficient access for
   * an error term by lookup from a value gotten in getFitZeros, which will be sorted into an order
   * that may not match the original response.
   *
   * @return Map of complex zero values to the respective real/imaginary error terms, each as a
   * single complex value
   */
  public Map<Complex, Complex> getZeroErrors() {
    return zeroErrors;
  }

  /**
   * Get poles used in input response, for reference against best-fit poles
   *
   * @return poles taken from initial response file
   */
  public List<Complex> getInitialPoles() {
    List<Complex> polesOut = new ArrayList<>();
    Set<Complex> retain = new HashSet<>(initialPoles);
    retain.removeAll(fitPoles);
    for (Complex c : initialPoles) {
      if (retain.contains(c)) {
        polesOut.add(c);
      }
    }
    complexRealsFirstSorter(polesOut);
    return polesOut;
  }

  /**
   * Get initial zeros from (nominal) response file used as input
   *
   * @return zeros taken from initial response file
   */
  public List<Complex> getInitialZeros() {
    List<Complex> zerosOut = new ArrayList<>();
    Set<Complex> retain = new HashSet<>(initialZeros);
    retain.removeAll(fitZeros);
    for (Complex c : initialZeros) {
      if (retain.contains(c)) {
        zerosOut.add(c);
      }
    }
    complexRealsFirstSorter(zerosOut);
    return zerosOut;
  }

  /**
   * Get the residual value of the initial response parameters
   *
   * @return the residual of the initial poles from fed-in response
   */
  public double getInitResidual() {
    return initialResidual;
  }

  /**
   * Get the number of times the algorithm iterated to produce the optimum response fit, from the
   * underlying least squares solver
   *
   * @return the number of iterations
   */
  public int getIterations() {
    return numIterations;
  }

  /**
   * Get the highest frequency value included in data set for the solver. This is used to mark where
   * on the chart to plot a vertical line to denote the solver bounds
   *
   * @return Max frequency over fit range for HF cals
   */
  public double getMaxFitFrequency() {
    return freqs[freqs.length - 1];
  }

  /**
   * Get the values used to weight the residual calculation function. The first value is the
   * magnitude weighting, the second is phase.
   *
   * @return Weighting values for least-squared error terms of
   */
  public double[] getWeights() {
    return new double[]{maxMagWeight, maxArgWeight};
  }

  @Override
  public boolean hasEnoughData(DataStore dataStore) {
    return (dataStore.blockIsSet(0) && dataStore.bothComponentsSet(1));
  }

  @Override
  public int[] listActiveResponseIndices() {
    // NOTE: not used by corresponding panel, overrides with active indices
    // of components in the combo-box
    return new int[]{1};
  }

  /**
   * Used to control the scaling of the calibration signals based on whether or not a given
   * calibration is capacitive. The default value is false (resistive).
   *
   * @param isCapacitive True if the calibration to be calculated is capacitive
   */
  public void setCapactiveCalibration(boolean isCapacitive) {
    this.isCapacitive = isCapacitive;
  }

  /**
   * Set the new peak multiplier for the data region under analysis. This should be a positive
   * value, and is bounded by 0.8 (@see NumericUtils.PEAK_MULTIPLIER) This may be useful when trying
   * to run the solver over a high-frequency calibration where the data is particularly noisy in
   * some of the higher-frequency bounds, causing a bad corner fit to occur. We also enforce this to
   * be over 0.3 as well.
   *
   * @param newMultiplier New maximum fraction of nyquist rate to fit data over (should be from 0.3
   * to 0.8).
   */
  public void setNyquistMultiplier(double newMultiplier) {
    nyquistMultiplier = Math.min(newMultiplier, PEAK_MULTIPLIER);
    nyquistMultiplier = Math.max(MIN_MULTIPLIER, nyquistMultiplier);
  }

  /**
   * Auto-determines if a calibration is a long-period calibration or not based on the length of the
   * data being inputted. Because most high-frequency cals are around 15 minutes and most
   * low-frequency cals are several hours, we use a cutoff of one hour to make this determination.
   *
   * @param ds Data store which testing is to be done on
   */
  public void autoDetermineCalibrationStatus(DataStore ds) {
    if (!ds.blockIsSet(0)) {
      return;
    }
    long start = ds.getBlock(0).getStartTime();
    long end = ds.getBlock(0).getEndTime();
    long timeDiff = end - start;
    // this will be true if the calibration has a duration of more than one hour (3.6E6 ms)
    // i.e., most high-frequency calibrations will be
    isLowFrequencyCalibration = timeDiff >= 3.6E6;
  }

  /**
   * Get status if poles/zeros being fit are in low (<1Hz) or high (>1Hz) frequency band
   *
   * @return true if cal being run is a low-frequency cal
   */
  public boolean isLowFrequencyCalibration() {
    return isLowFrequencyCalibration;
  }

  /**
   * Determines which poles to fit when doing the response curve fitting; low frequency calibrations
   * set the first two poles; high frequency calibrations set the remaining poles
   *
   * @param isLowFrequencyCalibration True if a low frequency calibration is to be used
   */
  public void setLowFrequencyCalibration(boolean isLowFrequencyCalibration) {
    this.isLowFrequencyCalibration = isLowFrequencyCalibration;
  }

  /**
   * Set whether or not to plot in units of frequency (Hz) or period (s)
   *
   * @param setFreq true if plots should be in frequency units (Hz)
   */
  public void setPlotUsingHz(boolean setFreq) {
    plotUsingHz = setFreq;
  }


  private static class PoleValidator implements ParameterValidator {

    int numZeros;

    public PoleValidator(int numZeros) {
      this.numZeros = numZeros;
    }


    /**
     * Simple validator method to enforce poles to be negative for their values (Since imaginary
     * values are stored in resps as their value and complex conjugate, we can mandate this for all
     * values in the vector, though only real components are strictly required to be negative).
     *
     * @param poleParams RealVector of parameters to be evaluated by solver
     * @return Vector of parameters but with components all negative
     */
    @Override
    public RealVector validate(RealVector poleParams) {
      for (int i = 0; i < poleParams.getDimension(); ++i) {
        double value = poleParams.getEntry(i);
        if (value > 0 && (i % 2) == 0) {
          // even index means this is a real-value vector entry
          // if it's above zero, put it back below zero
          // with delta offset to prevent it from going to zero on the iterative step.
          poleParams.setEntry(i, Math.nextDown(-Double.MIN_VALUE));
        } else if (value > 0) {
          // this means the value is complex, we can multiply it by -1
          // this is ok for complex values since their conjugate is implied
          // to be part of the set of poles being fit
          poleParams.setEntry(i, value * -1);
        }
      }
      return poleParams;
    }
  }

}