package asl.sensor.experiment;

import asl.sensor.input.DataStore;
import asl.utils.FFTResult;
import asl.utils.NumericUtils;
import asl.utils.TimeSeriesUtils;
import asl.utils.input.DataBlock;
import asl.utils.input.InstrumentResponse;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.complex.ComplexFormat;
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
 * This experiment takes in a randomized calibration signal and the
 * corresponding output from a seismic sensor. It calculates the implicit
 * response by deconvolving the calibration signal from the sensor output, and
 * then finds the best-fit poles (lowest 2 for low-frequency calibrations,
 * all remaining poles for high-frequency calibrations) to
 * match the magnitude and rotation angle of the calculated response curve
 * produced from the deconvolution.
 * Plottable data includes the sensor's response curve
 * (Bode plot), the calculated response from deconvolution, and the plot
 * of the response from the best-fit parameters.
 * These plots are returned in a list of two datasets: the first holds
 * the absolute values per-frequency of those curves, and the second holds
 * the angle of each such point in complex space.
 * For more details on the algorithm, see Ringler, Hutt, et al.,
 * "Estimating Pole-Zero Errors in GSN-IRIS/USGS Network Calibration Metadata",
 * Bulletin of the Seismological Society of America, Vol 102 (Apr. 2012).
 *
 * @author akearns - KBRWyle
 */
public class RandomizedExperiment extends Experiment {
  
  /**
   * Maximum possible frequency value as a multiple of nyquist (0.9).
   * The solver will still default to 0.8 as results above that are very unstable for noisy cals
   */
  public static final double PEAK_MULTIPLIER = 0.9;

  /**
   * Sets the default normalization point for curves (0.02 Hz)
   */
  private static final double ZERO_TARGET = 0.02;
  private double initialResidual, fitResidual;
  private List<Complex> initialPoles;
  private List<Complex> fitPoles;
  private Map<Complex, Complex> poleErrors;
  private List<Complex> initialZeros;
  private List<Complex> fitZeros;
  private Map<Complex, Complex> zeroErrors;

  /**
   * 5-place precision data to match the level of precision in a typical response file
   * This decimal formatter is used for formatting error term real and imaginary parts
   * as well as the basis for the complex formatter for returned pole/zero values
   */
  public static final ThreadLocal<DecimalFormat> HIGH_PRECISION_DF =
      ThreadLocal.withInitial(() -> new DecimalFormat("#.#####"));


  /**
   * 5-place precision formatter for pole/zero values
   */
  public static final ThreadLocal<ComplexFormat> COMPLEX_FORMAT =
      // we want to use a 5-place precision for our data if it's at all possible
      ThreadLocal.withInitial(() -> new ComplexFormat(HIGH_PRECISION_DF.get()));



  /**
   * True if calibration is low frequency.
   * This affects which poles are fitted, either low or high frequencies.
   */
  private boolean isLowFrequencyCalibration;
  /**
   * Used to ensure that scaling is done correctly on calibrations
   * that use a capacitive setting. Default value will be false.
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
    nyquistMultiplier = 0.8; // defaults to 0.8
  }

  private static String complexListToString(List<Complex> complexList) {
    final int MAX_LINE = 2; // maximum number of entries per line

    ComplexFormat complexFormat = new ComplexFormat(DECIMAL_FORMAT.get());
    StringBuilder stringBuilder = new StringBuilder();
    int numInLine = 0;

    for (Complex number : complexList) {

      double initPrd = NumericUtils.TAU / number.abs();

      stringBuilder.append(complexFormat.format(number));
      // p/z value as period isn't as important to get so precisely so we use the standard formatter
      stringBuilder.append(" (");
      stringBuilder.append(DECIMAL_FORMAT.get().format(initPrd));
      stringBuilder.append(")");
      ++numInLine;
      // want to fit two to a line for paired values
      if (numInLine >= MAX_LINE) {
        stringBuilder.append("\n");
        numInLine = 0;
      } else {
        stringBuilder.append(", ");
      }
    }

    return stringBuilder.toString();
  }

  private static String complexListToStringWithErrorTerms(
      List<Complex> complexList, Map<Complex, Complex> errorTerms) {
    final int MAX_LINE = 2; // maximum number of entries per line

    StringBuilder stringBuilder = new StringBuilder();
    int numInLine = 0;

    for (Complex number : complexList) {

      double initPrd = NumericUtils.TAU / number.abs();


      stringBuilder.append(COMPLEX_FORMAT.get().format(number));
      stringBuilder.append(" (");
      stringBuilder.append(HIGH_PRECISION_DF.get().format(initPrd));
      stringBuilder.append(")");

      if (errorTerms.containsKey(number)) {
        Complex error = errorTerms.get(number);
        stringBuilder.append(" [Err: +/- ");
        stringBuilder.append(HIGH_PRECISION_DF.get().format(error.getReal()));
        stringBuilder.append(", ");
        stringBuilder.append(HIGH_PRECISION_DF.get().format(error.getImaginary()));
        stringBuilder.append("i]");
      }

      ++numInLine;
      // want to fit two to a line for paired values
      if (numInLine >= MAX_LINE) {
        stringBuilder.append("\n");
        numInLine = 0;
      } else {
        stringBuilder.append(", ");
      }
    }

    return stringBuilder.toString();
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

  /**
   * Backend function to set instrument response according to current
   * test variables (for best-fit calculation / backward difference) and
   * produce a response from that result. The passed response is copied on
   * start and is not modified directly. Which values (poles) are modified
   * depends on high or low frequency calibration setting.
   *
   * @param variables values to set the instrument response to
   * @return Doubles representing new response curve evaluation
   */
  private static double[] evaluateResponse(double[] variables, double[] freqs, int numZeros,
      InstrumentResponse fitResponse, boolean isLowFrequencyCalibration) {

    InstrumentResponse testResp = new InstrumentResponse(fitResponse);

    // prevent terrible case where, say, only high-freq poles above nyquist rate
    if (variables.length > 0) {
      testResp = fitResponse.buildResponseFromFitVector(
          variables, isLowFrequencyCalibration, numZeros);
    } else {
      System.out.println("NO VARIABLES TO SET. THIS IS AN ERROR.");
    }

    Complex[] appliedCurve = testResp.applyResponseToInput(freqs);
    double[] curValue = new double[freqs.length * 2];

    for (int i = 0; i < freqs.length; ++i) {
      int argIdx = freqs.length + i;
      Complex c = appliedCurve[i];
      curValue[i] = c.abs();
      curValue[argIdx] = NumericUtils.atanc(c);
    }

    scaleValues(curValue, freqs, isLowFrequencyCalibration);

    return curValue;
  }

  /**
   * Function to run evaluation and backward difference for Jacobian
   * approximation given a set of points to set as response.
   * Mainly a wrapper for the evaluateResponse function.
   *
   * @param variables Values to set the response's poles to
   * @return RealVector with evaluation at current response value and
   * RealMatrix with backward difference of that response (Jacobian)
   */
  static Pair<RealVector, RealMatrix> jacobian(RealVector variables, double[] freqs,
      int numZeros, InstrumentResponse fitResponse, boolean isLowFreq) {
    int numVars = variables.getDimension();

    double[] currentVars = new double[numVars];

    for (int i = 0; i < numVars; ++i) {
      currentVars[i] = variables.getEntry(i);
    }

    double[] mag = evaluateResponse(currentVars, freqs, numZeros, fitResponse, isLowFreq);

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


      double diffX = 100 * Math.ulp(changedVars[i]);
      changedVars[i] = changedVars[i] + diffX;


      double[] diffY =
          evaluateResponse(changedVars, freqs, numZeros, fitResponse, isLowFreq);

      for (int j = 0; j < diffY.length; ++j) {
        jacobian[j][i] = diffY[j] - mag[j];
        jacobian[j][i] /= diffX;
      }

    }

    RealVector result = MatrixUtils.createRealVector(mag);
    RealMatrix jacobianMatrix = MatrixUtils.createRealMatrix(jacobian);

    return new Pair<>(result, jacobianMatrix);
  }

  static Pair<RealVector, RealMatrix> errorJacobian(RealVector variables, double[] freqs,
      int varIndex, InstrumentResponse fitResponse, boolean pole, boolean isLowFreq) {

    // variables should always be size 2 (fitting one pole or zero value at a time)
    Complex currentVar = new Complex(variables.getEntry(0), variables.getEntry(1));

    double[] mag = evaluateError(currentVar, freqs, varIndex, fitResponse, pole, isLowFreq);

    double[][] jacobian = new double[mag.length][2];

    double diff = 100 * Math.ulp(currentVar.getReal());
    double next = diff + currentVar.getReal();
    Complex diffX = new Complex(next, currentVar.getImaginary());
    double[] diffY = evaluateError(diffX, freqs, varIndex, fitResponse, pole, isLowFreq);

    for (int j = 0; j < diffY.length; ++j) {
        jacobian[j][0] = diffY[j] - mag[j];
        jacobian[j][0] /= diff;
    }

    if (currentVar.getImaginary() != 0) {
      diff = 100 * Math.ulp(currentVar.getImaginary());
      next = diff + currentVar.getImaginary();
      diffX = new Complex(currentVar.getReal(), next);
      diffY = evaluateError(diffX, freqs, varIndex, fitResponse, pole, isLowFreq);

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

  static double[] evaluateError(Complex currentVar, double[] freqs, int varIndex,
      InstrumentResponse fitResponse, boolean pole, boolean isLowFreq) {

    InstrumentResponse testResp = new InstrumentResponse(fitResponse);
    if (pole) {
      testResp.replaceFitPole(currentVar, varIndex, isLowFreq);
    } else {
      testResp.replaceFitZero(currentVar, varIndex, isLowFreq);
    }

    Complex[] appliedCurve = testResp.applyResponseToInput(freqs);
    double[] curValue = new double[freqs.length];

    for (int i = 0; i < freqs.length; ++i) {
      Complex c = appliedCurve[i];
      curValue[i] = c.abs();
    }

    scaleMagnitude(curValue, freqs);

    return curValue;
  }

  static void scaleValues(double[] unrot, double[] freqs, boolean isLowFrequencyCalibration) {

    int normalIdx = FFTResult.getIndexOfFrequency(freqs, ZERO_TARGET);
    int argStart = unrot.length / 2;

    scaleMagnitude(unrot, freqs);

    double unrotScaleArg = unrot[argStart + normalIdx];
    double phiPrev = 0;
    if (isLowFrequencyCalibration) {
      phiPrev = unrot[3 * unrot.length / 4];
    }
    for (int i = argStart; i < unrot.length; ++i) {
      double phi = unrot[i] - unrotScaleArg;
      phi = NumericUtils.unwrap(phi, phiPrev);
      phiPrev = phi;
      unrot[i] = Math.toDegrees(phi);
    }
  }


  static void scaleMagnitude(double[] unscaled, double[] freqs){
    int normalIdx = FFTResult.getIndexOfFrequency(freqs, ZERO_TARGET);
    double unrotScaleAmp = 20 * Math.log10(unscaled[normalIdx]);
    for (int i = 0; i < freqs.length; ++i) {
      double db = 20 * Math.log10(unscaled[i]);
      unscaled[i] = db - unrotScaleAmp;
    }
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
    numeratorPSD = FFTResult.spectralCalc(sensorOut, sensorOut);
    denominatorPSD = FFTResult.spectralCalc(calib, calib);
    crossPSD = FFTResult.spectralCalc(sensorOut, calib);

    double[] freqsUntrimmed = numeratorPSD.getFreqs(); // should be same for both results

    // store nyquist rate of data because freqs will be trimmed down later
    double nyquist = TimeSeriesUtils.ONE_HZ_INTERVAL / (double) sensorOut.getInterval();
    nyquist = nyquist / 2.;

    // trim frequency window in order to restrict range of response fits
    double minFreq, maxFreq, maxPlotFreq;
    // extfreq is how far out to extend data past range of fit
    // maxFreq is the largest frequency that we use in the solver
    // low frequency cal fits over a different range of data
    if (isLowFrequencyCalibration) {
      minFreq = 0.001; // 1000s period
      maxFreq = 0.05; // 20s period
      maxPlotFreq = maxFreq;
    } else {
      minFreq = .2; // lower bound of .2 Hz (5s period) due to noise
      // get factor of nyquist rate, again due to noise
      maxFreq = nyquistMultiplier * nyquist;
      maxPlotFreq = PEAK_MULTIPLIER * nyquist; // always plot up to 90% of nyquist
      // maxFreq = extFreq;
    }

    fireStateChange("Finding and trimming data to relevant frequency range");
    // now trim frequencies to in range
    // start and end are the region of the fit area, extIdx is the full plotted region
    int startIndex = FFTResult.getIndexOfFrequency(freqsUntrimmed, minFreq);
    int endIndex = FFTResult.getIndexOfFrequency(freqsUntrimmed, maxFreq);

    //Used for plotting full range
    int maxPlotIndex = FFTResult.getIndexOfFrequency(freqsUntrimmed, maxPlotFreq);

    double[] plottingFreqs = Arrays.copyOfRange(freqsUntrimmed, startIndex, maxPlotIndex);
    freqs = Arrays.copyOfRange(freqsUntrimmed, startIndex, endIndex);

    int normalIdx = FFTResult.getIndexOfFrequency(freqs, ZERO_TARGET);

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
      untrimmedPhase[i] = NumericUtils.atanc(phaseValue.multiply(scaleFactor));
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
      untrimmedAmplitude = NumericUtils.multipointMovingAverage(untrimmedAmplitude,
          smoothingPoints, !isLowFrequencyCalibration);
      // phase smoothing also includes an unwrapping step
      untrimmedPhase = NumericUtils.unwrapArray(untrimmedPhase);
      untrimmedPhase = NumericUtils.multipointMovingAverage(untrimmedPhase, smoothingPoints,
          !isLowFrequencyCalibration);
    }

    // experimentation with offsets to deal with the way the moving average shifts the data
    // since the plot is basically logarithmic this only matters due to the limited data
    // on the low-frequency corner -- shifting here by half the smoothing re-centers the data
    if (isLowFrequencyCalibration) {
      startIndex -= offset;
      maxPlotIndex -= offset;
    }

    fireStateChange("Trimming calculated resp data down to range of interest...");
    double[] plottedAmp = Arrays.copyOfRange(untrimmedAmplitude, startIndex, maxPlotIndex);
    double[] plottedPhs = Arrays.copyOfRange(untrimmedPhase, startIndex, maxPlotIndex);

    fireStateChange("Scaling & weighting calculated resp data in preparation for solving");
    // get the data at the normalized index, use this to scale the data
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
      calcArg.add(xAxis, NumericUtils.rewrapAngleDegrees(phs));
    }

    maxMagWeight = 1. / maxMagWeight;
    if (maxArgWeight != 0) {
      maxArgWeight = 1. / maxArgWeight;
    }

    // apply weights
    for (int i = 0; i < freqs.length; ++i) {
      int argIndex = i + freqs.length;
      weights[argIndex] = maxArgWeight;
      weights[i] = maxMagWeight; // / denominator;
    }

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
      calcArg.add(xAxis, NumericUtils.rewrapAngleDegrees(phs));
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

    MultivariateJacobianFunction jacobian = new MultivariateJacobianFunction() {

      final double[] freqsSet = freqs;
      final int numZerosSet = numZeros;
      final boolean isLowFrequency = isLowFrequencyCalibration;
      final InstrumentResponse fitSet = fitResponse;

      @Override
      public Pair<RealVector, RealMatrix> value(final RealVector point) {
        ++numIterations;
        fireStateChange("Fitting, iteration count " + numIterations);
        return jacobian(point, freqsSet, numZerosSet, fitSet, isLowFrequency);
      }

    };

    final double costTolerance = 1.0E-10;
    // probably acceptable tolerance for clean low-frequency cals BUT
    // high frequency cals are noisy and slow to converge
    // that said using these parameters in high-freq cases still seems to work

    LeastSquaresOptimizer optimizer = new LevenbergMarquardtOptimizer().
        withCostRelativeTolerance(costTolerance).
        withOrthoTolerance(1E-25).
        withParameterRelativeTolerance(1E-10);

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
      constructErrorTerms(observedResult, numZeros, fitParams);
    }

    fireStateChange("Getting extended resp curves for high-freq plots...");
    // we use the apply response method here to get the full range of plotted data, not just fit
    Complex[] init = initResponse.applyResponseToInput(plottingFreqs);
    Complex[] fit = fitResponse.applyResponseToInput(plottingFreqs);
    double[] initialValues = new double[plottingFreqs.length * 2];
    double[] fitValues = new double[plottingFreqs.length * 2];
    for (int i = 0; i < plottingFreqs.length; ++i) {
      int argIdx = plottingFreqs.length + i;
      initialValues[i] = init[i].abs();
      initialValues[argIdx] = NumericUtils.atanc(init[i]);
      fitValues[i] = fit[i].abs();
      fitValues[argIdx] = NumericUtils.atanc(fit[i]);
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
      double plotArg;
      initMag.add(xValue, initialValues[i]);
      plotArg = NumericUtils.rewrapAngleDegrees(initialValues[argIdx]);
      initArg.add(xValue, plotArg);
      fitMag.add(xValue, fitValues[i]);
      plotArg = NumericUtils.rewrapAngleDegrees(fitValues[argIdx]);
      fitArg.add(xValue, plotArg);

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

  private void constructErrorTerms(double[] observedResult, int numZeros, double[] fitParams) {
    int currentZeroIndex = 0; // where zero under analysis lies in the response
    int currentPoleIndex = 0; // as above for pole

    LeastSquaresOptimizer optimizer;
    LeastSquaresOptimizer.Optimum optimum; // Get error analysis for each pole/zero value

    for (int i = 0; i < fitParams.length; i += 2) {
      boolean pole = i >= numZeros;
      Complex fitTerm = new Complex(fitParams[i], fitParams[i+1]);
      double corner = fitTerm.abs() / NumericUtils.TAU;
      // get the frequency range over the octave centered on the p/z corner frequency
      int lowIndex = FFTResult.getIndexOfFrequency(freqs, corner / Math.sqrt(2.));
      int highIndex = FFTResult.getIndexOfFrequency(freqs, Math.sqrt(2.) * corner);

      double[] errorTermFreqsFull = Arrays.copyOfRange(freqs, lowIndex, highIndex);
      double[] observedMagnitudeFull = Arrays.copyOfRange(observedResult, lowIndex, highIndex);
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
      optimizer = new LevenbergMarquardtOptimizer().
          withCostRelativeTolerance(1E-5).
          withOrthoTolerance(1E-25).
          withParameterRelativeTolerance(1E-5);

      for (int j = 0; j < errorTermFreqsFull.length; ++j) {
        String message = "Estimating error for variable" + (i/2 + 1) +  " of " +
            fitParams.length/2 + " using frequency range " + (j + 1) + " of " +
            errorTermFreqsFull.length;
        fireStateChange(message);

        // get all but one frequency (and corresponding magnitude) term
        final double[] errorTermFreqs = new double[errorTermFreqsFull.length - 1];
        System.arraycopy(errorTermFreqsFull, 0, errorTermFreqs, 0, j);
        if (j + 1 < errorTermFreqsFull.length) {
          System.arraycopy(errorTermFreqsFull, j + 1,
              errorTermFreqs, j, errorTermFreqs.length - j);
        }
        final double[] observedMagnitude = new double[errorTermFreqsFull.length - 1];
        System.arraycopy(observedMagnitudeFull, 0, observedMagnitude, 0, j);
        if (j + 1 < observedMagnitudeFull.length) {
          System.arraycopy(observedMagnitudeFull, j + 1,
              observedMagnitude, j, observedMagnitude.length - j);
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
            return errorJacobian(point, freqsSet, variableIndex, fitSet, isLowFrequency, pole);
          }
        };

        RealVector initialError = MatrixUtils.createRealVector(
            new double[]{fitParams[i], fitParams[i+1]});
        RealVector observed = MatrixUtils.createRealVector(observedMagnitude);

        LeastSquaresProblem errorLsq = new LeastSquaresBuilder().
            start(initialError).
            target(observed).
            model(errorJacobian).
            parameterValidator(new PoleValidator(numZeros)).
            lazyEvaluation(false).
            maxEvaluations(Integer.MAX_VALUE).
            maxIterations(Integer.MAX_VALUE).
            build();

        optimum = optimizer.optimize(errorLsq);
        RealVector errorVector = optimum.getPoint();
        Complex c = new Complex(errorVector.getEntry(0), errorVector.getEntry(1));
        bestFits.add(c);
      } // end loop over frequency range (error term estimation for a given point)

      // now that we have a list of best-fit p/z over range, we get the standard deviation
      Complex threeSigma =
          NumericUtils.getComplexSDev(bestFits.toArray(new Complex[]{})).multiply(3);

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
    }
  }

  @Override
  public int blocksNeeded() {
    return 2;
  }

  /**
   * Calculates the scaling/rotation factor to be done on the calibration
   * signal. If the calibration is capacitive this value is always 1.
   * @param freq Given frequency to apply the calculation to
   * @return Either 1 or 2*pi*i*freq depending on calibration type
   */
  private Complex getSignalScalingFactor(double freq) {
    if (isCapacitive) {
      return new Complex(1.);
    }
    return new Complex(0., NumericUtils.TAU * freq);
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
    NumericUtils.complexRealsFirstSorter(polesOut);
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
    NumericUtils.complexRealsFirstSorter(zerosOut);
    return zerosOut;
  }


  /**
   * Get the error terms of all the fitted poles.
   * A map is used so as to allow efficient access for an error term by lookup from a value gotten
   * in getFitPoles, which will be sorted into an order that may not match the original response.
   * @return Map of complex pole values to the respective real/imaginary error terms, each as a
   * single complex value
   */
  public Map<Complex, Complex> getPoleErrors() {
    return poleErrors;
  }

  /**
   * Get the error terms of all the fitted zeros.
   * A map is used so as to allow efficient access for an error term by lookup from a value gotten
   * in getFitZeros, which will be sorted into an order that may not match the original response.
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
    NumericUtils.complexRealsFirstSorter(polesOut);
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
    NumericUtils.complexRealsFirstSorter(zerosOut);
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
   * Get the number of times the algorithm iterated to produce the optimum
   * response fit, from the underlying least squares solver
   *
   * @return the number of iterations
   */
  public int getIterations() {
    return numIterations;
  }

  public double getMaxFitFrequency() {
    return freqs[freqs.length - 1];
  }

  /**
   * Get the values used to weight the residual calculation function.
   * The first value is the magnitude weighting, the second is phase.
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
   * Used to control the scaling of the calibration signals based on whether or
   * not a given calibration is capacitive. The default value is false (resistive).
   * @param isCapacitive True if the calibration to be calculated is capacitive
   */
  public void setCapactiveCalibration(boolean isCapacitive) {
    this.isCapacitive = isCapacitive;
  }

  /**
   * Set the new peak multiplier for the data region under analysis.
   * This should be a positive value, and is bounded by 0.8 (@see NumericUtils.PEAK_MULTIPLIER)
   * This may be useful when trying to run the solver over a high-frequency calibration where
   * the data is particularly noisy in some of the higher-frequency bounds, causing a bad corner
   * fit to occur. We also enforce this to be over 0.3 as well.
   *
   * @param newMultiplier New maximum fraction of nyquist rate to fit data over (should be
   * from 0.3 to 0.8).
   */
  public void setNyquistMultiplier(double newMultiplier) {
    nyquistMultiplier = Math.min(newMultiplier, PEAK_MULTIPLIER);
    nyquistMultiplier = Math.max(0.3, nyquistMultiplier);
  }

  /**
   * Determines which poles to fit when doing the response curve fitting;
   * low frequency calibrations set the first two poles; high frequency
   * calibrations set the remaining poles
   *
   * @param isLowFrequencyCalibration True if a low frequency calibration is to be used
   */
  public void setLowFrequencyCalibration(boolean isLowFrequencyCalibration) {
    this.isLowFrequencyCalibration = isLowFrequencyCalibration;
  }

  /**
   * Auto-determines if a calibration is a long-period calibration or not based on the length of
   * the data being inputted. Because most high-frequency cals are around 15 minutes and most
   * low-frequency cals are several hours, we use a cutoff of one hour to make this determination.
   *
   * @param ds
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
   * @return true if cal being run is a low-frequency cal
   */
  public boolean isLowFrequencyCalibration() {
    return isLowFrequencyCalibration;
  }

  /**
   * Set whether or not to plot in units of frequency (Hz) or period (s)
   *
   * @param setFreq true if plots should be in frequency units (Hz)
   */
  public void setPlotUsingHz(boolean setFreq) {
    plotUsingHz = setFreq;
  }


  private class PoleValidator implements ParameterValidator {

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
