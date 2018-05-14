package asl.sensor.experiment;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
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
import asl.sensor.input.DataBlock;
import asl.sensor.input.DataStore;
import asl.sensor.input.InstrumentResponse;
import asl.sensor.utils.FFTResult;
import asl.sensor.utils.NumericUtils;
import asl.sensor.utils.TimeSeriesUtils;

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
public class RandomizedExperiment extends Experiment implements ParameterValidator {

  static final double DELTA = 1E-12;
  private static final double PEAK_MULTIPLIER = InstrumentResponse.PEAK_MULTIPLIER;

  private double initialResidual, fitResidual;
  private List<Complex> initialPoles;
  private List<Complex> fitPoles;
  private List<Complex> initialZeros;
  private List<Complex> fitZeros;

  /**
   * True if calibration is low frequency.
   * This affects which poles are fitted, either low or high frequencies.
   */
  private boolean isLowFrequencyCalibration;

  private InstrumentResponse fitResponse;

  private double[] freqs;

  private boolean plotUsingHz;

  private double maxMagWeight, maxArgWeight; // max values of magnitude, phase
  private double nyquistMultiplier; // region up to nyquist to take for data

  private static final double ZERO_TARGET = 0.02; // location of value to set to 0 in curves for scaling
  private int numZeros; // how many entries in parameter vector define zeros
  private int numIterations; // how much the solver ran

  public RandomizedExperiment() {
    super();
    isLowFrequencyCalibration = false;
    numIterations = 0;
    plotUsingHz = true;
    nyquistMultiplier = PEAK_MULTIPLIER;
  }

  /*
   * (non-Javadoc)
   * BACKEND FUNCTION BEGINS HERE
   * @see asl.sensor.experiment.Experiment#backend(asl.sensor.input.DataStore)
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
    double nyquist = TimeSeriesUtils.ONE_HZ_INTERVAL / sensorOut.getInterval();
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
      maxPlotFreq = InstrumentResponse.PEAK_MULTIPLIER * nyquist; // i.e., 80% of nyquist
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

    double zeroTarget = 0.02; // frequency to set all curves to zero at
    int normalIdx = FFTResult.getIndexOfFrequency(freqs, zeroTarget);

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
      Complex scaleFactor = new Complex(0., NumericUtils.TAU * freqsUntrimmed[i]);
      // convert from displacement to velocity
      ampValue = ampValue.multiply(scaleFactor.pow(2));
      untrimmedAmplitude[i] = 10 * Math.log10(ampValue.abs());
      Complex phaseValue = phaseNumer.divide(denom);
      untrimmedPhase[i] = NumericUtils.atanc(phaseValue.multiply(scaleFactor));
    }

    fireStateChange("Smoothing calculated resp data...");
    // smoothingPoints should be an even number to prevent biasing on the half-sample
    // since the averaging isn't done from the center of a sample but from either the left or right
    int offset = 3; // how far to shift the data after the smoothing has been done
    int smoothingPoints = 2 * offset;
    // now smooth the data
    // scan starting at the high-frequency range for low-freq data (will be trimmed) & vice-versa
    untrimmedAmplitude = NumericUtils.multipointMovingAverage(untrimmedAmplitude,
        smoothingPoints, !isLowFrequencyCalibration);
    // phase smoothing also includes an unwrapping step
    untrimmedPhase = NumericUtils.unwrapArray(untrimmedPhase);
    untrimmedPhase = NumericUtils.multipointMovingAverage(untrimmedPhase, smoothingPoints,
        !isLowFrequencyCalibration);

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
    maxArgWeight = 1;
    maxMagWeight = 0;
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

    maxMagWeight = 1000 / maxMagWeight;
    if (maxArgWeight != 0) {
      maxArgWeight = 1. / maxArgWeight;
    }

    // apply weights
    for (int i = 0; i < freqs.length; ++i) {
      int argIndex = i + freqs.length;
      double denominator;
      if (!isLowFrequencyCalibration) {
        if (freqs[i] < 1) {
          denominator = 1; // weight everything up to 1Hz equally
        } else {
          denominator = freqs[i]; // set everything (else) to 1/f weighting
        }
      } else {
        if (freqs[i] < .01) {
          denominator = freqs[i];
        } else {
          denominator = .01;
        }
      }
      weights[argIndex] = maxArgWeight / denominator;
      weights[i] = maxMagWeight / denominator;
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
      // add to plot (re-wrap phase)
      calcMag.add(xAxis, amp);
      calcArg.add(xAxis, NumericUtils.rewrapAngleDegrees(phs) );
    }


    DiagonalMatrix weightMat = new DiagonalMatrix(weights);

    fireStateChange("Getting estimate and setting up solver...");

    // now to set up a solver for the params -- first, get the input variables
    // complex values are technically two variables, each a double
    // so, let's split up their real and im components and make each one a
    // variable. (we also need to ignore conjugate values, for constraints)
    RealVector initialGuess, initialPoleGuess, initialZeroGuess;

    initialPoleGuess = fitResponse.polesToVector(isLowFrequencyCalibration, nyquistMultiplier * nyquist);
    initialZeroGuess = fitResponse.zerosToVector(isLowFrequencyCalibration, nyquistMultiplier * nyquist);
    numZeros = initialZeroGuess.getDimension();
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

    double costTolerance = 1.0E-15;
    double paramTolerance = 1.0E-10;
    // probably acceptable tolerance for clean low-frequency cals BUT
    // high frequency cals are noisy and slow to converge
    // so this branch is to enable using higher tolerance to deal with that
    if (!isLowFrequencyCalibration) {
      costTolerance = 1.0E-15;
      paramTolerance = 1.0E-10;
    }

    LeastSquaresOptimizer optimizer = new LevenbergMarquardtOptimizer().
        withCostRelativeTolerance(costTolerance).
        withOrthoTolerance(1E-25).
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
        parameterValidator(this).
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
    XYSeries initResidPhase = new XYSeries("Diff. with init phase");
    XYSeries fitResidMag = new XYSeries("Percent error of fit amplitude");
    XYSeries fitResidPhase = new XYSeries("Diff with fit phase");

    fitResponse = fitResponse.buildResponseFromFitVector(
        fitParams, isLowFrequencyCalibration, numZeros);
    fitPoles = fitResponse.getPoles();
    fitZeros = fitResponse.getZeros();

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
    scaleValues(initialValues, freqs, isLowFrequencyCalibration);
    scaleValues(fitValues, freqs, isLowFrequencyCalibration);

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

        initResidPhase.add(xValue, Math.abs(initialValues[argIdx] - observedResult[obsArgIdx]));
        fitResidPhase.add(xValue, Math.abs(fitValues[argIdx] - observedResult[obsArgIdx]));
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

  @Override
  public int blocksNeeded() {
    return 2;
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

      double[] changedVars = new double[currentVars.length];
      System.arraycopy(currentVars, 0, changedVars, 0, currentVars.length);

      double diffX = changedVars[i] - DELTA;
      changedVars[i] = diffX;
      double[] diffY =
          evaluateResponse(changedVars, freqs, numZeros, fitResponse, isLowFreq);

      for (int j = 0; j < diffY.length; ++j) {
        if (changedVars[i] - currentVars[i] == 0.) {
          jacobian[j][i] = 0.;
        } else {
          jacobian[j][i] = mag[j] - diffY[j];
          jacobian[j][i] /= currentVars[i] - changedVars[i];
        }
      }

    }

    RealVector result = MatrixUtils.createRealVector(mag);
    RealMatrix jacobianMatrix = MatrixUtils.createRealMatrix(jacobian);

    return new Pair<>(result, jacobianMatrix);
  }

  @Override
  public int[] listActiveResponseIndices() {
    // NOTE: not used by corresponding panel, overrides with active indices
    // of components in the combo-box
    return new int[]{1};
  }

  /**
   * Set the new peak multiplier for the data region under analysis.
   * This should be a positive value, and is bounded by 0.8 (@see NumericUtils.PEAK_MULTIPLIER)
   * This may be useful when trying to run the solver over a high-frequency calibration where
   * the data is particularly noisy in some of the higher-frequency bounds, causing a bad corner
   * fit to occur. We also enforce this to be over 0.3 as well.
   * @param newMultiplier New maximum fraction of nyquist rate to fit data over (should be
   * from 0.3 to 0.8).
   */
  public void setNyquistMultiplier(double newMultiplier) {
    nyquistMultiplier = Math.min(newMultiplier, PEAK_MULTIPLIER);
    nyquistMultiplier = Math.max(0.3, nyquistMultiplier);
  }

  static void scaleValues(double[] unrot, double[] freqs, boolean isLowFrequencyCalibration) {
    int normalIdx = FFTResult.getIndexOfFrequency(freqs, ZERO_TARGET);
    int argStart = unrot.length / 2;
    double unrotScaleAmp = 20 * Math.log10(unrot[normalIdx]);
    double unrotScaleArg = unrot[argStart + normalIdx];
    double phiPrev = 0;
    if (isLowFrequencyCalibration) {
      phiPrev = unrot[3 * unrot.length / 4];
    }
    for (int i = 0; i < argStart; ++i) {
      int argIdx = argStart + i;
      double db = 20 * Math.log10(unrot[i]);
      unrot[i] = db - unrotScaleAmp;
      double phi = unrot[argIdx] - unrotScaleArg;
      phi = NumericUtils.unwrap(phi, phiPrev);
      phiPrev = phi;
      unrot[argIdx] = Math.toDegrees(phi);
    }
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
   * Set whether or not to plot in units of frequency (Hz) or period (s)
   *
   * @param setFreq true if plots should be in frequency units (Hz)
   */
  public void setPlotUsingHz(boolean setFreq) {
    plotUsingHz = setFreq;
  }

  /**
   * Simple validator method to enforce poles to be negative for their values
   * (Since imaginary values are stored in resps as their value and complex
   * conjugate, we can mandate this for all values in the vector, though only
   * real components are strictly required to be negative).
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
        poleParams.setEntry(i, -DELTA - Double.MIN_VALUE);
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
