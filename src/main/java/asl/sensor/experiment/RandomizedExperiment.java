package asl.sensor.experiment;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.fitting.leastsquares.EvaluationRmsChecker;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresBuilder;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresOptimizer;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresProblem;
import 
org.apache.commons.math3.fitting.leastsquares.LevenbergMarquardtOptimizer;
import 
org.apache.commons.math3.fitting.leastsquares.MultivariateJacobianFunction;
import org.apache.commons.math3.fitting.leastsquares.ParameterValidator;
import org.apache.commons.math3.linear.DiagonalMatrix;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.optim.ConvergenceChecker;
import org.apache.commons.math3.util.Pair;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import asl.sensor.input.DataBlock;
import asl.sensor.input.DataStore;
import asl.sensor.input.InstrumentResponse;
import asl.sensor.utils.FFTResult;
import asl.sensor.utils.FFTResult.TaperType;
import asl.sensor.utils.NumericUtils;

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
 * @author akearns
 *
 */
public class RandomizedExperiment 
extends Experiment implements ParameterValidator {

  private static final double DELTA = 1E-7;
  public static final double PEAK_MULTIPLIER = 
      NumericUtils.PEAK_MULTIPLIER; // max pole-fit frequency
  
  // To whomever has to maintain this code after I'm gone:
  // I'm sorry, I'm so so sorry
  // I suppose it's a little neater now that some functions are part of the
  // response class? It's still inherently nasty due to issues relating to
  // converting complex lists into arrays of doubles in order to use the solver
  
  private double initialResidual, fitResidual;
  private List<Complex> initialPoles;
  private List<Complex> fitPoles;
  private List<Complex> initialZeros;
  private List<Complex> fitZeros;
  
  // when true, doesn't run solver, in event parameters have an issue
  // (does the solver seem to have frozen? try rebuilding with this as true,
  // and then run the plot -- show nominal resp. and estimated curves)
  public final boolean SKIP_SOLVING = false;
  
  private boolean lowFreq; // fit the low- or high-frequency poles?
  
  private InstrumentResponse fitResponse;
  
  private double[] freqs;
  private double nyquist;
  
  private boolean freqSpace;
  
  private double maxMagWeight, maxArgWeight; // max values of magnitude, phase
  
  private int normalIdx; // location of value to set to 0 in curves for scaling
  private int numZeros; // how many entries in parameter vector define zeros
  private int sensorOutIdx; // location to load response from?
  private int numIterations; // how much the solver ran
  
  public RandomizedExperiment() {
    super();
    lowFreq = false;
    normalIdx = 0;
    numIterations = 0;
    freqSpace = true;
  }
  
  /**
   * Set whether or not to plot in units of frequency (Hz) or period (s)
   * @param setFreq true if plots should be in frequency units (Hz)
   */
  public void useFreqUnits(boolean setFreq) {
    freqSpace = setFreq;
  }
  
  /*
   * (non-Javadoc)
   * BACKEND FUNCTION BEGINS HERE
   * @see asl.sensor.experiment.Experiment#backend(asl.sensor.input.DataStore)
   */
  @Override
  protected void backend(DataStore ds) {
    
    // construct response plot
    DataBlock calib = ds.getXthLoadedBlock(1);
    sensorOutIdx = ds.getXthFullyLoadedIndex(1);
    
    // if first data has response loaded erroneously, load in next data set
    if (sensorOutIdx == 0) {
      sensorOutIdx = ds.getXthFullyLoadedIndex(2);
    }
    
    /*
    if ( ds.getBlock(sensorOutIdx).getName().equals( calib.getName() ) ) {
      sensorOutIdx = ds.getXthFullyLoadedIndex(2);
    }
    */

    DataBlock sensorOut = ds.getBlock(sensorOutIdx);
    fitResponse = new InstrumentResponse( ds.getResponse(sensorOutIdx) );
    
    dataNames.add( calib.getName() );
    dataNames.add( sensorOut.getName() );
    dataNames.add( fitResponse.getName() );
    
    initialPoles = new ArrayList<Complex>( fitResponse.getPoles() );
    initialZeros = new ArrayList<Complex>( fitResponse.getZeros() );
    
    // get the plots of the calculated response from deconvolution
    // PSD(out, in) / PSD(in, in) gives us PSD(out) / PSD(in) while removing
    // imaginary terms from the denominator due to multiplication with the
    // complex conjugate
    // PSD(out) / PSD(in) is the response curve (i.e., deconvolution)
    int windowSize, change;
    TaperType taper;
    // also, use those frequencies to get the applied response to input
    if (lowFreq) {
      int maxLen = Math.max( sensorOut.size(), calib.size() );
      windowSize = 2;
      while (windowSize <= maxLen) {
        windowSize *= 2;
      }
      windowSize *= 2;
      change = windowSize;
      taper = TaperType.MULT;
    } else {
      windowSize = sensorOut.size() / 4;
      change = windowSize / 4;
      taper = TaperType.COS;
    }
    
    FFTResult numeratorPSD = 
        FFTResult.spectralCalc(sensorOut, calib, windowSize, change, taper);
    FFTResult denominatorPSD = 
        FFTResult.spectralCalc(calib, calib, windowSize, change, taper);
    freqs = numeratorPSD.getFreqs(); // should be same for both results
    
    // store nyquist rate of data because freqs will be trimmed down later
    nyquist = sensorOut.getSampleRate() / 2.;
    nyquist += .5; // increasing to prevent issues with pole frequency rounding 
    
    // trim frequency window in order to restrict range of response fits
    double minFreq, maxFreq;
    
    // low frequency cal fits over a different range of data
    if (lowFreq) {
      minFreq = 0.001; // 1000s period
      maxFreq = 0.05; // 20s period
    } else {
      minFreq = .2; // lower bound of .2 Hz (5s period) due to noise
      // get up to .8 of nyquist rate, again due to noise
      maxFreq = 0.8 * nyquist;
    }
    
    // now trim frequencies to in range
    // use list because bounds are by frequency rather than index
    // use variable-size data structures to prevent issues with rounding
    // based on calculation of where minimum index should exist
    List<Double> freqList = new LinkedList<Double>();
    Map<Double, Complex> numPSDMap = new HashMap<Double, Complex>();
    Map<Double, Complex> denomPSDMap = new HashMap<Double, Complex>();
    for (int i = 0; i < freqs.length; ++i) {
      
      if (freqs[i] < minFreq) {
        continue;
      }
      if (freqs[i] > maxFreq) {
        break;
      }
      
      freqList.add(freqs[i]);
      numPSDMap.put(freqs[i], numeratorPSD.getFFT()[i]);
      denomPSDMap.put(freqs[i], denominatorPSD.getFFT()[i]);
    }
    
    double zeroTarget; // frequency to set all curves to zero at
    if (lowFreq) {
      zeroTarget = 0.02;
    } else {
      zeroTarget = 1.0;
    }
    
    // Collections.sort(freqList); // done mostly for peace of mind
    
    int len = freqList.size(); // length of trimmed frequencies
    freqs = new double[len];
    // trim the PSDs to the data in the trimmed frequency range
    Complex[] numeratorPSDVals = new Complex[len];
    Complex[] denominatorPSDVals = new Complex[len];
    
    for (int i = 0; i < len; ++i) {
      freqs[i] = freqList.get(i);
      
      numeratorPSDVals[i] = numPSDMap.get(freqs[i]);
      denominatorPSDVals[i] = denomPSDMap.get(freqs[i]);
      
      if ( freqs[i] == 1.0 || 
          (freqs[i] > zeroTarget && freqs[i - 1] < zeroTarget) ) {
        normalIdx = i;
      }
    }
    
    // applied response. make sure to use the correct units (velocity)
    Complex[] appResponse = fitResponse.applyResponseToInput(freqs);
    for (int i = 0; i < appResponse.length; ++i) {
      Complex scaleFactor = new Complex(0., NumericUtils.TAU * freqs[i]);
      appResponse[i] = appResponse[i].divide(scaleFactor);
    }
    
    // calculated response from deconvolving calibration from signal
    // (this will be in displacement and need to be integrated)
    Complex[] estResponse = new Complex[len];
    for (int i = 0; i < estResponse.length; ++i) {
      Complex numer = numeratorPSDVals[i];
      Complex denom = denominatorPSDVals[i];
      estResponse[i] = numer.divide(denom);
      // convert from displacement to velocity
      Complex scaleFactor = new Complex(0., NumericUtils.TAU * freqs[i]);
      estResponse[i] = estResponse[i].multiply(scaleFactor);
    }
    
    // next, normalize estimated response
    String name = sensorOut.getName();
    XYSeries calcMag = new XYSeries("Calc. resp. (" + name + ") magnitude");
    XYSeries calcArg = new XYSeries("Calc. resp. (" + name + ") phase");
    
    // scaling values, used to set curve values to 0 at 1Hz
    Complex scaleValue = estResponse[normalIdx];
    double subtractBy = 10 * Math.log10( scaleValue.abs() );
    double rotateBy = NumericUtils.atanc(scaleValue);
    
    // data to fit poles to; first half of data is magnitudes of resp (dB)
    // second half of data is angles of resp (radians, scaled)
    double[] observedResult = new double[2 * estResponse.length];
    
    // prevent discontinuities in angle plots
    double phiPrev = 0.;
    
    double[] obsdAmps = new double[estResponse.length];
    
    for (int i = 0; i < estResponse.length; ++i) {
      
      int argIdx = estResponse.length + i;
      
      Complex estValue = estResponse[i];
      // estValue = estValue.subtract(scaleValue);
      double estValMag = estValue.abs();
      double phi = NumericUtils.atanc(estValue);
      phi -= rotateBy;
      
      phi = NumericUtils.unwrap(phi, phiPrev);
      // iterative step
      phiPrev = phi;

      phi = Math.toDegrees(phi);
      
      if ( Double.isNaN(estValMag) ) {
        observedResult[i] = 0;
        observedResult[argIdx] = 0;
      } else {
        
        obsdAmps[i] = estValMag / scaleValue.abs();
        observedResult[i] = 10 * Math.log10(estValMag);
        observedResult[i] -= subtractBy;
        
        double argument = phi;
        // argument /= rotateBy;
        // argument *= -1;
        observedResult[argIdx] = argument;
        
      }
      
      double xAxis;
      if (freqSpace) {
        xAxis = freqs[i];
      } else {
        xAxis = 1. / freqs[i];
      }
      
      calcMag.add(xAxis, observedResult[i]);
      calcArg.add(xAxis, observedResult[argIdx]);
    }
    
    // want to set up weight-scaling for the input so rotation doesn't dominate
    // solver's residual calculations, i.e., so no phase overfitting
    
    fireStateChange("Getting weighting....");
    
    maxArgWeight = 1.; maxMagWeight = 0.;
    Complex weightScaler = estResponse[normalIdx];
    double subtractWeight = 10 * Math.log10( weightScaler.abs() );
    double rotateWeight = NumericUtils.atanc(weightScaler);
    for (int i = 0; i < estResponse.length; ++i) {
      // int argIdx = i + appResponse.length;
      double magCandidate = 10 * Math.log10( estResponse[i].abs() );
      magCandidate -= subtractWeight;
      double phiCandidate = Math.abs( NumericUtils.atanc(estResponse[i]) );
      phiCandidate -= rotateWeight;
      if ( magCandidate > maxMagWeight ) {
        maxMagWeight = magCandidate;
      }
      if ( phiCandidate > maxArgWeight ) {
        maxArgWeight = phiCandidate;
      }
    }
    
    fireStateChange("Setting weight matrix...");
    // System.out.println(maxMagWeight);
    
    // we have the candidate mag and phase, now to turn them into weight values
    maxMagWeight = 1000. / maxMagWeight; // scale factor to weight over
    if (maxArgWeight != 0.) {
      maxArgWeight = 1./ maxArgWeight;
    }
    
    // weight matrix
    double[] weights = new double[observedResult.length];
    for (int i = 0; i < estResponse.length; ++i) {
      int argIdx = i + estResponse.length;
      double denom;
      if (!lowFreq) {
        denom = 100.;
        // give frequencies below 1 less weight in high-freq calibrations
        if (freqs[i] > 10.) {
          // for high enough freqs, make weighting (100/f^3) rather than 1/f;
          denom *= Math.pow(freqs[i], 2) / 100.;
        } else if (freqs[i] > 1.) {
          denom = freqs[i];
        }
      } else {
        if (freqs[i] < .01) {
          denom = freqs[i];
        } else {
          denom = .01;
        }
        
      }

      weights[i] = maxMagWeight / denom;
      weights[argIdx] = maxArgWeight / denom;
    }
    
    DiagonalMatrix weightMat = new DiagonalMatrix(weights);
    
    fireStateChange("Getting estimate and setting up solver...");
    
    // now to set up a solver for the params -- first, get the input variables
    // complex values are technically two variables, each a double
    // so, let's split up their real and im components and make each one a
    // variable. (we also need to ignore conjugate values, for constraints)
    RealVector initialGuess, initialPoleGuess, initialZeroGuess;
    
    initialPoleGuess = fitResponse.polesToVector(lowFreq, nyquist);
    initialZeroGuess = fitResponse.zerosToVector(lowFreq, nyquist);
    numZeros = initialZeroGuess.getDimension();
    initialGuess = initialZeroGuess.append(initialPoleGuess);
    
    //System.out.println(nyquist);
    //for (int i = 0; i < initialPoles.size(); ++i) {
    //  System.out.println( initialPoles.get(i).abs() / NumericUtils.TAU );
    //}
    
    // now, solve for the response that gets us the best-fit response curve
    // RealVector initialGuess = MatrixUtils.createRealVector(responseVariables);
    RealVector obsResVector = MatrixUtils.createRealVector(observedResult);
    
    MultivariateJacobianFunction jacobian = new MultivariateJacobianFunction() {
      
      int numIterations = 0;
      
      public Pair<RealVector, RealMatrix> value(final RealVector point) {
        ++numIterations;
        fireStateChange("Fitting, iteration count " + numIterations);
        Pair<RealVector, RealMatrix> pair = 
            jacobian(point);
        return pair;
      }
      
    };
    
    ConvergenceChecker<LeastSquaresProblem.Evaluation> svc = 
        new EvaluationRmsChecker(1.0E-14, 1.0E-14);
    
    LeastSquaresOptimizer optimizer = new LevenbergMarquardtOptimizer().
        withCostRelativeTolerance(1.0E-14).
        withParameterRelativeTolerance(1.0E-14);
    
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
        checker(svc).
        build();
    
    fireStateChange("Built least-squares problem; evaluating intial guess...");

    // residuals used to determine quality of solution convergence
    
    LeastSquaresProblem.Evaluation initEval = lsp.evaluate(initialGuess);
    initialResidual = initEval.getCost();
    
    fireStateChange("Got initial evaluation; running solver...");
    
    double[] initialValues =
        jacobian.value(initialGuess).getFirst().toArray();
    
    RealVector finalResultVector;

    boolean dontSolve = getSolverState(); // true if we should NOT run solver
    
    if (!dontSolve) {
      LeastSquaresOptimizer.Optimum optimum = optimizer.optimize(lsp);
      finalResultVector = optimum.getPoint();
      numIterations = optimum.getIterations();
    } else {
      finalResultVector = initialGuess;
    }
    
    LeastSquaresProblem.Evaluation optimum = lsp.evaluate(finalResultVector);
    fitResidual = optimum.getCost();
    double[] fitParams = optimum.getPoint().toArray();
    // get results from evaluating the function at the two points
    double[] fitValues = 
        jacobian.value( optimum.getPoint() ).getFirst().toArray();
    
    // double[] initResidList = initEval.getResiduals().toArray();
    // double[] fitResidList = optimum.getResiduals().toArray();
    XYSeries initResidMag = new XYSeries("Amplitude of init. residual");
    XYSeries initResidPhase = new XYSeries("Phase of init. residual");
    XYSeries fitResidMag = new XYSeries("Amplitude of fit residual");
    XYSeries fitResidPhase = new XYSeries("Phase of fit residual");
    
    // InstrumentResponse init = ds.getResponse(sensorOutIdx);
    
    fitResponse = fitResponse.buildResponseFromFitVector(
        fitParams, lowFreq, numZeros);
    fitPoles = fitResponse.getPoles();
    fitZeros = fitResponse.getZeros();
    
    fireStateChange("Compiling data...");
    
    for (int i = 0; i < freqs.length; ++i) {
      double xValue;
      if (freqSpace) {
        xValue = freqs[i];
      } else {
        xValue = 1. / freqs[i];
      }
      
      int argIdx = freqs.length + i;
      initMag.add(xValue, initialValues[i]);
      initArg.add(xValue, initialValues[argIdx]);
      fitMag.add(xValue, fitValues[i]);
      fitArg.add(xValue, fitValues[argIdx]);
      
      // Complex scaledInit = initTerms[i].subtract(init1Hz);
      // Complex scaledFit = fitTerms[i].subtract(fit1Hz);
      
      double initAmpNumer = Math.pow(10, initialValues[i]/10);
      double fitAmpNumer = Math.pow(10, fitValues[i]/10); 
      
      double obsAmpDbl = obsdAmps[i];
      if (obsAmpDbl == 0.) {
        obsAmpDbl = Double.MIN_VALUE;
      }
      
      double errInitMag = 100. * (initAmpNumer - obsAmpDbl) / obsAmpDbl; 
      double errFitMag = 100. * (fitAmpNumer - obsAmpDbl) / obsAmpDbl; 
      initResidMag.add(xValue, errInitMag);
      fitResidMag.add(xValue, errFitMag);
      
      double observedPhase = observedResult[argIdx];
      initResidPhase.add(xValue, initialValues[argIdx] - observedPhase);
      fitResidPhase.add(xValue, fitValues[argIdx] - observedPhase);
    }
    
    XYSeriesCollection xysc = new XYSeriesCollection();
    xysc.addSeries(initMag);
    xysc.addSeries(calcMag);
    if (!dontSolve) {
      xysc.addSeries(fitMag);
    }

    xySeriesData.add(xysc);
    
    xysc = new XYSeriesCollection();
    xysc.addSeries(initArg);
    xysc.addSeries(calcArg);    
    if (!dontSolve) {
      xysc.addSeries(fitArg);
    }
    xySeriesData.add(xysc);
    
    
    xysc = new XYSeriesCollection();
    xysc.addSeries(initResidMag);
    if (!dontSolve) {
      xysc.addSeries(fitResidMag);
    }
    xySeriesData.add(xysc);
    
    xysc = new XYSeriesCollection();    
    xysc.addSeries(initResidPhase);
    if (!dontSolve) {
      xysc.addSeries(fitResidPhase);
    }
    xySeriesData.add(xysc);
    
  }
  
  @Override
  public int blocksNeeded() {
    return 2;
  }
  
  /**
   * Backend function to set instrument response according to current
   * test variables (for best-fit calculation / forward difference) and
   * produce a response from that result. The passed response is copied on
   * start and is not modified directly. Which values (poles) are modified
   * depends on high or low frequency calibration setting.
   * @param variables values to set the instrument response to 
   * @return Doubles representing new response curve evaluation
   */
  private double[] evaluateResponse(double[] variables) {
    
    InstrumentResponse testResp = new InstrumentResponse(fitResponse);
    
    // prevent terrible case where, say, only high-freq poles above nyquist rate
    if ( variables.length > 0) {
      testResp = fitResponse.buildResponseFromFitVector(
          variables, lowFreq, numZeros);
    } else {
      System.out.println("NO VARIABLES TO SET. THIS IS AN ERROR.");
    }
    
    Complex[] appliedCurve = testResp.applyResponseToInput(freqs);
    
    // array is magnitudes, then arguments of complex number
    double[] curValue = new double[appliedCurve.length * 2];
    curValue[0] = 0.;
    curValue[appliedCurve.length] = 0.;
    
    
    Complex scaleFactor = new Complex(0., NumericUtils.TAU * freqs[normalIdx]);
    Complex scaleBy = 
        appliedCurve[normalIdx].divide(scaleFactor);
    double magScale = 10 * Math.log10( scaleBy.abs() );
    double argScale = NumericUtils.atanc(scaleBy);
    
    double phiPrev = 0.;
    if (lowFreq) {
      int argIdx = appliedCurve.length / 2;
      double startAngle = NumericUtils.atanc(appliedCurve[argIdx]);
      // note that we are still in units of radians at this point
      phiPrev = startAngle; // hopefully fix issue with phase scaling
    }
    
    // System.out.println(appliedCurve[0]);
    // now, do scaling and create the result vector (or, rather, array)
    for (int i = 0; i < appliedCurve.length; ++i) {
      
      int argIdx = appliedCurve.length + i;
      
      if (freqs[i] == 0.) {
        // this would be a divide by 0 error, let's just call the result 0
        // (because freqs is trimmed in main function, this is likely not
        // necessary, but better safe than sorry when it comes down to it)
        curValue[i] = 0.;
        curValue[argIdx] = 0.;
        continue;
      }
      
      // from acceleration to velocity
      Complex value = appliedCurve[i];
      scaleFactor = new Complex(0., NumericUtils.TAU * freqs[i]);
      value = value.divide(scaleFactor);
      
      // value = value.subtract(scaleBy);
      
      if ( value.equals(Complex.NaN) ) {
        // this shouldn't happen, but just in case, make sure it's 0;
        System.out.println("It's NaN: " + i+"; freq: "+freqs[i]);
        curValue[i] += 0;
        curValue[argIdx] = 0;
      } else {
        // System.out.println(value);
        double temp = 10 * Math.log10( value.abs() );
        temp -= magScale;
        curValue[i] = temp;
        
        double phi = NumericUtils.atanc(value) - argScale;
        phi = NumericUtils.unwrap(phi, phiPrev);
        phiPrev = phi; // iterative step for unwrap function
        
        // phi /= argScale;
        // phi /= TAU;
        curValue[argIdx] = Math.toDegrees(phi);
        
        // if line below is uncommented, we're not fitting angle
        // curValue[argIdx] = 0.;
      }
    }
    
    return curValue;
  }
  
  /**
   * Get the poles that the solver has found to best-fit the est. response
   * @return new poles that should improve fit over inputted response, as a list
   */
  public List<Complex> getFitPoles() {
    return getPoleSubList(fitPoles);
  }
  
  /**
   * Get the residual value from the solved response parameters
   * @return the residual of the solved-for poles (best-fit response)
   */
  public double getFitResidual() {
    return fitResidual;
  }
  
  /**
   * Get the zeros fitted from the experiment
   * @return List of zeros (complex numbers) that are used in best-fit curve
   */
  public List<Complex> getFitZeros() {
    return getZeroSubList(fitZeros);
  }
  
  /**
   * Get poles used in input response, for reference against best-fit poles 
   * @return poles taken from initial response file
   */
  public List<Complex> getInitialPoles() {
    return getPoleSubList(initialPoles);
  }
  
  /**
   * Get initial zeros from (nominal) response file used as input
   * @return zeros taken from initial response file
   */
  public List<Complex> getInitialZeros() {
    return getZeroSubList(initialZeros);
  }
  
  /**
   * Get the residual value of the initial response parameters
   * @return the residual of the initial poles from fed-in response
   */
  public double getInitResidual() {
    return initialResidual;
  }

  /**
   * Get the number of times the algorithm iterated to produce the optimum
   * response fit, from the underlying least squares solver
   * @return the number of iterations
   */
  public int getIterations() {
    return numIterations;
  }
  
  /**
   * Trim down the poles to those within the range of those being fit
   * @param polesToTrim Either fit or input poles, sorted by frequency
   * @return Sublist of data to be fed to output reports
   */
  private List<Complex> getPoleSubList(List<Complex> polesToTrim) {
    List<Complex> subList = new ArrayList<Complex>();  
    
    double peak = PEAK_MULTIPLIER * nyquist;
    
    for (int i = 0; i < polesToTrim.size(); ++i) {
      double freq = initialPoles.get(i).abs() / NumericUtils.TAU;
      
      if ( ( lowFreq && freq > 1. ) || ( !lowFreq && freq > peak ) ) {
        break;
      }
      if (!lowFreq && freq < 1.) {
        continue; // ignore b
      }
      
      subList.add( polesToTrim.get(i) );
    }
    
    return subList;
  }

  /**
   * Used to determine whether to run the solver or not; disabling the solver
   * is useful for determining the quality of a given calibration function
   * @return True if the solver is to be run
   */
  public boolean getSolverState() {
    return SKIP_SOLVING;
  }
  
  /**
   * Get the values used to weight the residual calculation function.
   * The first value is the magnitude weighting, the second is phase.
   * @return Weighting values for least-squared error terms of
   */
  public double[] getWeights() {
    return new double[]{maxMagWeight, maxArgWeight};
  }
  
  private List<Complex> getZeroSubList(List<Complex> zerosToTrim) {
    List<Complex> subList = new ArrayList<Complex>();
    
    double peak = PEAK_MULTIPLIER * nyquist;
    
    for (int i = 0; i < zerosToTrim.size(); ++i) {
      double freq = initialZeros.get(i).abs() / NumericUtils.TAU;
      
      if ( ( lowFreq && freq > 1. ) || ( !lowFreq && freq > peak ) ) {
        break;
      }
      if (!lowFreq && freq < 1. || freq == 0.) {
        continue;
      }
      
      subList.add( zerosToTrim.get(i) );
    }
    
    return subList;
  }
  
  @Override
  public boolean hasEnoughData(DataStore ds) {
    return ( ds.blockIsSet(0) && ds.bothComponentsSet(1) );
  }
  
  /**
   * Function to run evaluation and forward difference for Jacobian 
   * approximation given a set of points to set as response. 
   * Mainly a wrapper for the evaluateResponse function.
   * @param variables Values to set the response's poles to
   * @return RealVector with evaluation at current response value and 
   * RealMatrix with forward difference of that response (Jacobian)
   */
  private Pair<RealVector, RealMatrix> 
  jacobian(RealVector variables) {
    
    // variables = validate(variables);
    
    int numVars = variables.getDimension();
    
    double[] currentVars = new double[numVars];
    
    for (int i = 0; i < numVars; ++i) {
      currentVars[i] = variables.getEntry(i);
    }
    
    double[] mag = evaluateResponse(currentVars);
    
    double[][] jacobian = new double[mag.length][numVars];
    // now take the forward difference of each value 
    for (int i = 0; i < numVars; ++i) {
      
      if (i % 2 == 1 && currentVars[i] == 0.) {
        // this is a zero pole. don't bother changing it
        for (int j = 0; j < mag.length; ++j) {
          jacobian[j][i] = 0.;
        }
        continue;
      }
      
      double[] changedVars = new double[currentVars.length];
      for (int j = 0; j < currentVars.length; ++j) {
        changedVars[j] = currentVars[j];
      }
      
      double diffX = changedVars[i] + DELTA;
      
      // real-value pole components must be less than zero
      if (diffX > 0. && (i % 2) == 0.) {
        diffX = 0.;
      }
      changedVars[i] = diffX;
      
      double[] diffY = 
          evaluateResponse(changedVars);
      
      for (int j = 0; j < diffY.length; ++j) {
        if (changedVars[i] - currentVars[i] == 0.) {
          jacobian[j][i] = 0.;
        } else {
          jacobian[j][i] = diffY[j] - mag[j];
          jacobian[j][i] /= changedVars[i] - currentVars[i];
        }
        /*
        if ( (i % 2) == 0 && currentVars[i] > 0) {
          // enforce that real values of poles must be negative
          jacobian[j][i] = -1;
        }
        */
      }
      
    }
    
    RealVector result = MatrixUtils.createRealVector(mag);
    RealMatrix jMat = MatrixUtils.createRealMatrix(jacobian);
    
    return new Pair<RealVector, RealMatrix>(result, jMat);
    
  }
  
  @Override
  public int[] listActiveResponseIndices() {
    // NOTE: not used by corresponding panel, overrides with active indices
    // of components in the combo-box
    return new int[]{sensorOutIdx};
  }
  
  /**
   * Determines which poles to fit when doing the response curve fitting;
   * low frequency calibrations set the first two poles; high frequency
   * calibrations set the remaining poles
   * @param lowFreq True if a low frequency calibration is to be used
   */
  public void setLowFreq(boolean lowFreq) {
    this.lowFreq = lowFreq;
  }
 
  /**
   * Simple validator method to enforce poles to be negative for their values
   * (Since imaginary values are stored in resps as their value and complex
   * conjugate, we can mandate this for all values in the vector, though only
   * real components are strictly required to be negative).
   * @param poleParams RealVector of parameters to be evaluated by solver
   * @return Vector of parameters but with components all negative
   */
  public RealVector validate(RealVector poleParams) {
    for (int i = 0; i < poleParams.getDimension(); ++i) {
      double value = poleParams.getEntry(i);
      if (value > 0 && (i % 2) == 0) {
        // even index means this is a real-value vector entry
        // if it's above zero, put it back below zero
        poleParams.setEntry(i, -Double.MIN_VALUE);
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
