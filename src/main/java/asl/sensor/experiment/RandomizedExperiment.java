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
  private static final double CUTOFF = 1. / 1000.; // used w/ KS-54000
  
  
  // To whomever has to maintain this code after I'm gone:
  // I'm sorry, I'm so so sorry
  // The vast bulk of this class is just adapter functions to get the
  // least-squares/LM solver to play nice with (constrained) complex values
  // and to get those to play nice with the instrument-response class
  
  /**
   * Take the parameters given by the fit function and build a response from 
   * them
   * @param fitParams Array of doubles representing alternately
   * complex real and imaginary parts of the results from the cal solver,
   * both pole and zero values
   * @param lowFreq True if the features represent the result from solving
   * for low frequency cal parameters
   * @param numZeros How much of input array is zero values, i.e. two zeros
   * being fit means this is 4 (each zero has real and complex value)
   * @param nyquist Nyquist rate of sensor data the response is associated with
   * @return Instrument response with fit parameters set as new poles, zeros
   */
  public static InstrumentResponse 
  fitResultToResp(double[] fitParams, InstrumentResponse ir, 
                  boolean lowFreq, int numZeros, double nyquist) {
    
    InstrumentResponse fitResp = new InstrumentResponse(ir);
    
    double[] zeros = new double[numZeros];
    for (int i = 0; i < zeros.length; ++i) {
      zeros[i] = fitParams[i];
    }
    
    fitResp = zerosToResp(zeros, fitResp, lowFreq, nyquist);
    
    double[] poles = new double[fitParams.length - zeros.length];
    for (int i = 0; i < poles.length; ++i) {
      int fitIdx = i + zeros.length;
      poles[i] = fitParams[fitIdx];
    }
    
    fitResp = polesToResp(poles, fitResp, lowFreq, nyquist);
    
    return fitResp;
    
  }
  
  /**
   * Quick check if the lowest-frequency pole is the KS54000, which should not
   * be fit in this experiment's operation, merely ignored
   * @param poles List of poles from response file, sorted
   * @return True if the lowest-frequency pole is too low to try to fit
   */
  private static boolean isKS54000(List<Complex> poles) {
    if ( ( poles.get(0).abs() / NumericUtils.TAU ) < CUTOFF ) {
      // first two poles are low-frequency
      return true;
    }
    
    return false;
  }
  
  /**
   * Convert a list of variables to a set of poles to be applied to an
   * instrument response file. This variable must have an even number of
   * entries or the function will return an out-of-bounds exception.
   * Variable pairs are combined to create complex numbers, with the first
   * entry representing the real component of the complex number and the second
   * representing the imaginary component. If this is a low-frequency
   * calibration, there should be 4 variables; the lowest two poles are the
   * ones which will be modified. If this is a high-frequency calibration,
   * then there should be 2*(pole.size() - 2) variables as all but the lowest
   * two poles will be modified.
   * Note that the passed response file is not changed or modified by this
   * process; a deep copy of the response is made, then that is modified and 
   * returned.
   * @param variables List of doubles representing complex number real and
   * imaginary components of poles
   * @param ir Response file to take as source
   * @param lowFreq True if the bottom low-frequency poles should be modified, 
   * false if high-frequency poles up to Nyquist rate should be
   * @param nyquist Nyquist rate of input data
   * @return New response file with the altered poles
   */
  public static InstrumentResponse 
  polesToResp(double[] variables, InstrumentResponse ir, 
              boolean lowFreq, double nyquist) {
    
    int numVars = variables.length;
    InstrumentResponse testResp = new InstrumentResponse(ir);
    
    List<Complex> poleList = new ArrayList<Complex>( testResp.getPoles() );
    List<Complex> builtPoles = new ArrayList<Complex>();
    
    
    if (!lowFreq) {
      // first, add the low-frequency poles if we're doing high-freq cal
      for (int i = 0; i < poleList.size(); ++i) {
        Complex pole = poleList.get(i);
        if ( pole.abs() / NumericUtils.TAU > 1. ) {
          break;
        }
        builtPoles.add(pole);
      }
      
    } else if ( isKS54000(poleList) ) {
      // in the low-frequency cal we've ignored this low-freq damping pole
      // so we need to add it here
      builtPoles.add( poleList.get(0) );
    }
    
    // now add the poles under consideration for fit
    // these are the high-frequency poles if we're doing high-frequency cal
    // or the low-frequency poles
    for (int i = 0; i < numVars; i += 2) {
      Complex c = new Complex(variables[i], variables[i+1]);
      builtPoles.add(c);
      if ( variables[i+1] != 0. ) {
        builtPoles.add( c.conjugate() );
      }
    }
    
    if (lowFreq) {
      // add ALL poles above 1Hz
      for (int i = 0; i < poleList.size(); ++i) {
        Complex pole = poleList.get(i);
        if ( pole.abs() / NumericUtils.TAU >= 1.) {
          builtPoles.add(pole);
        }
      }
    } else {
      // add the poles above the Nyquist rate when doing high-frequency cal
      for (int i = 0; i < poleList.size(); ++i) {
        Complex pole = poleList.get(i);
        if ( pole.abs() / NumericUtils.TAU >= nyquist ) {
          builtPoles.add(pole);
        }
      }
    }
    
    // System.out.println(poleList);
    // get the result for the input value
    testResp.setPoles(builtPoles);
    return testResp;
  }
  
  /**
   * Converts pole values into variables to be fit by the solver, ignoring the
   * inclusion of a pole's complex conjugate pair in the list. If the imaginary
   * component of a pole is 0, then no such conjugate pole exists.
   * For low frequency calibrations, only poles below 1Hz are used. For high
   * frequency calibrations, only poles between 1Hz and the Nyquist rate are
   * used.
   * Note that even vector indices (0, 2, 4...) are the real part of each pole
   * and the odd vector indices (1, 3, 5...) are the corresponding imaginary
   * part of each pole. That is, the vector index pair (0, 1) defines the first
   * high frequency pole's complex value.
   * @param poles List of poles to extract variables from (from input response).
   * These are sorted when loaded in from a RESP file, which this method expects
   * @param lowFreq True if a low frequency [long period] cal is being analysed
   * @param nyquist Nyquist rate of cal signal
   * @return RealVector containing each non-conjugate pole's components
   */
  public static RealVector 
  polesToVector(List<Complex> poles, boolean lowFreq, double nyquist) {
    
    // shame we can't use complex numbers in the LSP solver, huh
    
    // create a list of doubles that are the non-conjugate elements from list
    // of poles, to convert to array and then vector format
    List<Double> componentList = new ArrayList<Double>();
    
    // WARNING! Pole list should always be sorted by magnitude!
    
    int start = 0;
    
    if ( isKS54000(poles) ) {
      start = 1;
    }
    
    for (int i = start; i < poles.size(); ++i) {
      
      if ( !lowFreq && ( poles.get(i).abs() / NumericUtils.TAU < 1. ) ) {
        // don't include poles below 1Hz in high-frequency calibration
        continue;
      }
      
      if ( lowFreq && ( poles.get(i).abs() / NumericUtils.TAU > 1. ) ) {
        // only do low frequency calibrations on poles up to 
        break;
      }
      
      if ( !lowFreq && ( poles.get(i).abs() / NumericUtils.TAU >= nyquist ) ) {
        // don't fit poles above nyquist rate of sensor output
        break;
      }
      
      double realPart = poles.get(i).getReal();
      double imagPart = poles.get(i).getImaginary();
      
      componentList.add(realPart);
      componentList.add(imagPart);
      
      if (imagPart != 0.) {
        // next value is complex conjugate of this one, so skip it
        ++i;
      }
      
    }

    // turn into array to be turned into vector
    // can't use toArray because List doesn't use primitive double objects
    double[] responseVariables = new double[componentList.size()];
    for (int i = 0; i < responseVariables.length; ++i) {
      responseVariables[i] = componentList.get(i);
    }
    
    return MatrixUtils.createRealVector(responseVariables);
  
  }
  
  /**
   * Converts zero values from the solver into parameters as part of a
   * response object, and creates an instrument response from them. Zeros above
   * the nyquist rate have not been fit, and so will need to come from response
   * if a high-frequency cal was done (this is not an issue for low-frequency
   * cals as they already take the high-frequency values from the response)
   * @param variables List of doubles representing complex number real and
   * imaginary components of zeros
   * @param ir Response file to take as source
   * @param lowFreq True if the bottom low-frequency zeros should be modified, 
   * false if high-frequency poles up to Nyquist rate should be
   * @param nyquist Nyquist rate of input data
   * @return New response file with the altered zeros
   */
  public static InstrumentResponse 
  zerosToResp(double[] variables, InstrumentResponse ir, 
              boolean lowFreq, double nyquist) {
   
    int numVars = variables.length;
    InstrumentResponse testResp = new InstrumentResponse(ir);
    
    List<Complex> zeroList = new ArrayList<Complex>( testResp.getZeros() );
    List<Complex> builtZeros = new ArrayList<Complex>();
    
    // first, add the literally zero values (no more than 2)
    for (int i = 0; i < 2; ++i) {
      Complex zero = zeroList.get(i);
      if ( zero.abs() > 0. ) {
        break;
      }
      builtZeros.add(zero);
    }
    
    // now add the low-frequency zeros if they're not being fit
    if (!lowFreq) {
      // start by ignoring the zeros we just added
      for (int i = builtZeros.size(); i < zeroList.size(); ++i) {
        Complex zero = zeroList.get(i);
        // same cutoff criterias as in vector conversions and such
        if ( zero.abs() / NumericUtils.TAU > 1. ) {
          break;
        }
        builtZeros.add(zero);
      }
    }
    
    // now add the poles under consideration for fit
    // these are the high-frequency poles if we're doing high-frequency cal
    // or the low-frequency poles
    for (int i = 0; i < numVars; i += 2) {
      Complex c = new Complex(variables[i], variables[i+1]);
      builtZeros.add(c);
      if ( variables[i+1] != 0. ) {
        builtZeros.add( c.conjugate() );
      }
    }
    
    if (lowFreq) {
      // now for the high-frequency zeros
      // we'll start from zero and only add those above 1Hz
      for (int i = 0; i < zeroList.size(); ++i) {
        Complex zero = zeroList.get(i);
        if ( zero.abs() / NumericUtils.TAU >= 1. ) {
          builtZeros.add(zero);
        }
      }
    } else {
      // add the zeros above the Nyquist rate for high-freq cal case
      for (int i = 0; i < zeroList.size(); ++i) {
        Complex pole = zeroList.get(i);
        if ( pole.abs() / NumericUtils.TAU >= nyquist ) {
          builtZeros.add(pole);
        }
      }
    }
    
    // now add any zeros in the response above the nyquist cut-off
    
    testResp.setZeros(builtZeros);
    return testResp;
  }
  
  /**
   * Create a vector of the zeros to be fit by the function as a vector of
   * alternating real and imaginary values. We ignore literally zero-valued
   * zeros, zeros outside of the range of fit for frequency, and constrain
   * zeros with complex conjugates
   * @param zeros List of zeros to be made into a vector by function, sorted
   * @param lowFreq True if a low-frequency cal is used
   * @param nyquist Nyquist rate of data associated with sensor's zeros
   * @return RealVector of zeros under analysis
   */
  public static RealVector 
  zerosToVector(List<Complex> zeros, boolean lowFreq, double nyquist) {
    
    // also a shame we can't use multiple vectors of variables in the solver?
    
    // create a list of doubles that are the non-conjugate elements from list
    // of poles, to convert to array and then vector format
    List<Double> componentList = new ArrayList<Double>();
    
    // WARNING! Zeros list should always be sorted!
    
    for (int i = 0; i < zeros.size(); ++i) {
      
      if ( zeros.get(i).abs() == 0. ) {
        // ignore zeros that are literally zero-valued
        continue;
      }
      
      if ( !lowFreq && ( zeros.get(i).abs() / NumericUtils.TAU < 1. ) ) {
        // don't include zeros 1Hz in high-frequency calibration
        continue;
      }
      
      if ( lowFreq && ( zeros.get(i).abs() / NumericUtils.TAU > 1. ) ) {
        // only do low frequency calibrations on zeros up to 1Hz
        break;
      }
      
      if ( !lowFreq && ( zeros.get(i).abs() / NumericUtils.TAU > nyquist) ) {
        // don't fit zeros above nyquist rate of sensor output
        break;
      }
      
      double realPart = zeros.get(i).getReal();
      double imagPart = zeros.get(i).getImaginary();
      
      componentList.add(realPart);
      componentList.add(imagPart);
      
      if (imagPart != 0.) {
        // next value is complex conjugate of this one, so skip it
        ++i;
      }
      
    }

    // turn into array to be turned into vector
    // can't use toArray because List doesn't use primitive double objects
    double[] responseVariables = new double[componentList.size()];
    for (int i = 0; i < responseVariables.length; ++i) {
      responseVariables[i] = componentList.get(i);
    }
    
    return MatrixUtils.createRealVector(responseVariables);
    
  }
  
  // yes, folks, appx. half this class is static methods designed to serve as
  // adapter functions to convert from poles and zeros to fittable vectors and
  // vice-versa
  
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
  
  private double maxMagWeight, maxArgWeight; // max values of magnitude, phase
  
  private int normalIdx; // location of value to set to 0 in curves for scaling
  private int numZeros; // how many entries in parameter vector define zeros
  private int sensorOutIdx; // location to load response from?
  
  public RandomizedExperiment() {
    super();
    lowFreq = false;
    normalIdx = 0;
  }
  
  public boolean getSolverState() {
    return SKIP_SOLVING;
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
    
    if ( ds.getBlock(sensorOutIdx).getName().equals( calib.getName() ) ) {
      sensorOutIdx = ds.getXthFullyLoadedIndex(2);
    }

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
    
    // also, use those frequencies to get the applied response to input
    
    FFTResult numeratorPSD = FFTResult.spectralCalc(sensorOut, calib);
    FFTResult denominatorPSD = FFTResult.spectralCalc(calib, calib);
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
      appResponse[i] = appResponse[i].divide(NumericUtils.TAU * freqs[i]);
    }
    
    // calculated response from deconvolving calibration from signal
    // (this will be in displacement and need to be integrated)
    Complex[] estResponse = new Complex[len];
    for (int i = 0; i < estResponse.length; ++i) {
      Complex numer = numeratorPSDVals[i];
      Complex denom = denominatorPSDVals[i];
      estResponse[i] = numer.divide(denom);
      // convert from displacement to velocity
      estResponse[i] = estResponse[i].multiply(NumericUtils.TAU * freqs[i]);
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
        
        observedResult[i] = 10 * Math.log10(estValMag);
        observedResult[i] -= subtractBy;
        
        double argument = phi;
        // argument /= rotateBy;
        // argument *= -1;
        observedResult[argIdx] = argument;
        
      }
      
      calcMag.add(freqs[i], observedResult[i]);
      calcArg.add(freqs[i], observedResult[argIdx]);
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
    maxMagWeight = 10. / maxMagWeight;
    maxArgWeight = 1./ maxArgWeight;
    
    // weight matrix
    double[] weights = new double[observedResult.length];
    for (int i = 0; i < estResponse.length; ++i) {
      int argIdx = i + estResponse.length;
      // weights[i] = 1 / Math.pow(10, maxMagWeight);
      // weights[i] = 10000;
      weights[i] = maxMagWeight; // scale by 100 due to peak adjustment
      weights[argIdx] = maxArgWeight;
    }
    
    DiagonalMatrix weightMat = new DiagonalMatrix(weights);
    
    fireStateChange("Getting estimate and setting up solver...");
    
    // now to set up a solver for the params -- first, get the input variables
    // complex values are technically two variables, each a double
    // so, let's split up their real and im components and make each one a
    // variable. (we also need to ignore conjugate values, for constraints)
    RealVector initialGuess, initialPoleGuess, initialZeroGuess;
    
    initialPoleGuess = polesToVector(initialPoles, lowFreq, nyquist);
    initialZeroGuess = zerosToVector(initialZeros, lowFreq, nyquist);
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
      
      public Pair<RealVector, RealMatrix> value(final RealVector point) {
        Pair<RealVector, RealMatrix> pair = 
            jacobian(point);
        return pair;
      }
    };
    
    ConvergenceChecker<LeastSquaresProblem.Evaluation> svc = 
        new EvaluationRmsChecker(1.0E-7, 1.0E-7);
    
    LeastSquaresOptimizer optimizer = new LevenbergMarquardtOptimizer().
        withCostRelativeTolerance(1.0E-7).
        withParameterRelativeTolerance(1.0E-7);
    
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
    
    XYSeries initResidMag = new XYSeries("Input resp. mag residual");
    XYSeries initResidPhase = new XYSeries("Input resp. phase residual");
    
    RealVector finalResultVector;

    boolean dontSolve = getSolverState(); // true if we should NOT run solver
    
    if (!dontSolve) {
      LeastSquaresOptimizer.Optimum optimum = optimizer.optimize(lsp);
      finalResultVector = optimum.getPoint();
    } else {
      finalResultVector = initialGuess;
    }

    LeastSquaresProblem.Evaluation optimum = lsp.evaluate(finalResultVector);
    fitResidual = optimum.getCost();
    double[] fitParams = optimum.getPoint().toArray();
    // get results from evaluating the function at the two points
    double[] fitValues = 
        jacobian.value( optimum.getPoint() ).getFirst().toArray();
    
    double[] initResidList = initEval.getResiduals().toArray();
    double[] fitResidList = optimum.getResiduals().toArray();
    

    XYSeries fitResidMag = new XYSeries("Fit resp. mag sqd. error");
    XYSeries fitResidPhase = new XYSeries("Fit resp. phase sqd. error");
    
    fitResponse = 
        fitResultToResp(fitParams, fitResponse, lowFreq, numZeros, nyquist);
    fitPoles = fitResponse.getPoles();
    fitZeros = fitResponse.getZeros();
    
    fireStateChange("Compiling data...");
    
    for (int i = 0; i < freqs.length; ++i) {
      int argIdx = freqs.length + i;
      initMag.add(freqs[i], initialValues[i]);
      initArg.add(freqs[i], initialValues[argIdx]);
      fitMag.add(freqs[i], fitValues[i]);
      fitArg.add(freqs[i], fitValues[argIdx]);
      
      initResidMag.add( freqs[i], Math.pow(initResidList[i], 2) );
      initResidPhase.add( freqs[i], Math.pow(initResidList[argIdx], 2) );
      fitResidMag.add( freqs[i], Math.pow(fitResidList[i], 2) );
      fitResidPhase.add( freqs[i], Math.pow(fitResidList[argIdx], 2) );
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
   * @param ir InstrumentResponse that will be copied 
   * @return Doubles representing new response curve evaluation
   */
  private double[] 
  evaluateResponse(double[] variables) {
    
    InstrumentResponse testResp = new InstrumentResponse(fitResponse);
    
    // prevent terrible case where, say, only high-freq poles above nyquist rate
    if (variables.length > 0) {
      testResp = 
          fitResultToResp(variables, testResp, lowFreq, numZeros, nyquist);
    }
    
    Complex[] appliedCurve = testResp.applyResponseToInput(freqs);
    
    // array is magnitudes, then arguments of complex number
    double[] curValue = new double[appliedCurve.length * 2];
    curValue[0] = 0.;
    curValue[appliedCurve.length] = 0.;
    
    Complex scaleBy = 
        appliedCurve[normalIdx].divide(NumericUtils.TAU * freqs[normalIdx]);
    double magScale = 10 * Math.log10( scaleBy.abs() );
    double argScale = NumericUtils.atanc(scaleBy);
    
    double phiPrev = 0.;
    
    // System.out.println(appliedCurve[0]);
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
      value = value.divide(NumericUtils.TAU * freqs[i]);
      
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
  
  public List<Complex> getFitZeros() {
    return getZeroSubList(fitZeros);
  }
  
  /**
   * Get poles used in input response, for reference against best-fit poles 
   * @return original poles being modified by the solver in calib. processing,
   * as a list
   */
  public List<Complex> getInitialPoles() {
    return getPoleSubList(initialPoles);
  }
  
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
   * Trim down the poles to those within the range of those being fit
   * @param polesToTrim Either fit or input poles, sorted by frequency
   * @return Sublist of data to be fed to output reports
   */
  private List<Complex> getPoleSubList(List<Complex> polesToTrim) {
    
    List<Complex> subList = new ArrayList<Complex>();  
    
    int start = 0;
    
    if ( isKS54000(polesToTrim) ) {
      start = 1; // ignore first pole of KS54000
    }
    
    for (int i = start; i < polesToTrim.size(); ++i) {
      
      double freq = initialPoles.get(i).abs() / NumericUtils.TAU;
      
      if ( ( lowFreq && freq > 1. ) || ( !lowFreq && freq > nyquist ) ) {
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
   * Get the values used to weight the residual calculation function.
   * The first value is the magnitude weighting, the second is phase.
   * @return Weighting values for least-squared error terms of
   */
  public double[] getWeights() {
    return new double[]{maxMagWeight, maxArgWeight};
  }
  
  private List<Complex> getZeroSubList(List<Complex> zerosToTrim) {
    
    List<Complex> subList = new ArrayList<Complex>();
    
    for (int i = 0; i < zerosToTrim.size(); ++i) {
      
      double freq = initialZeros.get(i).abs() / NumericUtils.TAU;
      
      if ( ( lowFreq && freq > 1. ) || ( !lowFreq && freq > nyquist ) ) {
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
   * @param numZeros How much of input vector is zeros of response
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
      if (i > numZeros && diffX > 0. && (i % 2) == 0.) {
        diffX = 0.;
      }
      changedVars[i] = diffX;
      
      double[] diffY = 
          evaluateResponse(changedVars);
      
      for (int j = 0; j < diffY.length; ++j) {
        jacobian[j][i] = diffY[j] - mag[j];
        jacobian[j][i] /= changedVars[i] - currentVars[i];
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
    for (int i = numZeros; i < poleParams.getDimension(); ++i) {
      double value = poleParams.getEntry(i);
      if (value > 0 && (i % 2) == 0) {
        // even index means this is a real-value vector entry
        // if it's above zero, set it back to zero
        poleParams.setEntry(i, 0.);
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
