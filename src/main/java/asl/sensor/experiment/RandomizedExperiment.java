package asl.sensor.experiment;

import java.util.ArrayList;
import java.util.Collections;
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
 * @author akearns
 *
 */
public class RandomizedExperiment extends Experiment {

  private static final double DELTA = 1E-7;
  
  /**
   * Convert high-frequency pole values into variables to be fit by the solver,
   * ignoring the inclusion of a pole's complex conjugate pair in the list.
   * There may be one or more poles in the list that have imaginary component 0.
   * If this is the case, then no such complex conjugate exists in the pole
   * list, and so additional handling is done to prevent skipping over the next
   * valid pole. 
   * Note that even vector indices (0, 2, 4...) are the real part of each pole
   * and the odd vector indices (1, 3, 5...) are the corresponding imaginary
   * part of each pole. That is, the vector index pair (0, 1) defines the first
   * high frequency pole's complex value. 
   * @param poles List of poles to extract variables from (from input response)
   * @return
   */
  public static RealVector highFreqPolesToVector(List<Complex> poles) {
    
    // shame we can't use complex numbers in the LSP solver, huh
    
    int startIdx = 2; // starting index for non-conjugate pole counting
    if ( poles.get(0).getImaginary() == 0. ) {
      // start at 1 if the 0th pole has no conjugate
      startIdx = 1;
    }
    
    // create a list of doubles that are the non-conjugate elements from list
    // of poles, to convert to array and then vector format
    List<Double> componentList = new ArrayList<Double>();
    
    for (int i = startIdx; i < poles.size(); ++i) {
      
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
   * Convert low-frequency pole values into variables to be fit by the solver.
   * Because there are only two such poles at most, and they are complex
   * conjugates of each other, we need only get the first real and imaginary
   * values of the function, stored as a vector of doubles.
   * Note that the first vector index (0) is the real part, 
   * and the second (1) is the imaginary part.
   * @param poles List of poles to extract variables from (from input response)
   * @return Vectorized representation of complex-value real and imaginary 
   * components
   */
  public static RealVector lowFreqPolesToVector(List<Complex> poles) {
    // we only need to get one pole; second low-freq pole is complex conjugate
    double[] responseVariables = new double[2];
    
    Complex cpx = poles.get(0);
    responseVariables[0] = cpx.getReal();
    responseVariables[1] = cpx.getImaginary();
    // cell [2] would be same as [0]
    // and cell[3] would be [-1]
    
    return MatrixUtils.createRealVector(responseVariables);
   
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
   * @param ir Response file which the list of pole variables will modify.
   * @param lowFreq True if the bottom 2 poles should be modified, false if
   * all other poles should be modified
   * @return New response file with the altered poles
   */
  public static InstrumentResponse polesToResp(double[] variables, 
      InstrumentResponse ir, boolean lowFreq) {
    
    int numVars = variables.length;
    InstrumentResponse testResp = new InstrumentResponse(ir);
    
    List<Complex> poleList = new ArrayList<Complex>( testResp.getPoles() );
    List<Complex> builtPoles = new ArrayList<Complex>();
    
    for (int i = 0; i < numVars; i += 2) {
      Complex c = new Complex(variables[i], variables[i+1]);
      builtPoles.add(c);
      if ( variables[i+1] != 0. ) {
        builtPoles.add( c.conjugate() );
      }
    }
    
    if (lowFreq) {
      for (int i = 0; i < builtPoles.size(); ++i) {
        poleList.set( i, builtPoles.get(i) );
      }
    } else {
      int offset = 2; // assume first two poles are complex conjugates
      if ( poleList.get(0).getImaginary() == 0. ) {
        // the first two poles are not complex conjugates
        offset = 1;
      }
      for (int i = 0; i < poleList.size() - offset; ++i) {
        poleList.set( i + offset, builtPoles.get(i) );
      }
    }
    // System.out.println(poleList);
    // get the result for the input value
    testResp.setPoles(poleList);
    return testResp;
    
  }
  private double initialResidual, fitResidual;
  private List<Complex> inputPoles;
  private List<Complex> fitPoles;
  
  private boolean lowFreq; // fit the low- or high-frequency poles?
  
  private InstrumentResponse fitResponse;
  
  private double[] freqs;
  
  private int oneHzIdx;
  
  private String responseName;
  
  public RandomizedExperiment() {
    super();
    lowFreq = false;
    oneHzIdx = 0;
  }
  
  @Override
  protected void backend(DataStore ds) {
    
    // construct response plot
    DataBlock calib = ds.getXthLoadedBlock(1);
    int sensorOutIndex = ds.getXthFullyLoadedIndex(1);
    
    if ( ds.getBlock(sensorOutIndex).getName().equals( calib.getName() ) ) {
      sensorOutIndex = ds.getXthFullyLoadedIndex(2);
    }
    
    DataBlock sensorOut = ds.getBlock(sensorOutIndex);
    fitResponse = new InstrumentResponse( ds.getResponse(sensorOutIndex) );
    responseName = fitResponse.getName();
    
    inputPoles = new ArrayList<Complex>( fitResponse.getPoles() );
    fitPoles = new ArrayList<Complex>( fitResponse.getPoles() );
    
    // get the plots of the calculated response from deconvolution
    // PSD(out, in) / PSD(in, in) gives us PSD(out) / PSD(in) while removing
    // imaginary terms from the denominator due to multiplication with the
    // complex conjugate
    // PSD(out) / PSD(in) is the response curve (i.e., deconvolution)
    
    // also, use those frequencies to get the applied response to input
    
    FFTResult numeratorPSD = FFTResult.spectralCalc(sensorOut, calib);
    FFTResult denominatorPSD = FFTResult.spectralCalc(calib, calib);
    
    // TODO: shift back to .5 if this is too extreme a cut
    double proportion = 0.4; // get up to .5 of Nyquist freq due to noise
    // (used to be .8 of frequency but this included too much noise)
    double minFreq = .2; // lower bound of .2 Hz (5s period) due to noise
    
    freqs = numeratorPSD.getFreqs();    

    int len = (int) (freqs.length * proportion); // index of peak freq to check
    
    // trim down the frequency array to the specified range
    // we use a list because the lower bound is not fixed by index
    
    // use variable-size data structures to prevent issues with rounding
    // based on calculation of where minimum index should exist
    List<Double> freqList = new LinkedList<Double>();
    Map<Double, Complex> numPSDMap = new HashMap<Double, Complex>();
    Map<Double, Complex> denomPSDMap = new HashMap<Double, Complex>();
    for (int i = 0; i < len; ++i) {
      if (freqs[i] < minFreq) {
        continue;
      }
      
      freqList.add(freqs[i]);
      numPSDMap.put(freqs[i], numeratorPSD.getFFT()[i]);
      denomPSDMap.put(freqs[i], denominatorPSD.getFFT()[i]);
    }
    
    Collections.sort(freqList); // done mostly for peace of mind
    
    len = freqList.size(); // now len is length of trimmed frequencies
    freqs = new double[len];
    // trim the PSDs to the data in the trimmed frequency range
    Complex[] numeratorPSDVals = new Complex[len];
    Complex[] denominatorPSDVals = new Complex[len];
    
    for (int i = 0; i < len; ++i) {
      freqs[i] = freqList.get(i);
      
      numeratorPSDVals[i] = numPSDMap.get(freqs[i]);
      denominatorPSDVals[i] = denomPSDMap.get(freqs[i]);
      
      if ( freqs[i] == 1.0 || (freqs[i] > 1.0 && freqs[i - 1] < 1.0) ) {
        oneHzIdx = i;
      }
    }
    
    // applied response. make sure to use the correct units (velocity)
    Complex[] appResponse = fitResponse.applyResponseToInput(freqs);
    for (int i = 0; i < appResponse.length; ++i) {
      appResponse[i] = appResponse[i].divide( 2 * Math.PI * freqs[i] );
      
    }
    
    // calculated response from deconvolving calibration from signal
    Complex[] estimatedResponse = new Complex[len];
    for (int i = 0; i < estimatedResponse.length; ++i) {
      Complex numer = numeratorPSDVals[i];
      Complex denom = denominatorPSDVals[i];
      estimatedResponse[i] = numer.divide(denom);
      // convert from displacement to velocity
      estimatedResponse[i] = 
          estimatedResponse[i].multiply(2 * Math.PI * freqs[i]);
    }
    
    
    // next, normalize estimated response
    XYSeries calcMag = new XYSeries("Calc. resp. magnitude");
    XYSeries calcArg = new XYSeries("Calc. resp. phase");
    
    // scaling values, used to set curve values to 0 at 1Hz
    Complex scaleValue = estimatedResponse[oneHzIdx];
    double subtractBy = 10 * Math.log10( scaleValue.abs() );
    double rotateBy = scaleValue.getArgument();
    
    // data to fit poles to; first half of data is magnitudes of resp (dB)
    // second half of data is angles of resp (radians, scaled)
    double[] observedResult = new double[2 * estimatedResponse.length];
    
    for (int i = 0; i < estimatedResponse.length; ++i) {
      
      int argIdx = estimatedResponse.length + i;
      
      Complex estValue = estimatedResponse[i];
      // estValue = estValue.subtract(scaleValue);
      double estValMag = estValue.abs();
      double phi = estValue.getArgument();
      
      if ( Double.isNaN(estValMag) ) {
        observedResult[i] = 0;
        observedResult[argIdx] = 0;
      } else {
        
        observedResult[i] = 10 * Math.log10(estValMag);
        observedResult[i] -= subtractBy;
        
        double argument = phi;
        argument -= rotateBy;
        // argument *= -1;
        observedResult[argIdx] = argument;
        
      }
      
      calcMag.add(freqs[i], observedResult[i]);
      calcArg.add(freqs[i], observedResult[argIdx]);
    }
    
    
    // now to set up a solver for the params -- first, get the input variables
    // complex values are technically two variables, each a double
    // so, let's split up their real and im components and make each one a
    // variable. (we also need to ignore conjugate values, for constraints)
    RealVector initialGuess;
    
    if (lowFreq) {
      initialGuess = lowFreqPolesToVector(fitPoles);
    } else {
      initialGuess =  highFreqPolesToVector(fitPoles);
    }
    
    // now, solve for the response that gets us the best-fit response curve
    
    // RealVector initialGuess = MatrixUtils.createRealVector(responseVariables);
    RealVector obsResVector = MatrixUtils.createRealVector(observedResult);
    
    MultivariateJacobianFunction jacobian = new MultivariateJacobianFunction() {
      
      public Pair<RealVector, RealMatrix> value(final RealVector point) {
        Pair<RealVector, RealMatrix> pair = 
            jacobian(point, fitResponse);
        return pair;
      }
    };
    
    ConvergenceChecker<LeastSquaresProblem.Evaluation> svc = 
        new EvaluationRmsChecker(1.0E-7, 1.0E-7);
    
    LeastSquaresOptimizer optimizer = new LevenbergMarquardtOptimizer().
        withCostRelativeTolerance(1.0E-7).
        withParameterRelativeTolerance(1.0E-7);
    
    LeastSquaresProblem lsp = new LeastSquaresBuilder().
        start(initialGuess).
        target(obsResVector).
        model(jacobian).
        lazyEvaluation(false).
        maxEvaluations(Integer.MAX_VALUE).
        maxIterations(Integer.MAX_VALUE).
        checker(svc).
        build();
    
    System.out.println("lsp built");
    
    // LeastSquaresProblem.Evaluation initEval = lsp.evaluate(initialGuess);
    // initResid = initEval.getRMS() * 100;
    // System.out.println("INITIAL GUESS RESIDUAL: " +  initEval.getRMS() );
    
    LeastSquaresOptimizer.Optimum optimum = optimizer.optimize(lsp);
    // residuals used to determine quality of solution convergence
    LeastSquaresProblem.Evaluation initEval = lsp.evaluate(initialGuess);
    initialResidual = initEval.getCost();
    fitResidual = optimum.getCost();
    
    double[] poleParams = optimum.getPoint().toArray();
    double[] initialValues =
        jacobian.value( initialGuess).getFirst().toArray();
    double[] fitValues = 
        jacobian.value( optimum.getPoint() ).getFirst().toArray();
    
    fitResponse = polesToResp(poleParams, fitResponse, lowFreq);
    fitPoles = fitResponse.getPoles();
    
    // System.out.println(inputPoles);
    // System.out.println(fitPoles);
    
    // System.out.println("FIT PARAMS RESIDUAL: " +  optimum.getRMS() );
    // fitResid = optimum.getRMS() * 100;
    
    // was initially used in plot before values were scaled
    // Complex[] fitRespCurve = fitResponse.applyResponseToInput(freqs);
    
    String name = fitResponse.getName();
    XYSeries initMag = new XYSeries("Initial param (" + name + ") magnitude");
    XYSeries initArg = new XYSeries("Initial param (" + name + ") phase");
    
    XYSeries fitMag = new XYSeries("Fit resp. magnitude");
    XYSeries fitArg = new XYSeries("Fit resp. phase");
    
    for (int i = 0; i < freqs.length; ++i) {
      int argIdx = freqs.length + i;
      initMag.add(freqs[i], initialValues[i]);
      initArg.add(freqs[i], initialValues[argIdx]);
      fitMag.add(freqs[i], fitValues[i]);
      fitArg.add(freqs[i], fitValues[argIdx]);
    }
    
    /*
    for (int i = 0; i < freqs.length; ++i) {
      // int argIdx = freqs.length + i;
      if (freqs[i] != 0) {
        Complex fitRespInteg = fitRespCurve[i].divide(2 * Math.PI * freqs[i]);
        fitMag.add( freqs[i], 10 * Math.log10( fitRespInteg.abs() ) );
        double argument = 
            Math.toDegrees( fitRespCurve[i].getArgument() );
        fitArg.add(freqs[i], argument);
      }
    }
    */
    
    
    // XYSeries fitMag = new XYSeries("Dummy plot a");
    // XYSeries fitArg = new XYSeries("Dummy plot b");
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
    
    System.out.println("Done!");
    
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
  private double[] evaluateResponse(double[] variables, InstrumentResponse ir) {
    
    InstrumentResponse testResp = polesToResp(variables, ir, lowFreq);
    
    Complex[] appliedCurve = testResp.applyResponseToInput(freqs);
    
    // array is magnitudes, then arguments of complex number
    double[] curValue = new double[appliedCurve.length * 2];
    curValue[0] = 0.;
    curValue[appliedCurve.length] = 0.;
    
    Complex scaleBy = appliedCurve[oneHzIdx].divide(2 * Math.PI);
    double magScale = 10 * Math.log10( scaleBy.abs() );
    double argScale = scaleBy.getArgument();
    
    // System.out.println(appliedCurve[0]);
    for (int i = 0; i < appliedCurve.length; ++i) {
      
      int argIdx = appliedCurve.length + i;
      
      if (freqs[i] == 0.) {
        // this would be a divide by 0 error, let's just call the result 0;
        curValue[i] = 0.;
        curValue[argIdx] = 0.;
        continue;
      }
      
      // from acceleration to velocity
      Complex value = appliedCurve[i];
      value = value.divide(2 * Math.PI * freqs[i]);
      
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

        
        // TODO: deal with potential rotation issues with parameters?
        double argument = value.getArgument();
        argument -= argScale;
        // argument /= TAU;
        curValue[argIdx] = argument;
        // curValue[argIdx] = Math.toDegrees(curValue[argIdx]);
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
    // for dealing with complex conjugate issues
    int highFreqStartIdx = 2;
    if ( fitPoles.get(0).getImaginary() == 0. ) {
      highFreqStartIdx = 1;
    }

    if (lowFreq) {
      return fitPoles.subList(0, highFreqStartIdx);
    } else {
      return fitPoles.subList( highFreqStartIdx, inputPoles.size() );
    }
  }

  /**
   * Get the residual value from the solved response parameters
   * @return the residual of the solved-for poles (best-fit response)
   */
  public double getFitResidual() {
    return fitResidual;
  }

  /**
   * Get poles used in input response, for reference against best-fit poles 
   * @return original poles being modified by the solver in calib. processing,
   * as a list
   */
  public List<Complex> getInitialPoles() {
    
    int highFreqStartIdx = 2;
    if ( fitPoles.get(0).getImaginary() == 0. ) {
      highFreqStartIdx = 1;
    }
    
    if (lowFreq) {
      return inputPoles.subList(0, highFreqStartIdx);
    } else {
      return inputPoles.subList( highFreqStartIdx, inputPoles.size() );
    }
  }

  /**
   * Get the residual value of the initial response parameters
   * @return the residual of the initial poles from fed-in response
   */
  public double getInitResidual() {
    return initialResidual;
  }
  
  /**
   * Get name of response file used in the calculation
   */
  public String getResponseName() {
    return responseName;
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
   * @param ir Response whose poles will be modified
   * @return RealVector with evaluation at current response value and 
   * RealMatrix with forward difference of that response (Jacobian)
   */
  private Pair<RealVector, RealMatrix> 
  jacobian(RealVector variables, InstrumentResponse ir) {
    
    int numVars = variables.getDimension();
    
    double[] currentVars = new double[numVars];
    
    for (int i = 0; i < numVars; ++i) {
      currentVars[i] = variables.getEntry(i);
    }
    
    double[] mag = evaluateResponse(currentVars, ir);
    
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
      
      double diffX = changedVars[i] * (1 + DELTA);
      changedVars[i] = diffX;
      
      double[] diffY = 
          evaluateResponse(changedVars, ir);
      
      
      for (int j = 0; j < diffY.length; ++j) {
        jacobian[j][i] = diffY[j] - mag[j];
        jacobian[j][i] /= changedVars[i] - currentVars[i];
      }
      
    }
    
    RealVector result = MatrixUtils.createRealVector(mag);
    RealMatrix jMat = MatrixUtils.createRealMatrix(jacobian);
    
    return new Pair<RealVector, RealMatrix>(result, jMat);
    
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
 
}
