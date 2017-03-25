package asl.sensor.experiment;

import java.util.ArrayList;
import java.util.Arrays;
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

import asl.sensor.input.DataBlock;
import asl.sensor.input.DataStore;
import asl.sensor.input.InstrumentResponse;
import asl.sensor.utils.FFTResult;

public class RandomizedExperiment extends Experiment {

  private List<Complex> fitPoles;
  private boolean lowFreq; // fit the low- or high-frequency poles?
  private InstrumentResponse fitResponse;
  private double[] freqs;
  
  private static final double DELTA = 1E-7;
  
  public RandomizedExperiment() {
    super();
    lowFreq = false;
  }
  
  public void setLowFreq(boolean lowFreq) {
    this.lowFreq = lowFreq;
  }
  
  @Override
  protected void backend(DataStore ds) {
    // TODO Auto-generated method stub
    
    // construct response plot
    DataBlock calib = ds.getXthLoadedBlock(1);
    int sensorOutIndex = ds.getXthFullyLoadedIndex(1);
    if ( ds.getBlock(sensorOutIndex).getName().equals( calib.getName() ) ) {
      sensorOutIndex = ds.getXthFullyLoadedIndex(2);
    }
    
    DataBlock sensorOut = ds.getBlock(sensorOutIndex);
    fitResponse = new InstrumentResponse( ds.getResponse(sensorOutIndex) );
    
    fitPoles = new ArrayList<Complex>( fitResponse.getPoles() );
    
    // get the plots of the calculated response from deconvolution
    // PSD(out, in) / PSD(in, in) gives us PSD(out) / PSD(in) while removing
    // imaginary terms from the denominator due to multiplication with the
    // complex conjugate
    // PSD(out) / PSD(in) is the response curve (i.e., deconvolution)
    
    // also, use those frequencies to get the applied response to input
    
    FFTResult numeratorPSD = FFTResult.spectralCalc(sensorOut, calib);
    FFTResult denominatorPSD = FFTResult.spectralCalc(calib, calib);
    
    double proportion = 0.8; // get up to .8 of Nyquist freq due to noise
    double minFreq = .2; // lower bound of .2 Hz (5s period) due to noise
    
    freqs = numeratorPSD.getFreqs();    

    int len = (int) (freqs.length * proportion); // index of peak freq to check
    
    // trim down the frequency array to the specified range
    // we use a list because the lower bound is not fixed by index
    int oneHzIdx = 0;
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
    
    Collections.sort(freqList); // probably not necessay, but for peace of mind
    
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
    
    Complex[] appResponse = fitResponse.applyResponseToInput(freqs);
    for (int i = 0; i < appResponse.length; ++i) {
      appResponse[i] = appResponse[i].divide( 2 * Math.PI * freqs[i] );
      
    }
    
    double scaleBy = appResponse[oneHzIdx].abs(); // value at 1Hz, to normalize
    
    Complex[] estimatedResponse = new Complex[len];
    for (int i = 0; i < estimatedResponse.length; ++i) {
      Complex numer = numeratorPSDVals[i];
      Complex denom = denominatorPSDVals[i];
      estimatedResponse[i] = numer.divide(denom);
      estimatedResponse[i] = 
          estimatedResponse[i].multiply(2 * Math.PI * freqs[i]);
      
      if (i == oneHzIdx) {
        double scaleDenom = estimatedResponse[oneHzIdx].abs();
        scaleBy /= scaleDenom;
      }
      
    }
    
    // next, normalize estimated response
    
    XYSeries respMag = new XYSeries(fitResponse.getName() + " magnitude");
    XYSeries respArg = new XYSeries(fitResponse.getName() + " arg. [phi]");
    XYSeries calcMag = new XYSeries("Calc. resp. magnitude");
    XYSeries calcArg = new XYSeries("Calc. resp. arg. [phi]");
    
    double[] observedResult = new double[estimatedResponse.length];
    
    for (int i = 0; i < estimatedResponse.length; ++i) {
      Complex estValue = estimatedResponse[i];
      double estValMag = estValue.abs() * scaleBy;
     
      double phi = Math.toDegrees( estValue.getArgument() );
      phi = (phi + 360) % 360; // keeps angles positive
      
      double respPhi = Math.toDegrees( appResponse[i].getArgument() );
      respPhi = (respPhi + 360) % 360;
      
      if (freqs[i] != 0) {
        respMag.add( freqs[i], 10 * Math.log10( appResponse[i].abs() ) );
        calcMag.add( freqs[i], 10 * Math.log10(estValMag) );
        respArg.add(freqs[i], respPhi);
        calcArg.add(freqs[i], phi);
      }
      
      // int argIdx = i + estimatedResponse.length;
      if ( Double.isNaN(estValMag) ) {
        observedResult[i] = 0;
      } else {
        observedResult[i] = 10 * Math.log10(estValMag);

      }
      // initResult[argIdx] = phi;
    }
    
    // now to set up a solver for the params
    double[] responseVariables;
    if (lowFreq) {
      responseVariables = new double[2 * 2];
      for (int i = 0; i < responseVariables.length; i += 2) {
        int realIdx = i;
        int imagIdx = realIdx + 1;
        int poleIdx = i / 2;
        responseVariables[realIdx] = fitPoles.get(poleIdx).getReal();
        responseVariables[imagIdx] = fitPoles.get(poleIdx).getImaginary();
      }
    } else {
      responseVariables = new double[( fitPoles.size() - 2 ) * 2];
      for (int i = 0; i < responseVariables.length; i += 2) {
        int realIdx = i;
        int imagIdx = realIdx + 1;
        int poleIdx = (i / 2) + 2; // don't include first two poles
        responseVariables[realIdx] = fitPoles.get(poleIdx).getReal();
        responseVariables[imagIdx] = fitPoles.get(poleIdx).getImaginary();
      }
    }
    
    
    // now, solve for the response that gets us the best-fit response curve
    
    RealVector initialGuess = MatrixUtils.createRealVector(responseVariables);
    RealVector obsResVector = MatrixUtils.createRealVector(observedResult);
    
    MultivariateJacobianFunction jacobian = new MultivariateJacobianFunction() {
      double[] observed = observedResult; // used for difference evaluation
      public Pair<RealVector, RealMatrix> value(final RealVector point) {
        Pair<RealVector, RealMatrix> pair = jacobian(point, observed);
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
    
    // LeastSquaresProblem.Evaluation initEval = lsp.evaluate(initialGuess);
    // initResid = initEval.getRMS() * 100;
    // System.out.println("INITIAL GUESS RESIDUAL: " +  initEval.getRMS() );
    
    LeastSquaresOptimizer.Optimum optimum = optimizer.optimize(lsp);
    
    double[] poleParams = optimum.getPoint().toArray();
    List<Complex> poles = new ArrayList<Complex>(fitPoles);
    List<Complex> builtPoles = new ArrayList<Complex>();
    // time to populate the poles
    for (int i = 0; i < poleParams.length; i += 2) {
      Complex c = new Complex(poleParams[i], poleParams[i + 1]);
      builtPoles.add(c);
    }
    
    System.out.println(fitPoles);
    
    if (lowFreq) {
      for (int i = 0; i < 2; ++i) {
        poles.set( i, builtPoles.get(i) );
      }
    } else {
      for (int i = 0; i < builtPoles.size() - 2; ++i) {
        poles.set( i + 2, builtPoles.get(i) );
      }
    }
    
    fitResponse.setPoles(poles);
    fitPoles = poles;
    
    System.out.println(fitPoles);
    
    // System.out.println("FIT PARAMS RESIDUAL: " +  optimum.getRMS() );
    
    // fitResid = optimum.getRMS() * 100;
    
    Complex[] fitRespCurve = fitResponse.applyResponseToInput(freqs);
    XYSeries fitMag = new XYSeries("Fit resp. magnitude");
    XYSeries fitArg = new XYSeries("Fit resp. arg. [phi]");
    for (int i = 0; i < freqs.length; ++i) {
      // int argIdx = freqs.length + i;
      if (freqs[i] != 0) {
        Complex fitRespInteg = fitRespCurve[i].divide(2 * Math.PI * freqs[i]);
        fitMag.add( freqs[i], 10 * Math.log10( fitRespInteg.abs() ) );
        double argument = Math.toDegrees( fitRespCurve[i].getArgument() );
        fitArg.add(freqs[i], ( argument + 360 ) % 360 );
      }
    }
    
    
    // XYSeries fitMag = new XYSeries("Dummy plot a");
    // XYSeries fitArg = new XYSeries("Dummy plot b");
    
    xySeriesData.addSeries(respMag);
    xySeriesData.addSeries(calcMag);
    xySeriesData.addSeries(fitMag);
    
    xySeriesData.addSeries(respArg);
    xySeriesData.addSeries(calcArg);
    xySeriesData.addSeries(fitArg);
    
    System.out.println("Done!");
    
  }
  
  private InstrumentResponse polesToResponse(double[] variables) {
    
    int numVars = variables.length;
    InstrumentResponse testResp = new InstrumentResponse(fitResponse);
    
    List<Complex> poleList = new ArrayList<Complex>( testResp.getPoles() );
    List<Complex> builtPoles = new ArrayList<Complex>();
    
    for (int i = 0; i < numVars; i += 2) {
      // System.out.println(variables[i]+","+variables[i+1]);
      Complex c = new Complex(variables[i], variables[i+1]);
      builtPoles.add(c);
    }
    
    if (lowFreq) {
      for (int i = 0; i < 2; ++i) {
        poleList.set( i, builtPoles.get(i) );
      }
    } else {
      for (int i = 0; i < builtPoles.size() - 2; ++i) {
        poleList.set( i + 2, builtPoles.get(i) );
      }
    }
    
    // System.out.println(poleList);
    
    // get the result for the input value
    testResp.setPoles(poleList);
    return testResp;
  }
  
  private double[] evaluateResponse(double[] variables, double[] observed) {
    
    int numVars = variables.length;
    
    InstrumentResponse testResp = polesToResponse(variables);
    
    Complex[] appliedCurve = testResp.applyResponseToInput(freqs);
    
    // array is magnitudes, then arguments of complex number
    double[] mag = new double[appliedCurve.length];
    mag[0] = 0.;
    // System.out.println(appliedCurve[0]);
    for (int i = 0; i < appliedCurve.length; ++i) {
      if (freqs[i] == 0.) {
        continue; // this would cause a div by 0 error, don't use it 
      }
      Complex value = appliedCurve[i];
      value = value.divide(2 * Math.PI * freqs[i]);
      if ( value.equals(Complex.NaN) ) {
        System.out.println("It's NaN: " + i+"; freq: "+freqs[i]);
        mag[0] += 0;
      } else {
        // System.out.println(value);
        double temp = 10 * Math.log10( value.abs() );
        mag[i] = temp;
      }

      /*
      double argument = Math.toDegrees( value.getArgument() );
      magAndAngle[argIdx] = argument; // taking (% 360) bad for solver
      */
    }
    
    return mag;
  }
  
  private Pair<RealVector, RealMatrix> 
  jacobian(RealVector variables, double[] observed) {
    
    int numVars = variables.getDimension();
    
    double[] currentVars = new double[numVars];
    
    for (int i = 0; i < numVars; ++i) {
      currentVars[i] = variables.getEntry(i);
    }
    
    double[] mag = evaluateResponse(currentVars, observed);
    
    double[][] jacobian = new double[mag.length][numVars];
    // now take the forward difference of each value 
    for (int i = 0; i < numVars; ++i) {
      
      double[] changedVars = new double[currentVars.length];
      for (int j = 0; j < currentVars.length; ++j) {
        changedVars[j] = currentVars[j];
      }
      // RealVector dx = variables.copy();
      // dx.setEntry(i, change);
      
      double diffX = Math.pow(changedVars[i], 9);
      changedVars[i] = diffX;
      
      System.out.println(changedVars[i]);
      
      double[] diffY = evaluateResponse(changedVars, observed);
      
      /*
      System.out.println( Arrays.toString(mag) + "," + Arrays.toString(diffY) );
      */
      
      for (int j = 0; j < diffY.length; ++j) {
        jacobian[j][i] = diffY[j] - mag[j];
        if (jacobian[j][i] != 0.) {
          System.out.println("It's NOT zero!!");
        }
        jacobian[j][i] /= changedVars[i] - currentVars[i];
      }
      
    }
    
    RealVector result = MatrixUtils.createRealVector(mag);
    RealMatrix jMat = MatrixUtils.createRealMatrix(jacobian);
    
    return new Pair<RealVector, RealMatrix>(result, jMat);
    
  }

  @Override
  public boolean hasEnoughData(DataStore ds) {
    return ( ds.blockIsSet(0) && ds.bothComponentsSet(1) );
  }

}
