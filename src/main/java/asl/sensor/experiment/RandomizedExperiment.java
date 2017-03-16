package asl.sensor.experiment;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.fitting.leastsquares.EvaluationRmsChecker;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresBuilder;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresOptimizer;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresProblem;
import org.apache.commons.math3.fitting.leastsquares.LevenbergMarquardtOptimizer;
import org.apache.commons.math3.fitting.leastsquares.MultivariateJacobianFunction;
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

public class RandomizedExperiment extends Experiment {

  private List<Complex> fitPoles;
  private boolean lowFreq; // fit the low- or high-frequency poles?
  private InstrumentResponse fitResponse;
  private double[] freqs;
  
  private static double delta = 1 + 0.1;
  
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
    fitResponse = ds.getResponse(sensorOutIndex);
    
    fitPoles = new ArrayList<Complex>( fitResponse.getPoles() );
    
    // first, get the plots of the response
    DataStore irHolder = new DataStore();
    irHolder.setResponse(0, fitResponse);
    
    ResponseExperiment respExp = new ResponseExperiment();
    respExp.setFreqSpace(true);
    respExp.setData(irHolder);
    XYSeriesCollection respPlots = (XYSeriesCollection) respExp.getData();
    XYSeries respMag = respPlots.getSeries(0);
    XYSeries respArg = respPlots.getSeries(1);
    
    // second, get the plots of the calculated response from deconvolution
    // PSD(out, in) / PSD(in, in) gives us PSD(out) / PSD(in) while removing
    // imaginary terms from the denominator due to multiplication with the
    // complex conjugate
    // PSD(out) / PSD(in) is the response curve (i.e., deconvolution)
    
    FFTResult numeratorPSD = FFTResult.spectralCalc(sensorOut, calib);
    FFTResult denominatorPSD = FFTResult.spectralCalc(calib, calib);
    
    freqs = numeratorPSD.getFreqs();
    int len = numeratorPSD.getFFT().length;
    
    Complex[] estimatedResponse = new Complex[len];
    
    for (int i = 0; i < estimatedResponse.length; ++i) {
      Complex numer = numeratorPSD.getFFT()[i];
      double denom = denominatorPSD.getFFT()[i].getReal();
      estimatedResponse[i] = numer.divide(denom);
    }
    
    XYSeries calcMag = new XYSeries("Calc. resp. magnitude");
    XYSeries calcArg = new XYSeries("Calc. resp. arg. [phi]");
    
    double[] initResult = new double[estimatedResponse.length];
    
    for (int i = 0; i < estimatedResponse.length; ++i) {
      Complex estValue = estimatedResponse[i].multiply(2 * Math.PI * freqs[i]);
     
      double phi = Math.toDegrees( estValue.getArgument() );
      phi = (phi + 360) % 360; // keeps angles positive
      
      if (freqs[i] != 0) {
        calcMag.add( freqs[i], 10 * Math.log10( estValue.abs() ) );
        
        
        calcArg.add(freqs[i], phi);
      }
      
      // int argIdx = i + estimatedResponse.length;
      initResult[i] = estValue.abs();
      // initResult[argIdx] = phi;
    }
    
    // now to set up a solver for the params
    double[] responseVariables;
    if (lowFreq) {
      responseVariables = new double[4];
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
    RealVector observedComponents = MatrixUtils.createRealVector(initResult);
    
    ConvergenceChecker<LeastSquaresProblem.Evaluation> svc = 
        new EvaluationRmsChecker(1E-50, 1E-50);
        
    MultivariateJacobianFunction jacobian = new MultivariateJacobianFunction() {
      public Pair<RealVector, RealMatrix> value(final RealVector point) {
        return jacobian(point);
      }
    };
    
    LeastSquaresProblem lsp = new LeastSquaresBuilder().
        start(initialGuess).
        target(observedComponents).
        model(jacobian).
        lazyEvaluation(false).
        maxEvaluations(Integer.MAX_VALUE).
        maxIterations(Integer.MAX_VALUE).
        checker(svc).
        build();
    
    // LeastSquaresProblem.Evaluation initEval = lsp.evaluate(initialGuess);
    // initResid = initEval.getRMS() * 100;
    // System.out.println("INITIAL GUESS RESIDUAL: " +  initEval.getRMS() );
    
    LeastSquaresOptimizer optimizer = new LevenbergMarquardtOptimizer().
        withCostRelativeTolerance(1.0E-15).
        withParameterRelativeTolerance(1.0E-15);
    
    LeastSquaresOptimizer.Optimum optimum = optimizer.optimize(lsp);
    
    double[] poleParams = optimum.getPoint().toArray();
    List<Complex> poles = new ArrayList<Complex>(fitPoles);
    List<Complex> builtPoles = new ArrayList<Complex>();
    // time to populate the poles
    for (int i = 0; i < poleParams.length; i += 2) {
      Complex c = new Complex(poleParams[i], poleParams[i + 1]);
      System.out.println(c);
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
        Complex fitRespInteg = fitRespCurve[i].multiply(  
           2 * Math.PI * freqs[i]);
        fitMag.add( freqs[i], 10 * Math.log10( fitRespInteg.abs() ) );
        double argument = Math.toDegrees( fitRespCurve[i].getArgument() );
        fitArg.add(freqs[i], ( argument + 360 ) % 360 );
      }
    }
    
    xySeriesData.addSeries(respMag);
    xySeriesData.addSeries(calcMag);
    xySeriesData.addSeries(fitMag);
    
    xySeriesData.addSeries(respArg);
    xySeriesData.addSeries(calcArg);
    xySeriesData.addSeries(fitArg);
    
    System.out.println("Done!");
    
  }
  
  private double[] evaluateResponse(RealVector variables) {
    
    int numVars = variables.getDimension();
    
    InstrumentResponse testResp = new InstrumentResponse(fitResponse);
    
    List<Complex> poleList = new ArrayList<Complex>( fitResponse.getPoles() );
    List<Complex> builtPoles = new ArrayList<Complex>();
    
    for (int i = 0; i < numVars; i += 2) {
      Complex c = new Complex(variables.getEntry(i), variables.getEntry(i+1));
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
    
    // get the result for the input value
    testResp.setPoles(poleList);
    Complex[] appliedCurve = testResp.applyResponseToInput(freqs);
    
    // array is magnitudes, then arguments of complex number
    double[] mag = new double[appliedCurve.length];
    for (int i = 0; i < appliedCurve.length; ++i) {
      int argIdx = i + appliedCurve.length;
      Complex value = appliedCurve[i].multiply(2 * Math.PI * freqs[i]);
      mag[i] = value.abs();
      /*
      double argument = Math.toDegrees( value.getArgument() );
      magAndAngle[argIdx] = argument; // taking (% 360) bad for solver
      */
    }
    
    return mag;
  }
  
  private Pair<RealVector, RealMatrix> jacobian(RealVector variables) {
    
    int numVars = variables.getDimension();
    
    double[] mag = evaluateResponse(variables);
    
    double[][] jacobian = new double[mag.length][numVars];
    double start, change;
    // now take the forward difference of each value 
    for (int i = 0; i < numVars; ++i) {
      RealVector dx = variables.copy();
      start = variables.getEntry(i);
      change = start * delta;
      dx.setEntry(i, change);
      
      double[] result = evaluateResponse(dx);
      
      for (int j = 0; j < result.length; ++j) {
        double numerator = result[i] - mag[i];
        double denominator = change - start;
        jacobian[j][i] = numerator / denominator;
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
