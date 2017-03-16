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
  
  private static double delta = 1 + Double.MIN_VALUE;
  
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
    
    DataStore irHolder = new DataStore();
    irHolder.setResponse(0, fitResponse);
    
    ResponseExperiment respExp = new ResponseExperiment();
    respExp.setFreqSpace(true);
    respExp.setData(irHolder);
    XYSeriesCollection respPlots = (XYSeriesCollection) respExp.getData();
    XYSeries respMag = respPlots.getSeries(0);
    XYSeries respArg = respPlots.getSeries(1);
    
    FFTResult calibFFTResult = FFTResult.singleSidedFFT(calib);
    freqs = calibFFTResult.getFreqs();
    Complex[] calibFFT = calibFFTResult.getFFT();
    Complex[] outputFFT = FFTResult.singleSidedFFT(sensorOut).getFFT();
    
    Complex[] estimatedResponse = new Complex[outputFFT.length];
    
    for (int i = 0; i < estimatedResponse.length; ++i) {
      estimatedResponse[i] = outputFFT[i].divide(calibFFT[i]);
    }
    
    XYSeries calcMag = new XYSeries("Calc. resp. magnitude");
    XYSeries calcArg = new XYSeries("Calc. resp. arg. [phi]");
    
    for (int i = 0; i < estimatedResponse.length; ++i) {
      if (freqs[i] != 0) {
        calcMag.add( freqs[i], 10 * Math.log10( estimatedResponse[i].abs() ) );
        
        double phi = Math.toDegrees( estimatedResponse[i].getArgument() );
        phi = (phi + 360) % 360; // keeps angles positive
        calcArg.add(freqs[i], phi);
      }
    }
    
    // magnitude and result of 
    double[] initResult = new double[estimatedResponse.length * 2];
    for (int i = 0; i < estimatedResponse.length; ++i) {
      int argIdx = i + estimatedResponse.length;
      initResult[i] = /*10 * Math.log10(*/ estimatedResponse[i].abs() /*)*/;
      initResult[argIdx] = estimatedResponse[i].getArgument();
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
    
    // System.out.println("FIT PARAMS RESIDUAL: " +  optimum.getRMS() );
    
    // fitResid = optimum.getRMS() * 100;
    
    double[] magArgFit = evaluateResponse( optimum.getPoint() );
    XYSeries fitMag = new XYSeries("Fit resp. magnitude");
    XYSeries fitArg = new XYSeries("Fit resp. arg. [phi]");
    for (int i = 0; i < freqs.length; ++i) {
      int argIdx = freqs.length + i;
      if (freqs[i] != 0) {
        fitMag.add( freqs[i], 10 * Math.log10(magArgFit[i]) );
        fitArg.add(freqs[i], magArgFit[argIdx]);
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
    double[] magAndAngle = new double[2 * appliedCurve.length];
    for (int i = 0; i < appliedCurve.length; ++i) {
      int argIdx = i + appliedCurve.length;
      magAndAngle[i] = appliedCurve[i].abs();
      magAndAngle[argIdx] = ( appliedCurve[i].getArgument() + 360 ) % 360;
    }
    
    return magAndAngle;
  }
  
  private Pair<RealVector, RealMatrix> jacobian(RealVector variables) {
    
    System.out.println( variables.getDimension() );
    
    int numVars = variables.getDimension();
    
    double[] magAndAngle = evaluateResponse(variables);
    
    double[][] jacobian = new double[magAndAngle.length][numVars];
    // now take the forward difference of each value 
    for (int i = 0; i < numVars; ++i) {
      RealVector dx = variables.copy();
      double increment = dx.getEntry(i);
      increment *= delta;
      dx.setEntry(i, increment);
      
      double[] result = evaluateResponse(dx);
      
      for (int j = 0; j < result.length; ++j) {
        double numerator = result[i] - magAndAngle[i];
        double denominator = Double.MIN_VALUE;
        jacobian[j][i] = numerator / denominator;
      }
    }
    
    RealVector result = MatrixUtils.createRealVector(magAndAngle);
    RealMatrix jMat = MatrixUtils.createRealMatrix(jacobian);
    
    return new Pair<RealVector, RealMatrix>(result, jMat);
    
  }

  @Override
  public boolean hasEnoughData(DataStore ds) {
    return ( ds.blockIsSet(0) && ds.bothComponentsSet(1) );
  }

}
