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
import asl.sensor.utils.TimeSeriesUtils;

public class RandomizedExperiment extends Experiment {

  private List<Complex> fitPoles;
  private boolean lowFreq; // fit the low- or high-frequency poles?
  private InstrumentResponse fitResponse;
  private double[] freqs;
  
  private static final double delta = 1E-5;
  
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
    
    double nyquistFreq = 
        TimeSeriesUtils.ONE_HZ_INTERVAL / ( sensorOut.getInterval() * 2);
    double proportion = 0.8;
    
    freqs = numeratorPSD.getFreqs();
    // int len = nyquistIdx + 1;
    //int len = numeratorPSD.getFFT().length;
    int len = (int)  ( ( (freqs.length / 2.) + 1. ) * proportion);
    
    double[] resizedFreqs = new double[len];
    for (int i = 0; i < resizedFreqs.length; ++i) {
      resizedFreqs[i] = freqs[i];
    }
    
    freqs = resizedFreqs;
    
    Complex[] appResponse = fitResponse.applyResponseToInput(freqs);
    double maxScaled = Double.NEGATIVE_INFINITY;
    double minScaled = Double.POSITIVE_INFINITY;
    for (int i = 0; i < appResponse.length; ++i) {
      appResponse[i] = appResponse[i].divide( 2 * Math.PI * freqs[i] );
      
    }
    
    maxScaled = appResponse[1].abs();
    minScaled = appResponse[appResponse.length - 1].abs();
    
    Complex[] estimatedResponse = new Complex[len];
    double maxVal = Double.NEGATIVE_INFINITY;
    double minVal = Double.POSITIVE_INFINITY;
    for (int i = 0; i < estimatedResponse.length; ++i) {
      Complex numer = numeratorPSD.getFFT()[i];
      double denom = denominatorPSD.getFFT()[i].getReal();
      estimatedResponse[i] = numer.divide(denom);
      estimatedResponse[i].multiply(2 * Math.PI * freqs[i]);
      if (estimatedResponse[i].abs() > maxVal) {
        maxVal = estimatedResponse[i].abs();
      }
      if (estimatedResponse[i].abs() < minVal) {
        minVal = estimatedResponse[i].abs();
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
      double estValMag = ( estValue.abs() - minVal ) / (maxVal - minVal);
      estValMag *= (maxScaled - minScaled);
      estValMag += minScaled;
     
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
        observedResult[i] = estValMag;

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
    RealVector observedComponents = 
        MatrixUtils.createRealVector(observedResult);
    
    ConvergenceChecker<LeastSquaresProblem.Evaluation> svc = 
        new EvaluationRmsChecker(1E-20, 1E-20);
    
    MultivariateJacobianFunction jacobian = new MultivariateJacobianFunction() {
      int iterator = 0;
      public Pair<RealVector, RealMatrix> value(final RealVector point) {
        ++iterator;
        System.out.print("ITERS: " + iterator + "\n");
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
        withCostRelativeTolerance(1.0E-20).
        withParameterRelativeTolerance(1.0E-20);
    
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
        Complex fitRespInteg = fitRespCurve[i].divide(2 * Math.PI * freqs[i]);
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
    // System.out.println(appliedCurve[0]);
    for (int i = 0; i < appliedCurve.length; ++i) {
      int argIdx = i + appliedCurve.length;
      Complex value = appliedCurve[i];
      if ( value.equals(Complex.NaN) ) {
        // System.out.println("It's NaN");
        mag[i] = 0;
      } else {
        value = value.divide(2 * Math.PI * freqs[i]);
        mag[i] = value.abs();
      }

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
      change = start * (1 + delta);
      dx.setEntry(i, change);
      
      double[] diffY = evaluateResponse(dx);
      
      for (int j = 0; j < diffY.length; ++j) {
        double numerator = diffY[j] - mag[j];
        double denominator = change - start;
        jacobian[j][i] = numerator / denominator;
        if ( Double.isNaN(jacobian[j][i]) ) {
          jacobian[j][i] = 0;
        } else if (Math.abs(jacobian[j][i]) == Double.POSITIVE_INFINITY) {
          System.out.println("This would mean the change is zero?!");
        }
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
