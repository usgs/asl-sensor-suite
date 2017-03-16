package asl.sensor.experiment;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

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
import org.apache.commons.math3.util.Pair;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import asl.sensor.input.DataBlock;
import asl.sensor.input.DataStore;
import asl.sensor.utils.FFTResult;
import asl.sensor.utils.TimeSeriesUtils;

public class OrthogonalExperiment extends Experiment {

  final static double TAU = Math.PI * 2;
  
  double[] diffs;
  double angle;
  
  public OrthogonalExperiment() {
    super();

  }

  @Override
  protected void backend(final DataStore ds) {
    
    long interval = ds.getXthLoadedBlock(1).getInterval();
    
    // assume the first two are the reference and the second two are the test?
    
    // we just need four timeseries, don't actually care about response
    DataBlock refLH1Block = ds.getXthLoadedBlock(1);
    DataBlock refLH2Block = ds.getXthLoadedBlock(2);
    DataBlock testLH1Block = ds.getXthLoadedBlock(3);
    DataBlock testLH2Block = ds.getXthLoadedBlock(4);    
    
    List<Number> refLH1 = new ArrayList<Number>( refLH1Block.getData() );
    String refName = refLH1Block.getName();
    List<Number> refLH2 = new ArrayList<Number>( refLH2Block.getData() );
    List<Number> testLH1 = new ArrayList<Number>( testLH1Block.getData() );
    String testName = testLH1Block.getName();
    List<Number> testLH2 = new ArrayList<Number>( testLH2Block.getData() );
    
    FFTResult.detrend(refLH1);
    FFTResult.detrend(refLH2);
    FFTResult.detrend(testLH1);
    FFTResult.detrend(testLH2);
    
    refLH1 = TimeSeriesUtils.normalize(refLH1);
    refLH2 = TimeSeriesUtils.normalize(refLH2);
    testLH1 = TimeSeriesUtils.normalize(testLH1);
    testLH2 = TimeSeriesUtils.normalize(testLH2);
    
    double sps = TimeSeriesUtils.ONE_HZ_INTERVAL / refLH1Block.getInterval();
    double low = 1./8;
    double high = 1./4;
    
    refLH1 = FFTResult.bandFilter(refLH1, sps, low, high);
    refLH2 = FFTResult.bandFilter(refLH2, sps, low, high);
    testLH1 = FFTResult.bandFilter(testLH1, sps, low, high);
    testLH2 = FFTResult.bandFilter(testLH2, sps, low, high);
    
    int len = refLH1.size();
    
    double[] refYArr = new double[len];
    double[] refXArr = new double[len];
    double[] testYArr = new double[len];
    double[] testXArr = new double[len];
    
    for (int i = 0; i < len; ++i) {
      refYArr[i] = refLH1.get(i).doubleValue();
      refXArr[i] = refLH2.get(i).doubleValue();
      testYArr[i] = testLH1.get(i).doubleValue();
      testXArr[i] = testLH2.get(i).doubleValue();
    }
    
    RealVector refX = MatrixUtils.createRealVector(refXArr);
    RealVector refY = MatrixUtils.createRealVector(refYArr);
    RealVector testX = MatrixUtils.createRealVector(testXArr);
    RealVector testY = MatrixUtils.createRealVector(testYArr);
    
    MultivariateJacobianFunction jacobian = new MultivariateJacobianFunction() {
      public Pair<RealVector, RealMatrix> value(final RealVector point) {
        double theta = ( point.getEntry(0) % TAU );
        
        if (theta < 0) {
          theta += TAU; 
        }
        
        double sinTheta = Math.sin(theta);
        double cosTheta = Math.cos(theta);
        
        RealVector curValue = 
            refX.mapMultiply(sinTheta).add( refY.mapMultiply(cosTheta) );
        
        // we have only 1 variable, so jacobian is a matrix w/ single column
        RealMatrix jbn = MatrixUtils.createRealMatrix(len, 1);
        RealVector jbnValue = 
            refX.mapMultiply(cosTheta).add( refY.mapMultiply(-sinTheta) );
        jbn.setColumnVector(0, jbnValue);
        
        return new Pair<RealVector, RealMatrix>(curValue, jbn);
      }
    };
    
    LeastSquaresProblem findAngleX = new LeastSquaresBuilder().
        start(new double[] {Math.PI / 2}).
        model(jacobian).
        target(testX).
        maxEvaluations(Integer.MAX_VALUE).
        maxIterations(Integer.MAX_VALUE).
        lazyEvaluation(false).
        build();
    
    LeastSquaresProblem findAngleY = new LeastSquaresBuilder().
        start(new double[] {0}).
        model(jacobian).
        target(testY).
        maxEvaluations(Integer.MAX_VALUE).
        maxIterations(Integer.MAX_VALUE).
        lazyEvaluation(false).
        build();
        
    LeastSquaresOptimizer optimizer = new LevenbergMarquardtOptimizer().
        withCostRelativeTolerance(1.0E-15).
        withParameterRelativeTolerance(1.0E-15);
    
    LeastSquaresOptimizer.Optimum optimumX = optimizer.optimize(findAngleX);
    LeastSquaresOptimizer.Optimum optimumY = optimizer.optimize(findAngleY);
    
    RealVector solvedX = optimumX.getPoint();
    RealVector solvedY = optimumY.getPoint();
    RealVector angleBetween = solvedY.subtract(solvedX);
    
    angle = Math.toDegrees( angleBetween.getEntry(0) );
    angle = ( (angle % 360) + 360 ) % 360;
    // get the INTERNAL angle of the two components
    if (angle > 180) {
      angle = (360-angle) % 360;
    }
    diffs = new double[2];
    diffs[0] = Math.toDegrees( solvedX.getEntry(0) );
    diffs[0] = ( (diffs[0] % 360) + 360 ) % 360;
    diffs[1] = Math.toDegrees( solvedY.getEntry(0) );
    diffs[1] = ( (diffs[1] % 360) + 360 ) % 360;
    
    // if x-plot chart way above y-plot, plot negative angle
    if (diffs[1] > diffs[0]) {
      diffs[1] -= 360;
    }
    
    System.out.println( Arrays.toString(diffs) );
    
    /*
    // will assume all vectors are trimmed to the same size
    double[] rhsArray = new double[len * 2];
    double[][] lhsArray = new double[len * 2][4];
    
    
    for (int i = 0; i < len; ++i) {
      rhsArray[i] = refLH1.get(i).doubleValue();
      rhsArray[i + len] = refLH2.get(i).doubleValue();
      
      lhsArray[i][0] = testLH1.get(i).doubleValue();
      lhsArray[i][1] = testLH2.get(i).doubleValue();
      lhsArray[i + len][2] = testLH1.get(i).doubleValue();
      lhsArray[i + len][3] = testLH2.get(i).doubleValue();
    }
    
    RealVector rhsVector = MatrixUtils.createRealVector(rhsArray);
    RealMatrix lhsMatrix = MatrixUtils.createRealMatrix(lhsArray);
    
    MultivariateJacobianFunction mjf = new MultivariateJacobianFunction() {
      public Pair<RealVector, RealMatrix> value(final RealVector point) {
        
        RealVector result = lhsMatrix.operate(point);
        return new Pair<RealVector, RealMatrix>(result, lhsMatrix);
      }
    };
    
    LeastSquaresProblem lsp = new LeastSquaresBuilder().
        start(new double[] {1.,0.,0.,1.}).
        model(mjf).
        target(rhsVector).
        maxEvaluations(Integer.MAX_VALUE).
        maxIterations(Integer.MAX_VALUE).
        lazyEvaluation(false).
        build();
    
    LeastSquaresOptimizer.Optimum optimum = optimizer.optimize(lsp);
    
    RealVector fitSolution = optimum.getPoint();
    // double divisor = solution.getEntry(0) + solution.getEntry(1);
        // solution.getMaxValue();
    // solution = solution.mapDivide(divisor);
    
    RealVectorFormat rvf = new RealVectorFormat("[","]",", ");
    System.out.println(rvf.format(fitSolution));
    
    diffs = fitSolution.toArray();
    

    double scale1 = diffs[0];
    double scale2 = diffs[1];
    */
    
    double timeAtPoint = 0.;
    double tick = interval / TimeSeriesUtils.ONE_HZ_INTERVAL;
    
    XYSeries diffSrs = new XYSeries("Diff(" + testName + ", " + refName + ")");
    XYSeries diffRotSrs = new XYSeries("Diff(" + testName + ", Rotated Ref.)");
    
    RealVector diffLH1 = testY.subtract(refY);
    RealVector diffComponents = jacobian.value(solvedY).getFirst();
    
    for (int i = 0; i < len; ++i) {
      diffSrs.add (timeAtPoint, diffLH1.getEntry(i) );
      diffRotSrs.add( timeAtPoint, diffComponents.getEntry(i) );
      
      timeAtPoint += tick;
    }
    
    xySeriesData = new XYSeriesCollection();
    xySeriesData.addSeries(diffSrs);
    xySeriesData.addSeries(diffRotSrs);
    
  }

  public double[] getSolutionParams() {
    return diffs;
  }
  
  public double getFitAngle() {
    return angle;
  }

  @Override
  public boolean hasEnoughData(DataStore ds) {
    for (int i = 0; i < 4; ++i) {
      if ( !ds.blockIsSet(i) ) {
        return false;
      }
    }
    return true;
  }
  

}
