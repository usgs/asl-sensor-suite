package asl.sensor.experiment;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.math3.fitting.leastsquares.LeastSquaresBuilder;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresOptimizer;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresProblem;
import org.apache.commons.math3.fitting.leastsquares.LevenbergMarquardtOptimizer;
import org.apache.commons.math3.fitting.leastsquares.MultivariateJacobianFunction;
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

public class AzimuthExperiment extends Experiment {

  final static double TAU = Math.PI * 2;
  
  double offset = 0.;
  double angle;
  
  public AzimuthExperiment() {
    super();

  }
  
  public void setOffset(double newOffset) {
    offset = newOffset;
  }

  @Override
  protected void backend(final DataStore ds, final boolean freqSpace) {
    
    // assume the first two are the reference and the second two are the test?
    
    // we just need four timeseries, don't actually care about response
    DataBlock refLH1Block = ds.getXthLoadedBlock(1);
    DataBlock refLH2Block = ds.getXthLoadedBlock(2);
    DataBlock testLH1Block = ds.getXthLoadedBlock(3);
    
    List<Number> refLH1 = new ArrayList<Number>( refLH1Block.getData() );
    String refName = refLH1Block.getName();
    List<Number> refLH2 = new ArrayList<Number>( refLH2Block.getData() );
    List<Number> testLH1 = new ArrayList<Number>( testLH1Block.getData() );
    String testName = testLH1Block.getName();
    
    FFTResult.detrend(refLH1);
    FFTResult.detrend(refLH2);
    FFTResult.detrend(testLH1);

    
    refLH1 = TimeSeriesUtils.normalize(refLH1);
    refLH2 = TimeSeriesUtils.normalize(refLH2);
    testLH1 = TimeSeriesUtils.normalize(testLH1);
    
    double sps = TimeSeriesUtils.ONE_HZ_INTERVAL / refLH1Block.getInterval();
    double low = 1./8;
    double high = 1./4;
    
    refLH1 = FFTResult.bandFilter(refLH1, sps, low, high);
    refLH2 = FFTResult.bandFilter(refLH2, sps, low, high);
    testLH1 = FFTResult.bandFilter(testLH1, sps, low, high);
    
    int len = refLH1.size();
    
    double[] refYArr = new double[len];
    double[] refXArr = new double[len];
    double[] testYArr = new double[len];
    
    for (int i = 0; i < len; ++i) {
      refYArr[i] = refLH1.get(i).doubleValue();
      refXArr[i] = refLH2.get(i).doubleValue();
      testYArr[i] = testLH1.get(i).doubleValue();
    }
    
    RealVector refX = MatrixUtils.createRealVector(refXArr);
    RealVector refY = MatrixUtils.createRealVector(refYArr);
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
    
    LeastSquaresOptimizer.Optimum optimumY = optimizer.optimize(findAngleY);
    RealVector angleVector = optimumY.getPoint();
    
    angle = Math.toDegrees( angleVector.getEntry(0) );
    angle = ( (angle % 360) + 360 ) % 360;
    // get the INTERNAL angle of the two components
    if (angle > 180) {
      angle = (360-angle) % 360;
    }
    
    XYSeries ref = new XYSeries(refName);
    ref.add(0 + offset, 0);
    ref.add(0 + offset, 1);
    XYSeries set = new XYSeries(testName + " relative to north");
    set.add(angle, 1);
    set.add(angle, 0);
    XYSeries fromNorth = new XYSeries (testName + " relative to ref. sensor");
    fromNorth.add(angle + offset, 1);
    fromNorth.add(angle + offset, 0);
    
    XYSeriesCollection xysc = (XYSeriesCollection) xySeriesData;
    xysc.addSeries(ref);
    xysc.addSeries(set);
    xysc.addSeries(fromNorth);
    
    /*
    XYSeries diffSrs = new XYSeries("Diff(" + testName + ", " + refName + ")");
    XYSeries diffRotSrs = new XYSeries("Diff(" + testName + ", Rotated Ref.)");
    
    
    RealVector diffLH1 = testY.subtract(refY);
    RealVector diffComponents = jacobian.value(angleVector).getFirst();
    
    for (int i = 0; i < len; ++i) {
      diffSrs.add (timeAtPoint, diffLH1.getEntry(i) );
      diffRotSrs.add( timeAtPoint, diffComponents.getEntry(i) );
      
      timeAtPoint += tick;
    }
    
    XYSeriesCollection xysc = (XYSeriesCollection) xySeriesData;
    
    xysc.addSeries(diffSrs);
    xysc.addSeries(diffRotSrs);
    */
  }
  
  public double getFitAngle() {
    return angle;
  }
  

}