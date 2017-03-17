package asl.sensor.experiment;

import java.util.ArrayList;
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
  protected void backend(final DataStore ds) {
    
    // assume the first two are the reference and the second two are the test?
    
    // we just need four timeseries, don't actually care about response
    DataBlock refLH1Block = ds.getXthLoadedBlock(1);
    DataBlock refLH2Block = ds.getXthLoadedBlock(2);
    DataBlock testLH1Block = ds.getXthLoadedBlock(3);
    
    List<Number> testNorth = new ArrayList<Number>( refLH1Block.getData() );
    String northName = refLH1Block.getName();
    List<Number> testEast = new ArrayList<Number>( refLH2Block.getData() );
    String eastName = refLH2Block.getName();
    List<Number> refNorth = new ArrayList<Number>( testLH1Block.getData() );
    String refName = testLH1Block.getName();
    
    FFTResult.detrend(testNorth);
    FFTResult.detrend(testEast);
    FFTResult.detrend(refNorth);

    
    testNorth = TimeSeriesUtils.normalize(testNorth);
    testEast = TimeSeriesUtils.normalize(testEast);
    refNorth = TimeSeriesUtils.normalize(refNorth);
    
    double sps = TimeSeriesUtils.ONE_HZ_INTERVAL / refLH1Block.getInterval();
    double low = 1./8;
    double high = 1./4;
    
    testNorth = FFTResult.bandFilter(testNorth, sps, low, high);
    testEast = FFTResult.bandFilter(testEast, sps, low, high);
    refNorth = FFTResult.bandFilter(refNorth, sps, low, high);
    
    int len = testNorth.size();
    
    double[] testYArr = new double[len];
    double[] testXArr = new double[len];
    double[] refYArr = new double[len];
    
    for (int i = 0; i < len; ++i) {
      testYArr[i] = testNorth.get(i).doubleValue();
      testXArr[i] = testEast.get(i).doubleValue();
      refYArr[i] = refNorth.get(i).doubleValue();
    }
    
    RealVector testX = MatrixUtils.createRealVector(testXArr);
    RealVector testY = MatrixUtils.createRealVector(testYArr);
    RealVector refY = MatrixUtils.createRealVector(refYArr);
    
    MultivariateJacobianFunction jacobian = new MultivariateJacobianFunction() {
      public Pair<RealVector, RealMatrix> value(final RealVector point) {
        double theta = ( point.getEntry(0) % TAU );
        
        if (theta < 0) {
          theta += TAU; 
        }
        
        double sinTheta = Math.sin(theta);
        double cosTheta = Math.cos(theta);
        
        RealVector curValue = 
            testX.mapMultiply(sinTheta).add( testY.mapMultiply(cosTheta) );
        
        // we have only 1 variable, so jacobian is a matrix w/ single column
        RealMatrix jbn = MatrixUtils.createRealMatrix(len, 1);
        RealVector jbnValue = 
            testX.mapMultiply(cosTheta).add( testY.mapMultiply(-sinTheta) );
        jbn.setColumnVector(0, jbnValue);
        
        return new Pair<RealVector, RealMatrix>(curValue, jbn);
      }
    };
    
    LeastSquaresProblem findAngleY = new LeastSquaresBuilder().
        start(new double[] {0}).
        model(jacobian).
        target(refY).
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
    // allows us to have 
    if (angle > 180) {
      angle = (360-angle) % 360;
    }
    
    XYSeries ref = new XYSeries(northName + " rel. to reference");
    ref.add(offset - angle, 0);
    ref.add(offset - angle, 1);
    XYSeries set = new XYSeries(eastName + " rel. to reference");
    set.add(offset - angle + 90, 1);
    set.add(offset - angle + 90, 0);
    XYSeries fromNorth = new XYSeries (refName + " location");
    fromNorth.add(offset, 1);
    fromNorth.add(offset, 0);

    xySeriesData = new XYSeriesCollection();
    xySeriesData.addSeries(ref);
    xySeriesData.addSeries(set);
    xySeriesData.addSeries(fromNorth);
    
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

  @Override
  public boolean hasEnoughData(DataStore ds) {
    for (int i = 0; i < 3; ++i) {
      if ( !ds.blockIsSet(i) ) {
        return false;
      }
    }
    return true;
  }
  

}