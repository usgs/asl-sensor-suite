package asl.sensor.experiment;

import java.util.Arrays;

import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealVector;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import asl.sensor.input.DataBlock;
import asl.sensor.input.DataStore;
import asl.sensor.utils.FFTResult;
import asl.sensor.utils.NumericUtils;
import asl.sensor.utils.TimeSeriesUtils;

/**
 * Finds the interior angle between two sensors of unknown orientation using
 * input from two sensors known to be orthogonal and at north and east
 * orientation. The result returns the relative orientation between angles
 * using the azimuth calculation as an intermediate step.
 * (See AzimuthExperiment for details on how the best-fit angles are found)
 * @author akearns
 *
 */
public class OrthogonalExperiment extends Experiment {
  
  /**
   * Return the rotated signal given an angle and orthogonal components
   * @param refX reference signal along the x-axis
   * @param refY reference signal along the y-axis
   * @param point angle (radians) to get as rotated signal
   * @return signal rotated in the direction of the given angle
   */
  public 
  static RealVector value(RealVector refX, RealVector refY, double point) {
    double theta = point % NumericUtils.TAU;
    
    if (theta < 0) {
      theta += NumericUtils.TAU; 
    }
    
    double sinTheta = Math.sin(theta);
    double cosTheta = Math.cos(theta);
    
    RealVector curValue = 
        refX.mapMultiply(sinTheta).add( refY.mapMultiply(cosTheta) );
    
    return curValue;
  }
  
  private double[] diffs;
  
  private double angle;

  public OrthogonalExperiment() {
    super();

  }

  @Override
  protected void backend(final DataStore ds) {
    
    // TODO: refactor using faster access point for azimuth?
    long interval = ds.getXthLoadedBlock(1).getInterval();
    
    // assume the first two are the reference and the second two are the test?
    
    // we just need four timeseries, don't actually care about response
    DataBlock refLH1Block = ds.getXthLoadedBlock(1);
    String refName = refLH1Block.getName();
    dataNames.add(refName);
    DataBlock refLH2Block = ds.getXthLoadedBlock(2);
    dataNames.add( refLH2Block.getName() );
    DataBlock testLH1Block = ds.getXthLoadedBlock(3);
    String testName = testLH1Block.getName();
    dataNames.add(testName);
    DataBlock testLH2Block = ds.getXthLoadedBlock(4);
    dataNames.add( testLH2Block.getName() );
    
    double[] refLH1 = refLH1Block.getData();
    double[] refLH2 = refLH2Block.getData();
    double[] testLH1 = testLH1Block.getData();
    double[] testLH2 = testLH2Block.getData();
    
    TimeSeriesUtils.detrend(refLH1);
    TimeSeriesUtils.detrend(refLH2);
    TimeSeriesUtils.detrend(testLH1);
    TimeSeriesUtils.detrend(testLH2);
    
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
    
    int len = refLH1.length;
    
    double[] refYArr = Arrays.copyOfRange(refLH1, 0, len);
    double[] refXArr = Arrays.copyOfRange(refLH2, 0, len);
    double[] testYArr = Arrays.copyOfRange(testLH1, 0, len);
    
    /*
    for (int i = 0; i < len; ++i) {
      refYArr[i] = refLH1.get(i).doubleValue();
      refXArr[i] = refLH2.get(i).doubleValue();
      testYArr[i] = testLH1.get(i).doubleValue();
    }
    */
    
    // since refLH1 and refLH2 are orthogonal, can use them with azimuth logic
    // to find angle between the other two datasets
    
    DataStore findTestX = new DataStore();
    findTestX.setBlock(0, refLH1Block);
    findTestX.setBlock(1, refLH2Block);
    findTestX.setBlock(2, testLH2Block);
    DataStore findTestY = new DataStore(findTestX);
    findTestY.setBlock(2, testLH1Block);
    
    // TODO: FIX THIS
    AzimuthExperiment azi = new AzimuthExperiment();
    azi.runExperimentOnData(findTestY);
    double angleY = -azi.getFitAngle();
    azi.runExperimentOnData(findTestX);
    double angleX = -azi.getFitAngle();
    
    
    angle = Math.abs(angleY - angleX);
    
    RealVector refX = MatrixUtils.createRealVector(refXArr);
    RealVector refY = MatrixUtils.createRealVector(refYArr);
    RealVector testY = MatrixUtils.createRealVector(testYArr);
    
    angle = ( (angle % 360) + 360 ) % 360;
    // get the INTERNAL angle of the two components
    if (angle > 180) {
      angle = (360-angle) % 360;
    }
    diffs = new double[2];
    diffs[0] = angleY; // north
    diffs[1] = angleX; // east
    diffs[0] = ( (diffs[0] % 360) + 360 ) % 360;
    diffs[1] = ( (diffs[1] % 360) + 360 ) % 360;
    
    // if x-plot chart way above y-plot, plot negative angle
    //if (diffs[1] > diffs[0]) {
    //  diffs[1] -= 360;
    //}
    
    double timeAtPoint = 0.;
    double tick = interval / TimeSeriesUtils.ONE_HZ_INTERVAL;
    
    XYSeries diffSrs = new XYSeries("Diff(" + testName + ", " + refName + ")");
    XYSeries diffRotSrs = new XYSeries("Diff(" + testName + ", Rotated Ref.)");
    
    RealVector diffLH1 = testY.subtract(refY);
    RealVector diffComponents = value(refX, refY, angleY);
    
    for (int i = 0; i < len; ++i) {
      diffSrs.add (timeAtPoint, diffLH1.getEntry(i) );
      diffRotSrs.add( timeAtPoint, diffComponents.getEntry(i) );
      
      timeAtPoint += tick;
    }
    
    XYSeriesCollection xysc = new XYSeriesCollection();
    xysc.addSeries(diffSrs);
    xysc.addSeries(diffRotSrs);
    xySeriesData.add(xysc);
    
  }
  
  @Override
  public int blocksNeeded() {
    return 4;
  }
  
  /**
   * Returns the difference of the best-fit angles for the unknown sensors
   * @return Angle, in degrees
   */
  public double getFitAngle() {
    return angle;
  };

  /**
   * Returns the intermediate result of the calculation, 
   * the azimuth angles of the unknown sensors 
   * @return Array of doubles (size 2), with the north and east azimuth
   * respectively
   */
  public double[] getSolutionParams() {
    return diffs;
  }

  @Override
  public boolean hasEnoughData(DataStore ds) {
    for (int i = 0; i < blocksNeeded(); ++i) {
      if ( !ds.blockIsSet(i) ) {
        return false;
      }
    }
    return true;
  }
  

}
