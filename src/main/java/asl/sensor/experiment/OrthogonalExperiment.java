package asl.sensor.experiment;

import asl.sensor.input.DataBlock;
import asl.sensor.input.DataStore;
import asl.sensor.utils.NumericUtils;
import asl.sensor.utils.TimeSeriesUtils;
import java.util.Arrays;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealVector;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

/**
 * Finds the interior angle between two sensors of unknown orientation using
 * input from two sensors known to be orthogonal and at north and east
 * orientation. The result returns the relative orientation between angles
 * using the (full, damped-windowed) azimuth calculation as an intermediate step.
 * (See AzimuthExperiment for details on how the best-fit angles are found)
 *
 * @author akearns
 */
public class OrthogonalExperiment extends Experiment {

  /**
   * Return the rotated signal given an angle and orthogonal components
   *
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
        refX.mapMultiply(sinTheta).add(refY.mapMultiply(cosTheta));

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

    // assume the first two are the test and the second two are the reference?
    // we just need four timeseries, don't actually care about response
    DataBlock refLH1Block = ds.getXthLoadedBlock(1);
    String refName = refLH1Block.getName();
    dataNames.add(refName);
    DataBlock refLH2Block = ds.getXthLoadedBlock(2);
    dataNames.add(refLH2Block.getName());
    DataBlock testLH1Block = ds.getXthLoadedBlock(3);
    String testName = testLH1Block.getName();
    dataNames.add(testName);
    DataBlock testLH2Block = ds.getXthLoadedBlock(4);
    dataNames.add(testLH2Block.getName());

    // this code is used to get the plotted difference between ref + test, ref + rotated test
    double[] refLH1 = refLH1Block.getData();
    double[] refLH2 = refLH2Block.getData();
    double[] testLH1 = testLH1Block.getData();
    double[] testLH2 = testLH2Block.getData();

    refLH1 = TimeSeriesUtils.demean(refLH1);
    refLH2 = TimeSeriesUtils.demean(refLH2);
    testLH1 = TimeSeriesUtils.demean(testLH1);
    testLH2 = TimeSeriesUtils.demean(testLH2);

    System.out.println(refLH1[0] + "," + testLH1[0]);

    refLH1 = TimeSeriesUtils.detrend(refLH1);
    refLH2 = TimeSeriesUtils.detrend(refLH2);
    testLH1 = TimeSeriesUtils.detrend(testLH1);
    testLH2 = TimeSeriesUtils.detrend(testLH2);

    // note that parent class preprocessing should have already downsampled all data to same rate
    // so this just takes it down to 1Hz if it's still above that

    refLH1 = TimeSeriesUtils.decimate(refLH1, interval, TimeSeriesUtils.ONE_HZ_INTERVAL);
    refLH2 = TimeSeriesUtils.decimate(refLH2, interval, TimeSeriesUtils.ONE_HZ_INTERVAL);
    testLH1 = TimeSeriesUtils.decimate(testLH1, interval, TimeSeriesUtils.ONE_HZ_INTERVAL);
    testLH2 = TimeSeriesUtils.decimate(testLH2, interval, TimeSeriesUtils.ONE_HZ_INTERVAL);

    interval = Math.max(interval, TimeSeriesUtils.ONE_HZ_INTERVAL);

    System.out.println(refLH1[0] + "," + testLH1[0]);

    int len = refLH1.length;
    double[] refYArr = Arrays.copyOfRange(refLH1, 0, len);
    double[] refXArr = Arrays.copyOfRange(refLH2, 0, len);
    double[] testYArr = Arrays.copyOfRange(testLH1, 0, len);
    double[] testXArr = Arrays.copyOfRange(testLH2, 0, len);

    AzimuthExperiment azi = new AzimuthExperiment();
    azi.setSimple(false); // set to see if damped window estimates are hurting our results
    fireStateChange("Getting y (north sensor) angle");
    azi.alternateEntryPoint(refYArr, refXArr, testYArr, interval, start, end);
    double angleY = -azi.getFitAngle(); // degrees
    fireStateChange("Getting x (east sensor) angle");
    azi.alternateEntryPoint(refYArr, refXArr, testXArr, interval, start, end);
    double angleX = -azi.getFitAngle();

    angle = Math.abs(angleY - angleX);

    RealVector refX = MatrixUtils.createRealVector(refXArr);
    RealVector refY = MatrixUtils.createRealVector(refYArr);
    RealVector testY = MatrixUtils.createRealVector(testYArr);

    angle = ((angle % 360) + 360) % 360;
    // get the INTERNAL angle of the two components
    if (angle > 180) {
      angle = (360 - angle) % 360;
    }
    diffs = new double[2];
    diffs[0] = angleY; // north
    diffs[1] = angleX; // east
    diffs[0] = ((diffs[0] % 360) + 360) % 360;
    diffs[1] = ((diffs[1] % 360) + 360) % 360;

    double timeAtPoint = 0.;
    double tick = interval / TimeSeriesUtils.ONE_HZ_INTERVAL;

    fireStateChange("Getting plottable data...");

    XYSeries diffSrs = new XYSeries("Diff(" + testName + ", " + refName + ")");
    XYSeries diffRotSrs = new XYSeries("Diff(" + testName + ", Rotated Ref.)");

    RealVector diffLH1 = testY.subtract(refY);
    RealVector diffComponents = testY.subtract(value(refX, refY, angleY));

    System.out.println(refY.getEntry(0) + "," + testY.getEntry(0));

    for (int i = 0; i < len; ++i) {
      diffSrs.add(timeAtPoint, diffLH1.getEntry(i));
      diffRotSrs.add(timeAtPoint, diffComponents.getEntry(i));

      timeAtPoint += tick;
    }

    XYSeriesCollection xysc = new XYSeriesCollection();
    xysc.addSeries(diffSrs);
    xysc.addSeries(diffRotSrs);
    xySeriesData.add(xysc);

    fireStateChange("Done!");

  }

  @Override
  public int blocksNeeded() {
    return 4;
  }

  /**
   * Returns the difference of the best-fit angles for the unknown sensors
   *
   * @return Angle, in degrees
   */
  public double getFitAngle() {
    return angle;
  }

  ;

  /**
   * Returns the intermediate result of the calculation,
   * the azimuth angles of the unknown sensors
   *
   * @return Array of doubles (size 2), with the north and east azimuth
   * respectively
   */
  public double[] getSolutionParams() {
    return diffs;
  }

  @Override
  public boolean hasEnoughData(DataStore ds) {
    for (int i = 0; i < blocksNeeded(); ++i) {
      if (!ds.blockIsSet(i)) {
        return false;
      }
    }
    return true;
  }


}
