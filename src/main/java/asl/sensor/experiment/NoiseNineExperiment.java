package asl.sensor.experiment;

import static asl.utils.TimeSeriesUtils.rotate;
import static asl.utils.TimeSeriesUtils.rotateX;

import asl.sensor.input.DataStore;
import asl.utils.input.DataBlock;
import java.util.stream.IntStream;
import org.jfree.data.xy.XYSeriesCollection;

/**
 * Enhanced version of the self-noise experiment using 9 inputs. These inputs
 * are the north, east, and vertical sensors from a seismometer setup.
 * While the first sensor is assumed facing in the correct directions and
 * thus orthogonal, the experiment solves for the azimuth for the remaining
 * sensors, which may not face north or have orthogonal orientations.
 * Once the sensors have been properly oriented (the vertical sensors are
 * assumed to be close enough that no rotation is necessary), then each
 * direction has its corresponding self-noise calculation performed. That is,
 * there is a self noise calculation for the north-aligned signals, for the
 * east-aligned signals, and for the vertical signals.
 * See also Ringler, Hutt, "Self-Noise Models of Seismic Instruments",
 * Seismological Research Letters 81 (SSA, Nov. 2010).
 *
 * @author akearns
 */
public class NoiseNineExperiment extends NoiseExperiment {

  private static final int DIMENSIONS = 3;
  private double[] northAngles, eastAngles;
  private int indexOfAngleRefData;

  public NoiseNineExperiment() {
    super();
    // indices are fixed since we need all 9 data points here
    respIndices = new int[]{0, 1, 2, 3, 4, 5, 6, 7, 8};
    indexOfAngleRefData = 0;
  }

  /**
   * Get text representation of angles used to rotate data
   *
   * @return String displaying angles of rotation for the 2nd, 3rd north sensors
   */
  private String getNorthChartString() {
    double[] angles = getNorthAnglesDegrees();

    return anglesToText(angles, "north");
  }

  /**
   * Driver for converting angles to human-readable text, used in plots and reports
   * @param angles Array of angles (north or south)
   * @return Text output of angles with identification of sensor that was rotated
   */
  private String anglesToText(double[] angles, String direction) {
    StringBuilder sb = new StringBuilder();
    for (int i = 0; i < angles.length; ++i) {
      if (i == indexOfAngleRefData) {
        continue;
      }
      sb.append("Angle of rotation of ")
          .append(direction).append(" sensor ").append(i + 1).append(" (deg): ");
      sb.append(DECIMAL_FORMAT.get().format(angles[i]));
      sb.append('\n');
    }
    sb.deleteCharAt(sb.length()-1);
    return sb.toString();
  }

  /**
   * Return the rotation of angles of the north-facing data in degrees, range (-180, 180]
   * @return Array of doubles representing degree rotation
   */
  public double[] getNorthAnglesDegrees() {
    return convertAnglesToDegrees(getNorthAngles());
  }

  /**
   * Return the rotation of angles of the east-facing data in degrees, range (-180, 180]
   * @return Array of doubles representing degree rotation
   */
  public double[] getEastAnglesDegrees() {
    return convertAnglesToDegrees(getEastAngles());
  }

  private double[] convertAnglesToDegrees(double[] angles) {
    for (int i = 0; i < angles.length; ++i) {
      double angle = Math.toDegrees(angles[i]);
      while (angle <= -180) {
        angle += 360;
      }
      while (angle > 180) {
        angle -= 360;
      }
      angles[i] = angle;
    }

    return angles;
  }

  /**
   * Get text representation of angles used to rotate data
   *
   * @return String displaying angles of rotation for the 2nd, 3rd east sensors
   */
  private String getEastChartString() {
    double[] angles = getEastAnglesDegrees();

    return anglesToText(angles, "east");
  }

  @Override
  String[] getDataStrings() {
    return new String[]{getNorthChartString(), getEastChartString()};
  }

  @Override
  public String[] getInsetStrings() {
    String[] result = new String[DIMENSIONS];
    String[] startingPoint = getDataStrings();
    for (int i = 0; i < result.length; ++i) {
      result[i] = "";
      if (i < startingPoint.length) {
        result[i] = startingPoint[i] + '\n';
      }
      result[i] = result[i] + getFormattedDateRange();
    }
    return result;
  }

  @Override
  protected void backend(DataStore dataStore) {

    northAngles = new double[DIMENSIONS];
    eastAngles = new double[DIMENSIONS];

    // NOTE: this may need to change in the event of a test using > 9 inputs
    for (int i = 0; i < 9; ++i) {
      // doing this loop here saves us time and significant lines of code
      dataNames.add(dataStore.getBlock(i).getName());
      dataNames.add(dataStore.getResponse(i).getName());
    }

    DataStore[] stores = new DataStore[DIMENSIONS];

    for (int i = 0; i < DIMENSIONS; ++i) {
      stores[i] = new DataStore();
      for (int j = 0; j < 3; ++j) {
        stores[i].setBlock(j, dataStore.getBlock(i + (j * DIMENSIONS)));
        stores[i].setResponse(j, dataStore.getResponse(i + (j * DIMENSIONS)));
      }
    }

    long interval = dataStore.getBlock(0).getInterval();
    long start = dataStore.getBlock(0).getStartTime();
    long end = dataStore.getBlock(0).getEndTime();

    fireStateChange("Getting data from rotation reference data...");
    DataBlock northRef = stores[0].getBlock(indexOfAngleRefData);
    DataBlock eastRef = stores[1].getBlock(indexOfAngleRefData);
    double[] northReference = northRef.getData();
    double[] eastReference = eastRef.getData();

    // bound here is the number of unknown angles -- one is assumed to be fixed at N & E
    // for (int i = 0; i < DATA_NEEDED; ++i) {
    IntStream.range(0, DATA_NEEDED).parallel().forEach(i -> {
      if (i == indexOfAngleRefData) {
        // unable to rotate the reference -- fix it at 0
        northAngles[i] = 0;
        eastAngles[i] = 0;
        return; // this should break out of the stream loop
      }

      fireStateChange("Getting orientation of data set "
          + (i + 1) + "...");
      DataBlock northRotate = stores[0].getBlock(i);
      DataBlock eastRotate = stores[1].getBlock(i);
      double[] northRotateData = northRotate.getData();
      double[] eastRotateData = eastRotate.getData();

      // angle is set negative because we are finding angle of reference input
      // which is what northRotateData is here
      fireStateChange("Getting north sensor " + (i + 1)
          + " orientation...");
      northAngles[i] = -AzimuthExperiment.getAzimuth(northReference, eastReference,
          northRotateData, interval, start, end);

      fireStateChange("Getting east sensor " + (i + 1)
          + " orientation...");
      // direction north angle should be if north and east truly orthogonal
      // then east component is x component of rotation in that direction
      // i.e., need to correct by 90 degrees to get rotation angle rather than
      // azimuth of east sensor
      // offset by 3Pi/2 is the same as offset Pi/2 (90 degrees) in other
      // rotation direction
      eastAngles[i] = -AzimuthExperiment.getAzimuth(northReference, eastReference,
          eastRotateData, interval, start, end) + (3 * Math.PI / 2);

      fireStateChange("Rotating data " + (i + 1) + "...");
      DataBlock northUnknownRotate =
          rotate(northRotate, eastRotate, northAngles[i]);
      stores[0].setBlock(i, northUnknownRotate);
      DataBlock eastUnknownRotate =
          rotateX(northRotate, eastRotate, eastAngles[i]);
      stores[1].setBlock(i, eastUnknownRotate);
    });

    // set components into N,E,Z directional subcomponents

    // get noise from each axis's data
    NoiseExperiment noiseExp = new NoiseExperiment();
    noiseExp.setFreqSpace(freqSpace);
    String[] directions = new String[]{"north", "east", "vertical"};
    for (int i = 0; i < DIMENSIONS; ++i) {
      String sb = "Calculating " + directions[i] + " noise components...";
      fireStateChange(sb);
      noiseExp.runExperimentOnData(stores[i]);
      XYSeriesCollection xys = noiseExp.getData().get(0);
      xySeriesData.add(xys);
    }

  }

  /**
   * Use the first set of inputs as north and east reference angles (default)
   */
  public void setFirstDataAsAngleReference() {
    indexOfAngleRefData = 0;
  }

  /**
   * Use the second set of inputs as north and east reference angles
   */
  public void setSecondDataAsAngleReference() {
    indexOfAngleRefData = 1;
  }

  /**
   * Use the third set of inputs as north and east reference angles
   */
  public void setThirdDataAsAngleReference() {
    indexOfAngleRefData = 2;
  }


  @Override
  public int blocksNeeded() {
    return 9;
  }

  /**
   * Return array of angles (radian-valued) which east components have been
   * rotated by, starting with the second component (1st east component is
   * assumed to have zero rotation)
   *
   * @return double array representing angles in degrees
   */
  public double[] getEastAngles() {
    return eastAngles;
  }

  /**
   * Return array of angles (radian-valued) which north components have been
   * rotated by, starting with the second component (1st north component is
   * assumed to have zero rotation)
   *
   * @return double array representing angles in degrees
   */
  public double[] getNorthAngles() {
    return northAngles;
  }

  @Override
  public boolean hasEnoughData(DataStore dataStore) {
    for (int i = 0; i < blocksNeeded(); ++i) {
      if (!dataStore.bothComponentsSet(i)) {
        return false;
      }
    }
    return true;
  }
}
