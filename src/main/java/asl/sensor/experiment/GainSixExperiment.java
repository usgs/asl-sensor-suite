package asl.sensor.experiment;

import asl.sensor.input.DataBlock;
import asl.sensor.input.DataStore;
import asl.sensor.utils.TimeSeriesUtils;
import org.jfree.data.xy.XYSeries;

/**
 * Augmented version of relative gain experiment that includes
 * calculation of components in all 3 main dimensions; this will
 * rotate north and east components into aligned orientation with
 * the first sensor specified in the data store.
 * This experiment uses the azimuth code to do alignment in these dimensions
 * and then calls a gain backend on the data in each dimension
 *
 * @author akearns - KBRWyle
 */
public class GainSixExperiment extends Experiment {

  private static final int DIMENSIONS = 3; // number of known space dimensions
  private final int[] indices;
  // used to store the intermediate result data of each of N,S,V components (in that order)
  private GainExperiment[] componentBackends;
  private double northAngle, eastAngle;
  private int indexOfAngleRefData;
  // indexOfAngleRefData represents which set of data to use as fixed angle reference
  // this can be either 0 (first NEZ set) or 1 (second)

  public GainSixExperiment() {
    super();

    indexOfAngleRefData = 0; // default to first set of data

    componentBackends = new GainExperiment[DIMENSIONS];
    for (int i = 0; i < componentBackends.length; i++) {
      componentBackends[i] = new GainExperiment();
    }

    indices = new int[6];
    for (int i = 0; i < indices.length; ++i) {
      indices[i] = i;
    }
  }

  private String getResultString(int componentIndex) {
    return componentBackends[componentIndex].getResultString();
  }

  @Override
  public String[] getDataStrings() {
    String[] dataStrings = new String[DIMENSIONS];
    dataStrings[0] = getResultString(0) +
        "\n" + "Estimated azimuth (deg): " +
        DECIMAL_FORMAT.get().format(Math.toDegrees(northAngle));
    dataStrings[1] = getResultString(1) +
        "\n" + "Estimated azimuth (deg): " +
        DECIMAL_FORMAT.get().format(Math.toDegrees(eastAngle));
    dataStrings[2] = getResultString(2);
    return dataStrings;
  }

  @Override
  public String[] getInsetStrings() {
    return getDataStrings();
  }

  @Override
  protected void backend(DataStore dataStore) {

    componentBackends = new GainExperiment[DIMENSIONS];
    for (int i = 0; i < componentBackends.length; i++) {
      componentBackends[i] = new GainExperiment();
    }

    long interval = dataStore.getBlock(0).getInterval();
    long start = dataStore.getBlock(0).getStartTime();
    long end = dataStore.getBlock(0).getEndTime();

    DataStore[] stores = new DataStore[DIMENSIONS];

    fireStateChange("Separating data into directional components...");
    for (int i = 0; i < DIMENSIONS; ++i) {
      stores[i] = new DataStore();
      for (int j = 0; j < 2; ++j) {
        stores[i].setBlock(j, dataStore.getBlock(i + (j * DIMENSIONS)));
        stores[i].setResponse(j, dataStore.getResponse(i + (j * DIMENSIONS)));
      }
    }

    // get the data to perform pre-process rotation step
    // note that data store index identifies dimension [N, E, Z]
    // and block index identifies either first or second group of data
    DataBlock northRef = stores[0].getBlock(indexOfAngleRefData);
    DataBlock eastRef = stores[1].getBlock(indexOfAngleRefData);
    double[] northRefSensor = northRef.getData();
    double[] eastRefSensor = eastRef.getData();

    // each dimension has 2 inputs, the rotation reference and to-rotate data
    // if ref is 0, this will be 1; if ref is 1, this will be 0
    int indexOfRotatingData = indexOfAngleRefData + 1 % 2;
    DataBlock northRotate = stores[0].getBlock(indexOfRotatingData);
    DataBlock eastRotate = stores[1].getBlock(indexOfRotatingData);
    double[] northRotateSensor = northRotate.getData();
    double[] east2Sensor = eastRotate.getData();

    // see also the rotation used in the 9-input self noise backend
    fireStateChange("Getting second north sensor orientation...");
    northAngle = -AzimuthExperiment.getAzimuth(northRefSensor, eastRefSensor,
        northRotateSensor, interval, start, end);

    fireStateChange("Getting second east sensor orientation...");
    // direction north angle should be if north and east truly orthogonal
    // then east component is x component of rotation in that direction
    // i.e., need to correct by 90 degrees to get rotation angle rather than
    // azimuth of east sensor
    // offset by 3Pi/2 is the same as offset Pi/2 (90 degrees) in other
    // rotation direction
    eastAngle = -AzimuthExperiment.getAzimuth(northRefSensor, eastRefSensor,
        east2Sensor, interval, start, end) + (3 * Math.PI / 2);

    // now to rotate the data according to these angles
    fireStateChange("Rotating data...");
    DataBlock north2Rotated =
        TimeSeriesUtils.rotate(northRotate, eastRotate, northAngle);
    stores[0].setBlock(indexOfRotatingData, north2Rotated);
    DataBlock east2Rotated =
        TimeSeriesUtils.rotateX(northRotate, eastRotate, eastAngle);
    stores[1].setBlock(indexOfRotatingData, east2Rotated);

    // now get the datasets to plug into the datastore
    String[] direction = new String[]{"north", "east", "vertical"};

    for (int i = 0; i < DIMENSIONS; ++i) {
      fireStateChange("Running calculations on " + direction[i] + " components...");
      componentBackends[i].runExperimentOnData(stores[i]);
    }

    for (Experiment exp : componentBackends) {
      // each backend only has one plot's worth of data
      // but is formatted as a list of per-plot data, so we use addAll
      xySeriesData.addAll(exp.getData());
      // also get the names of the data going in for use w/ PDF, metadata
    }

    for (GainExperiment componentBackend : componentBackends) {
      dataNames.addAll(componentBackend.getInputNames());
    }

  }

  public void setFirstDataAsAngleReference() {
    indexOfAngleRefData = 0;
  }

  public void setSecondDataAsAngleReference() {
    indexOfAngleRefData = 1;
  }

  @Override
  public int blocksNeeded() {
    return 6;
  }

  /**
   * Get the rotation angle used to rotate the second input set's east sensor
   * Ideally this should be close to the value used for the north azimuth
   *
   * @return Angle of second east sensor (radians) minus 90-degree offset
   * representing angle between north and east sensors; this is the angle sent
   * to the rotation function
   * @see TimeSeriesUtils#rotateX
   */
  public double getEastAzimuth() {
    return eastAngle;
  }


  /**
   * Get the frequency bounds of the data to be given to charts
   *
   * @return Array of form {low freq bound, high freq bound}
   */
  public double[] getMinMaxFrequencies() {
    XYSeries xys = xySeriesData.get(0).getSeries(0);
    if (xySeriesData.get(0).getSeriesKey(0).equals("NLNM")) {
      xys = xySeriesData.get(0).getSeries(0);
    }
    return new double[]{xys.getMinX(), xys.getMaxX()};
  }

  /**
   * Get the rotation angle used to rotate the second input set's north sensor
   *
   * @return Angle of second north sensor (radians)
   */
  public double getNorthAzimuth() {
    return northAngle;
  }

  @Override
  public boolean hasEnoughData(DataStore dataStore) {
    int needed = blocksNeeded();
    for (int i = 0; i < needed; ++i) {
      if (!dataStore.bothComponentsSet(i)) {
        return false;
      }
    }
    return true;
  }

  @Override
  public int[] listActiveResponseIndices() {
    return indices;
  }

  public void setRangeForStatistics(double lowPeriod, double highPeriod) {
    for (GainExperiment component : componentBackends) {
      component.setRangeForStatistics(lowPeriod, highPeriod);
    }
  }

  public void setReferenceIndex(int newIndex) {
    for (GainExperiment component : componentBackends) {
      component.setReferenceIndex(newIndex);
    }
  }

}
