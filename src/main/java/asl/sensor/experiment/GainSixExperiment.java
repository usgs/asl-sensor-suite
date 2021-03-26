package asl.sensor.experiment;

import static asl.utils.timeseries.TimeSeriesUtils.rotate;
import static asl.utils.timeseries.TimeSeriesUtils.rotateX;

import asl.sensor.input.DataStore;
import asl.utils.timeseries.DataBlock;
import java.util.Arrays;
import java.util.Collections;
import java.util.stream.IntStream;
import org.jfree.data.xy.XYSeries;

/**
 * Augmented version of relative gain experiment that includes
 * calculation of components in all 3 main dimensions; this will
 * rotate north and east components into aligned orientation with
 * the first sensor specified in the data store.
 * This experiment uses the azimuth code to do alignment in these dimensions
 * and then calls a {@link GainExperiment GainExperiment} backend on the data in each dimension.
 *
 * @author akearns - KBRWyle
 */
public class GainSixExperiment extends Experiment {

  /**
   * Number of known space dimensions (as in "3-dimensional")
   */
  static final int DIMENSIONS = 3; // number of known space dimensions
  private final int[] indices;
  // used to store the intermediate result data of each of N,S,V components (in that order)
  private GainExperiment[] componentBackends;
  private double northAngle, eastAngle;
  private int indexOfAngleRefData;
  // indexOfAngleRefData represents which set of data to use as fixed angle reference
  // this can be either 0 (first NEZ set) or 1 (second)
  private int indexOfGainRefData; // as above, but for which data to use as gain reference

  public GainSixExperiment() {
    super();

    indexOfAngleRefData = 0; // default to first set of data
    indexOfGainRefData = 0;

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
        DECIMAL_FORMAT.get().format(getNorthAzimuthDegrees());
    dataStrings[1] = getResultString(1) +
        "\n" + "Estimated azimuth (deg): " +
        DECIMAL_FORMAT.get().format(getEastAzimuthDegrees());
    dataStrings[2] = getResultString(2);
    return dataStrings;
  }

  @Override
  public String[] getInsetStrings() {
    return getDataStrings();
  }

  @Override
  protected void backend(DataStore dataStore) {

    // make sure that xyseriesdata is synchronized
    xySeriesData = Collections.synchronizedList(xySeriesData);

    assignReferenceIndex(); // make sure the expected reference is used
    // while GUI will start with reference index as 0, the server-side code needs to set this
    // before the backend is run

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
    int indexOfRotatingData = (indexOfAngleRefData + 1) % 2;
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
        rotate(northRotate, eastRotate, northAngle);
    stores[0].setBlock(indexOfRotatingData, north2Rotated);
    DataBlock east2Rotated =
        rotateX(northRotate, eastRotate, eastAngle);
    stores[1].setBlock(indexOfRotatingData, east2Rotated);

    // now get the datasets to plug into the datastore
    String[] direction = new String[]{"north", "east", "vertical"};

    //for (int i = 0; i < DIMENSIONS; ++i) {
    IntStream.range(0, DIMENSIONS).parallel().forEach(i -> {
      fireStateChange("Running calculations on " + direction[i] + " components...");
      componentBackends[i].runExperimentOnData(stores[i]);
    });

    // each backend only has one plot's worth of data
    // but is formatted as a list of per-plot data, so we use addAll
    // also get the names of the data going in for use w/ PDF, metadata
    Arrays.stream(componentBackends).parallel().forEachOrdered(
        exp -> xySeriesData.addAll(exp.getData())
    );

    for (GainExperiment componentBackend : componentBackends) {
      dataNames.addAll(componentBackend.getInputNames());
    }

  }

  /**
   * Use the first set of inputs as north and east reference angles (default).
   */
  public void setFirstDataAsAngleReference() {
    indexOfAngleRefData = 0;
  }

  /**
   * Use the second set of inputs as north and east reference angles.
   */
  public void setSecondDataAsAngleReference() {
    indexOfAngleRefData = 1;
  }

  @Override
  public int blocksNeeded() {
    return 6;
  }

  /**
   * Get the rotation angle used to rotate the second input set's east sensor.
   * Ideally this should be close to the value used for the north azimuth.
   *
   * @return Angle of second east sensor (radians) minus 90-degree offset
   * representing angle between north and east sensors; this is the angle sent
   * to the rotation function,
   * {@link asl.utils.timeseries.TimeSeriesUtils#rotateX(DataBlock, DataBlock, double)}
   */
  public double getEastAzimuthDegrees() {
    double eastDegrees = Math.toDegrees(eastAngle);
    while (eastDegrees <= -180) {
      eastDegrees += 360;
    }
    while (eastDegrees > 180) {
      eastDegrees -= 360;
    }
    return eastDegrees;
  }


  /**
   * Get the frequency bounds of the data to be given to charts.
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
   * Get the rotation angle used to rotate the second input set's north sensor.
   *
   * @return Angle of second north sensor (radians)
   */
  public double getNorthAzimuthDegrees() {
    double northDegrees = Math.toDegrees(northAngle);
    while (northDegrees <= -180) {
      northDegrees += 360;
    }
    while (northDegrees > 180) {
      northDegrees -= 360;
    }
    return northDegrees;
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

  /**
   * Select the range of data to be used as the range over which gain statistics are calculated.
   * This is by default the range from 3 to 9 seconds period.
   * This value will be used by each dimension's gain calculations, as in
   * {@link GainExperiment#setRangeForStatistics(double, double) setRangeForStatistics} in the
   * gain subcomponent.
   * @param lowPeriod New low period to get range over
   * @param highPeriod New high period to get range over
   */
  public void setRangeForStatistics(double lowPeriod, double highPeriod) {
    for (GainExperiment component : componentBackends) {
      component.setRangeForStatistics(lowPeriod, highPeriod);
    }
  }

  /**
   * Select which set of data should be used as reference for gain statistics.
   * This is either the first (0) or second(1) set, and will be used as reference for each
   * dimension's gain calculations, as given by
   * {@link asl.sensor.experiment.GainExperiment#setReferenceIndex(int) setReferenceIndex}
   * in the gain subcomponent.
   * @param newIndex Value of dataset to be chosen as reference for all gain estimations
   */
  public void setReferenceIndex(int newIndex) {
    indexOfGainRefData = newIndex;
    assignReferenceIndex();
  }

  private void assignReferenceIndex() {
    for (GainExperiment component : componentBackends) {
      component.setReferenceIndex(indexOfGainRefData);
    }
  }

  /**
   * Returns the set of statistics results of each dimension's experiment using the pre-set
   * frequency values.
   * This produces a 2D array. The first index selects which one of the dimensions to examine data
   * from; the second produces the statistics produced by the gain experiment in that dimension,
   * equivalent to calling {@link asl.sensor.experiment.GainExperiment#getStatsFromFreqs()
   * getStatsFromFreqs}
   * in one of the sub-components of the backend.
   * @return Array holding each dimension's statistical results
   */
  public double[][] getStatistics() {
    double[][] outer = new double[componentBackends.length][];
    for (int i = 0; i < outer.length; ++i) {
      outer[i] = componentBackends[i].getStatsFromFreqs();
    }
    return outer;
  }

  public boolean skipsFIRStages() {
    boolean skipsFIRStages = false;
    for (GainExperiment componentBackend : componentBackends) {
      skipsFIRStages |= componentBackend.skipsFIRStages();
    }
    return skipsFIRStages;
  }
}
