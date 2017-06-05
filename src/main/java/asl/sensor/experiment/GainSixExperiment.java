package asl.sensor.experiment;
 
import org.jfree.data.xy.XYSeries;

import asl.sensor.input.DataBlock;
import asl.sensor.input.DataStore;
import asl.sensor.input.InstrumentResponse;
import asl.sensor.utils.TimeSeriesUtils;

public class GainSixExperiment extends Experiment {

  private static final int DIMS = 3; // number of known space dimensions
  
  // used to store the intermediate result data of each of N,S,V components
  // (in that order)
  private GainExperiment[] componentBackends;
  private int[] indices;
  
  public GainSixExperiment() {
    super();
    
    componentBackends = new GainExperiment[DIMS];
    for (int i = 0; i < componentBackends.length; i++) {
      componentBackends[i] = new GainExperiment();
    }
    
    indices = new int[6]; // TODO: change to get set during backend
    for (int i = 0; i < indices.length; ++i) {
      indices[i] = i;
    }
  }
  

  /**
   * Get octave centered around peak of first (north) components. We assume
   * that all three component gain sets have their peaks at nearly the same
   * points.
   * @param idx Index of north component's data to use as reference
   * @return Upper and lower frequency bound of octave over given peak
   */
  public double[] getOctaveCenteredAtPeak(int idx) {
    // we assume the first sensor is a good enough 
    return componentBackends[0].getOctaveCenteredAtPeak(idx);
  }
  
  /**
   * Get the gain mean and deviation values from the peak frequency of the 
   * given dataset; the range for stats is the octave centered at the peak.
   * @param idx Index of north component's data to use as reference
   * @return Array of form {mean, standard deviation, ref. gain, calc. gain}
   */
  public double[][] getStatsFromPeak(int idx) {
    double[][] result = new double[DIMS][];
    for (int i = 0; i < result.length; ++i) {
      result[i] = componentBackends[i].getStatsFromPeak(idx);
    }
    return result;
  }
  

  /**
   * Get the gain mean and deviation values from a specified peak
   * frequency range.
   * @param idx Index of north component's data to use as reference
   * @param low Low frequency bound of range to get stats over
   * @param high High frequency bound of range to get stats over
   * @return Array of form {mean, standard deviation, ref. gain, calc. gain}
   */
  public double[][] getStatsFromFreqs(int idx, double low, double high) {
    double[][] result = new double[DIMS][];
    for (int i = 0; i < result.length; ++i) {
      result[i] = componentBackends[i].getStatsFromFreqs(idx, low, high);
    }
    return result;
  }
  
  @Override
  public boolean hasEnoughData(DataStore ds) {
    int needed = blocksNeeded();
    for (int i = 0; i < needed; ++i) {
      if ( !ds.bothComponentsSet(i) ) {
        return false;
      }
    }
    return true;
  }


  @Override
  protected void backend(DataStore ds) {
    
    componentBackends = new GainExperiment[DIMS];
    for (int i = 0; i < componentBackends.length; i++) {
      componentBackends[i] = new GainExperiment();
    }
    
    //first get the first set of data (NSV), then second
    DataBlock north1Sensor = ds.getBlock(0);
    InstrumentResponse north1Resp = ds.getResponse(0);
    DataBlock east1Sensor = ds.getBlock(1);
    InstrumentResponse east1Resp = ds.getResponse(1);
    DataBlock vert1Sensor = ds.getBlock(2);
    InstrumentResponse vert1Resp = ds.getResponse(2);
    
    DataBlock north2Sensor = ds.getBlock(3);
    InstrumentResponse north2Resp = ds.getResponse(3);
    DataBlock east2Sensor = ds.getBlock(4);
    InstrumentResponse east2Resp = ds.getResponse(4);
    DataBlock vert2Sensor = ds.getBlock(5);
    InstrumentResponse vert2Resp = ds.getResponse(5);
    
    // now, rotate the first data into the second
    AzimuthExperiment azi = new AzimuthExperiment();
    azi.setSimple(true); // use rough estimate of coherence, no windowing
    DataStore aziStore = new DataStore();
    aziStore.setData(0, north1Sensor);
    aziStore.setData(1, east1Sensor);
    
    // angle should be set negative -- rotate third sensor, not the opposite
    fireStateChange("Getting second north sensor orientation...");
    aziStore.setData(2, north2Sensor);
    azi.runExperimentOnData(aziStore);
    double north2Angle = -azi.getFitAngleRad();

    fireStateChange("Getting second east sensor orientation...");
    aziStore.setData(2, east2Sensor);
    azi.runExperimentOnData(aziStore);
    double east2Angle = -azi.getFitAngleRad() - (Math.PI / 2);
    // need to offset rotation by 90 degrees -- don't want it facing north
    
    // now to rotate the data according to these angles
    fireStateChange("Rotating data...");
    DataBlock north2Rotated =
        TimeSeriesUtils.rotate(north2Sensor, east2Sensor, north2Angle);
    DataBlock east2Rotated = 
        TimeSeriesUtils.rotate(east2Sensor, north2Sensor, east2Angle);
    
    // now get the datasets to plug into the datastore
    DataStore northComponents = new DataStore();
    northComponents.setData(0, north1Sensor);
    northComponents.setResponse(0, north1Resp);
    northComponents.setData(1, north2Rotated);
    northComponents.setResponse(1, north2Resp);
    
    DataStore eastComponents = new DataStore();
    eastComponents.setData(0, east1Sensor);
    eastComponents.setResponse(0, east1Resp);
    eastComponents.setData(1, east2Rotated);
    eastComponents.setResponse(1, east2Resp);
    
    DataStore vertComponents = new DataStore();
    vertComponents.setData(0, vert1Sensor);
    vertComponents.setResponse(0, vert1Resp);
    vertComponents.setData(1, vert2Sensor);
    vertComponents.setResponse(1, vert2Resp);
    
    fireStateChange("Running calculations on north components...");
    componentBackends[0].runExperimentOnData(northComponents);
    fireStateChange("Running calculations on east components...");
    componentBackends[1].runExperimentOnData(eastComponents);
    fireStateChange("Running calculations on vertical components...");
    componentBackends[2].runExperimentOnData(vertComponents);
    
    for (Experiment exp : componentBackends) {
      // each backend only has one plot's worth of data
      // but is formatted as a list of per-plot data, so we use addAll
      xySeriesData.addAll( exp.getData() );
    }
    
  }

  /**
   * Get the frequency bounds of the data to be given to charts
   * @return Array of form {low freq bound, high freq bound}
   */
  public double[] getMinMaxFrequencies() {
    XYSeries xys = xySeriesData.get(0).getSeries(0);
    if ( xySeriesData.get(0).getSeriesKey(0).equals("NLNM") ) {
      xys = xySeriesData.get(0).getSeries(0);
    }
    return new double[]{xys.getMinX(), xys.getMaxX()};
  }
  

  @Override
  public int blocksNeeded() {
    return 6;
  }
  
  @Override
  public int[] listActiveResponseIndices() {
    return indices;
  }
  
  
}
