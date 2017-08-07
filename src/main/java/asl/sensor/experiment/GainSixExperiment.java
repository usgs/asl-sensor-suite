package asl.sensor.experiment;
 
import java.util.List;

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
  private double north2Angle, east2Angle;
  
  public GainSixExperiment() {
    super();
    
    componentBackends = new GainExperiment[DIMS];
    for (int i = 0; i < componentBackends.length; i++) {
      componentBackends[i] = new GainExperiment();
    }
    
    indices = new int[6]; // TODO: change to get set during backend?
    for (int i = 0; i < indices.length; ++i) {
      indices[i] = i;
    }
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
    DataStore aziStore = new DataStore();
    aziStore.setData(0, north1Sensor);
    aziStore.setData(1, east1Sensor);
    
    // see also the rotation used in the 9-input self noise backend
    fireStateChange("Getting second north sensor orientation...");
    aziStore.setData(2, north2Sensor);
    north2Angle = -getAzimuth(aziStore);

    fireStateChange("Getting second east sensor orientation...");
    aziStore.setData(2, east2Sensor);
    // direction north angle should be if north and east truly orthogonal
    // then east component is x component of rotation in that direction
    // i.e., need to correct by 90 degrees to get rotation angle rather than
    // azimuth of east sensor
    // offset by 3Pi/2 is the same as offset Pi/2 (90 degrees) in other 
    // rotation direction
    east2Angle = -getAzimuth(aziStore) + (3 * Math.PI / 2);
    
    // now to rotate the data according to these angles
    fireStateChange("Rotating data...");
    DataBlock north2Rotated =
        TimeSeriesUtils.rotate(north2Sensor, east2Sensor, north2Angle);
    DataBlock east2Rotated = 
        TimeSeriesUtils.rotateX(north2Sensor, east2Sensor, east2Angle);
    
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
      // also get the names of the data going in for use w/ PDF, metadata
    }
    
    List<String> northNames = componentBackends[0].getInputNames();
    List<String> eastNames = componentBackends[1].getInputNames();
    List<String> vertNames = componentBackends[2].getInputNames();
    
    for (int i = 0; i < northNames.size(); i += 2) {
      dataNames.add( northNames.get(i) );
      dataNames.add( northNames.get(i + 1) );
      dataNames.add( eastNames.get(i) );
      dataNames.add( eastNames.get(i + 1) );
      dataNames.add( vertNames.get(i) );
      dataNames.add( vertNames.get(i + 1) );
    }
    
  }
  
  @Override
  public int blocksNeeded() {
    return 6;
  }

  /**
   * Get the rotation angle used to rotate the second input set's east sensor
   * Ideally this should be close to the value used for the north azimuth
   * @return Angle of second east sensor (radians) minus 90-degree offset 
   * representing angle between north and east sensors; this is the angle sent
   * to the rotation function
   * @see TimeSeriesUtils#rotateX
   */
  public double getEastAzimuth() {
    return east2Angle;
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
  

  /**
   * Get the rotation angle used to rotate the second input set's north sensor
   * @return Angle of second north sensor (radians)
   */
  public double getNorthAzimuth() {
    return north2Angle;
  }
  
  /**
   * Get octave centered around peak of vertical components. We assume
   * that all three component gain sets have their peaks at nearly the same
   * points. The vertical sensors are most likely to be useful, as their
   * data does not need to be rotated as with the north and east sensors.
   * @param idx Index of vertical component data to use as reference
   * @return Upper and lower frequency bound of octave over given peak
   */
  public double[] getOctaveCenteredAtPeak(int idx) {
    // we assume the vertical sensor, as it is not rotated, is
    return componentBackends[2].getOctaveCenteredAtPeak(idx);
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
    // vertical component requires no rotation
    
    for (int i = 0; i < DIMS; ++i) {
      result[i] = componentBackends[i].getStatsFromFreqs(idx, low, high);
    }
    
    /*
    // commented out if we don't actually need to revert from the calc values
    double meanAngle = (north2Angle + east2Angle) / 2;
    // do we need to rotate the north and east components?
    // get the values we need for doing rotation
    double northRef = result[0][2];
    double northCalc = result[0][3];
    double eastRef = result[1][2];
    double eastCalc = result[1][3];
    if ( idx == 0 ) {
      // the fixed data is being used as reference, un-rotate calc'd results
      double calcNish = (northCalc - eastCalc) / ( 2 * Math.cos(meanAngle) );
      double calcEish = (northCalc + eastCalc) / ( 2 * Math.sin(meanAngle) );
      result[0][3] = calcNish;
      result[1][3] = calcEish;
    } else {
      // the rotated data is being used as reference, un-rotate the reference
      double refNish = (northRef - eastRef) / ( 2 * Math.cos(meanAngle) );
      double refEish = (northRef + eastRef) / ( 2 * Math.sin(meanAngle) );
      result[0][2] = refNish;
      result[1][2] = refEish;
    }
    */
    return result;
  }

  /**
   * Get the gain mean and deviation values from the peak frequency of the 
   * given dataset; the range for stats is the octave centered at the peak.
   * @param idx Index of north component's data to use as reference
   * @return Array of form {mean, standard deviation, ref. gain, calc. gain}
   */
  public double[][] getStatsFromPeak(int idx) {
    
    double[] octave = getOctaveCenteredAtPeak(idx);
    return getStatsFromFreqs(idx, octave[0], octave[1]);
    
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
  public int[] listActiveResponseIndices() {
    return indices;
  }
  
  private double getAzimuth(DataStore ds) {
    AzimuthExperiment azi = new AzimuthExperiment();
    azi.setSimple(true); // do the faster angle calculation
    azi.runExperimentOnData(ds);
    return azi.getFitAngleRad();
  }
}
