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
    
    long interval = ds.getBlock(0).getInterval();
    long start = ds.getBlock(0).getStartTime();
    long end = ds.getBlock(0).getEndTime();
    
    DataStore[] stores = new DataStore[DIMS];
    
    fireStateChange("Separating data into directional components...");
    for (int i = 0; i < DIMS; ++i) {
      stores[i] = new DataStore();
      for (int j = 0; j < 2; ++j) {
        stores[i].setBlock( j, ds.getBlock( i + (j * DIMS) ) );
        stores[i].setResponse( j, ds.getResponse( i + (j * DIMS) ) );
      }
    }
    
    //first get the first set of data (NSV), then second
    double[] north1Sensor = ds.getBlock(0).getData();
    double[] east1Sensor = ds.getBlock(1).getData();
    
    double[] north2Sensor = ds.getBlock(3).getData();
    double[] east2Sensor = ds.getBlock(4).getData();
    
    // see also the rotation used in the 9-input self noise backend
    fireStateChange("Getting second north sensor orientation...");
    north2Angle = -getAzimuth(north1Sensor, east1Sensor, 
        north2Sensor, interval, start, end);

    fireStateChange("Getting second east sensor orientation...");
    // direction north angle should be if north and east truly orthogonal
    // then east component is x component of rotation in that direction
    // i.e., need to correct by 90 degrees to get rotation angle rather than
    // azimuth of east sensor
    // offset by 3Pi/2 is the same as offset Pi/2 (90 degrees) in other 
    // rotation direction
    east2Angle = -getAzimuth(north1Sensor, east1Sensor, 
        east2Sensor, interval, start, end) + (3 * Math.PI / 2);
    
    // now to rotate the data according to these angles
    fireStateChange("Rotating data...");
    DataBlock north2Rotated =
        TimeSeriesUtils.rotate(ds.getBlock(4), ds.getBlock(5), north2Angle);
    stores[0].setBlock(1, north2Rotated);
    DataBlock east2Rotated = 
        TimeSeriesUtils.rotateX(ds.getBlock(4), ds.getBlock(5), east2Angle);
    stores[1].setBlock(1, east2Rotated);
    
    // now get the datasets to plug into the datastore
    String[] direction = new String[]{"north", "east", "vertical"};
    
    for (int i = 0; i < DIMS; ++i) {
      
      StringBuilder state = new StringBuilder("Running calculations on ");
      state.append(direction[i]);
      state.append(" components...");
      fireStateChange( state.toString() );
      
      componentBackends[i].runExperimentOnData(stores[i]);
      
    }
    
    for (Experiment exp : componentBackends) {
      // each backend only has one plot's worth of data
      // but is formatted as a list of per-plot data, so we use addAll
      xySeriesData.addAll( exp.getData() );
      // also get the names of the data going in for use w/ PDF, metadata
    }    

    for (int j = 0; j < componentBackends.length; ++j) {
      List<String> names = componentBackends[j].getInputNames();
      for (int i = 0; i < names.size(); i += 2) {
        dataNames.add( names.get(i) );
        dataNames.add( names.get(i + 1) );
      }
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
  
  private double getAzimuth(double[] n, double[] e, double[] r, 
      long interval, long start, long end) {
    // TODO: FIX THIS
    AzimuthExperiment azi = new AzimuthExperiment();
    azi.setSimple(true); // do the faster angle calculation
    azi.alternateEntryPoint(n, e, r, interval, start, end);
    return azi.getFitAngleRad();
  }
}
