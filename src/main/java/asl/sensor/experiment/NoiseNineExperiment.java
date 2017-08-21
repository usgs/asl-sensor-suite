package asl.sensor.experiment;

import java.util.List;

import org.jfree.data.xy.XYSeriesCollection;

import asl.sensor.input.DataBlock;
import asl.sensor.input.DataStore;
import asl.sensor.input.InstrumentResponse;
import asl.sensor.utils.TimeSeriesUtils;

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
 * @author akearns
 *
 */
public class NoiseNineExperiment extends NoiseExperiment {

  private static final int DIMS = 3;
  private double[] northAngles, eastAngles;
  
  public NoiseNineExperiment() {
    super();
    // indices are fixed since we need all 9 data points here
    respIndices = new int[]{0, 1, 2, 3, 4, 5, 6, 7, 8};
  }
  
  @Override
  protected void backend(DataStore ds) {
    
    northAngles = new double[2];
    eastAngles = new double[2];
    
    // NOTE: this may need to change in the event of a test using > 9 inputs
    for (int i = 0; i < 9; ++i) {
      // doing this loop here saves us time and significant lines of code
      dataNames.add( ds.getBlock(i).getName() );
      dataNames.add( ds.getResponse(i).getName() );
    }
    
    DataStore[] stores = new DataStore[DIMS];
    
    for (int i = 0; i < DIMS; ++i) {
      stores[i] = new DataStore();
      for (int j = 0; j < 3; ++j) {
        stores[i].setBlock( j, ds.getBlock( i + (j * DIMS) ) );
        stores[i].setResponse( j, ds.getResponse( i + (j * DIMS) ) );
      }
    }
    
    // get the components
    // why unroll the datastore's contents like this?
    // since we need to do the azimuth calculations with subsets of the data
    // and then rotate some of those sensors to get new datablocks
    // we can only compact the code so hard, and this is an easier arrangement
    // than, say, trying to index into a series of arraylists
    List<Number> north1Sensor = ds.getBlock(0).getData();
    List<Number> east1Sensor = ds.getBlock(1).getData();
    
    List<Number> north2Sensor = ds.getBlock(3).getData();
    List<Number> east2Sensor = ds.getBlock(4).getData();
    
    List<Number> north3Sensor = ds.getBlock(6).getData();
    List<Number> east3Sensor = ds.getBlock(7).getData();
    
    long interval = ds.getBlock(0).getInterval();
    long start = ds.getBlock(0).getStartTime();
    long end = ds.getBlock(0).getEndTime();
    
    StringBuilder sb = new StringBuilder();
    sb.append("Beginning rotations (offset angle esimates)\n");
    sb.append("for second and third sets of horizontal (N, E) data...");
    fireStateChange( sb.toString() );
    
    // see also the rotation used in the 9-input self noise backend
    fireStateChange("Getting second north sensor orientation...");
    northAngles[0] = -getAzimuth(north1Sensor, east1Sensor, 
        north2Sensor, interval, start, end);

    fireStateChange("Getting second east sensor orientation...");
    // direction north angle should be if north and east truly orthogonal
    // then east component is x component of rotation in that direction
    // i.e., need to correct by 90 degrees to get rotation angle rather than
    // azimuth of east sensor
    // offset by 3Pi/2 is the same as offset Pi/2 (90 degrees) in other 
    // rotation direction
    eastAngles[0] = -getAzimuth(north1Sensor, east1Sensor, 
        east2Sensor, interval, start, end) + (3 * Math.PI / 2);
    
    // now to rotate the data according to these angles
    fireStateChange("Rotating data...");
    DataBlock north2Rotated =
        TimeSeriesUtils.rotate(ds.getBlock(4), ds.getBlock(5), northAngles[0]);
    stores[0].setBlock(1, north2Rotated);
    DataBlock east2Rotated = 
        TimeSeriesUtils.rotateX(ds.getBlock(4), ds.getBlock(5), eastAngles[0]);
    stores[1].setBlock(1, east2Rotated);
    
    // see also the rotation used in the 9-input self noise backend
    fireStateChange("Getting third north sensor orientation...");
    northAngles[1] = -getAzimuth(north1Sensor, east1Sensor, 
        north3Sensor, interval, start, end);

    fireStateChange("Getting third east sensor orientation...");
    // direction north angle should be if north and east truly orthogonal
    // then east component is x component of rotation in that direction
    // i.e., need to correct by 90 degrees to get rotation angle rather than
    // azimuth of east sensor
    // offset by 3Pi/2 is the same as offset Pi/2 (90 degrees) in other 
    // rotation direction
    eastAngles[1] = -getAzimuth(north1Sensor, east1Sensor, 
        east3Sensor, interval, start, end) + (3 * Math.PI / 2);
    
    
    // now to rotate the data according to these angles
    fireStateChange("Rotating data...");
    DataBlock north3Rotated =
        TimeSeriesUtils.rotate(ds.getBlock(7), ds.getBlock(8), northAngles[1]);
    stores[0].setBlock(2, north3Rotated);
    DataBlock east3Rotated = 
        TimeSeriesUtils.rotateX(ds.getBlock(7), ds.getBlock(8), eastAngles[1]);
    stores[1].setBlock(2, east3Rotated);
    fireStateChange("All offset horizontal data rotated!");
    
    // set components into N,E,Z directional subcomponents
    
    // get noise from each axis's data
    NoiseExperiment noiseExp = new NoiseExperiment();
    noiseExp.setFreqSpace(freqSpace);
    String[] directions = new String[]{"north","east","vertical"};
    for (int i = 0; i < DIMS; ++i) {
      sb = new StringBuilder("Calculating ");
      sb.append(directions[i]);
      sb.append(" noise components...");
      noiseExp.runExperimentOnData(stores[i]);
      XYSeriesCollection xys = noiseExp.getData().get(0);
      xySeriesData.add(xys);
    }
    
  }

  /**
   * Return array of angles (degree-valued) which north components have been
   * rotated by, starting with the second component (1st north component is
   * assumed to have zero rotation)
   * @return double array representing angles in degrees
   */
  public double[] getNorthAngles() {
    double[] returningNorth = new double[northAngles.length];
    for (int i = 0; i < returningNorth.length; ++i) {
      returningNorth[i] = Math.toDegrees(northAngles[i]) % 360;
    }
    return returningNorth;
  }
  
  /**
   * Return array of angles (degree-valued) which east components have been
   * rotated by, starting with the second component (1st east component is
   * assumed to have zero rotation)
   * @return double array representing angles in degrees
   */
  public double[] getEastAngles() {
    double[] returningEast = new double[eastAngles.length];
    for (int i = 0; i < returningEast.length; ++i) {
      returningEast[i] = Math.toDegrees(eastAngles[i]) % 360;
    }
    return returningEast;
  }
  
  @Override
  public int blocksNeeded() {
    return 9;
  }

  @Override
  public boolean hasEnoughData(DataStore ds) {
    for (int i = 0; i < blocksNeeded(); ++i) {
      if ( !ds.bothComponentsSet(i) ) {
        return false;
      }
    }
    return true;
  }

  
  private double getAzimuth(List<Number> n, List<Number> e, List<Number> r, 
      long interval, long start, long end) {
    // TODO: FIX THIS
    AzimuthExperiment azi = new AzimuthExperiment();
    azi.setSimple(true); // do the faster angle calculation
    azi.alternateEntryPoint(n, e, r, interval, start, end);
    return azi.getFitAngleRad();
  }
}
