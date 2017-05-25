package asl.sensor.experiment;

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

  
  public NoiseNineExperiment() {
    super();
    // indices are fixed since we need all 9 data points here
    respIndices = new int[]{1, 2, 3, 4, 5, 6, 7, 8, 9};
  }
  
  @Override
  protected void backend(DataStore ds) {
    
    // NOTE: this may need to change in the event of a test using > 9 inputs
    for (int i = 0; i < 9; ++i) {
      // doing this loop here saves us time and significant lines of code
      dataNames.add( ds.getBlock(0).getName() );
      dataNames.add( ds.getResponse(0).getName() );
    }
    
    // get the components
    // why unroll the datastore's contents like this?
    // since we need to do the azimuth calculations with subsets of the data
    // and then rotate some of those sensors to get new datablocks
    // we can only compact the code so hard, and this is an easier arrangement
    // than, say, trying to index into a series of arraylists
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
    
    DataBlock north3Sensor = ds.getBlock(6);
    InstrumentResponse north3Resp = ds.getResponse(6);
    DataBlock east3Sensor = ds.getBlock(7);
    InstrumentResponse east3Resp = ds.getResponse(7);
    DataBlock vert3Sensor = ds.getBlock(8);
    InstrumentResponse vert3Resp = ds.getResponse(8);
    
    
    StringBuilder sb = new StringBuilder();
    sb.append("Beginning rotations (offset angle esimates)\n");
    sb.append("for second and third sets of horizontal (N, E) data...");
    fireStateChange( sb.toString() );
    
    // set angles and then rotate data
    AzimuthExperiment azi = new AzimuthExperiment();
    azi.setSimple(true); // use rough estimate of coherence, no windowing
    DataStore aziStore = new DataStore();
    aziStore.setData(0, north1Sensor);
    aziStore.setData(1, east1Sensor);
    
    // angle should be set negative -- rotate third sensor, not the opposite
    aziStore.setData(2, north2Sensor);
    azi.runExperimentOnData(aziStore);
    double north2Angle = -azi.getFitAngleRad();
    fireStateChange("Found orientation of second north sensor!");
    
    aziStore.setData(2, east2Sensor);
    azi.runExperimentOnData(aziStore);
    double east2Angle = -azi.getFitAngleRad() - (Math.PI / 2);
    fireStateChange("Found orientation of second east sensor!");
    // need to offset rotation by 90 degrees -- don't want it facing north
    
    // same as above
    aziStore.setData(2, north3Sensor);
    azi.runExperimentOnData(aziStore);
    double north3Angle = -azi.getFitAngleRad();
    fireStateChange("Found orientation of third north sensor!");
    
    aziStore.setData(2, east2Sensor);
    azi.runExperimentOnData(aziStore);
    double east3Angle = -azi.getFitAngleRad() - (Math.PI / 2);
    fireStateChange("Found orientation of third east sensor!");
    
    // now to rotate the data according to these angles
    DataBlock north2Rotated =
        TimeSeriesUtils.rotate(north2Sensor, east2Sensor, north2Angle);
    DataBlock east2Rotated = 
        TimeSeriesUtils.rotate(east2Sensor, north2Sensor, east2Angle);
    DataBlock north3Rotated =
        TimeSeriesUtils.rotate(north3Sensor, east3Sensor, north3Angle);
    DataBlock east3Rotated =
        TimeSeriesUtils.rotate(east3Sensor, north3Sensor, east3Angle);
    fireStateChange("All offset horizontal data rotated!");
    
    // set components into N,E,Z directional subcomponents
    
    DataStore northComponents = new DataStore();
    northComponents.setData(0, north1Sensor);
    northComponents.setResponse(0, north1Resp);
    northComponents.setData(1, north2Rotated);
    northComponents.setResponse(1, north2Resp);
    northComponents.setData(2, north3Rotated);
    northComponents.setResponse(2, north3Resp);
    
    DataStore eastComponents = new DataStore();
    eastComponents.setData(0, east1Sensor);
    eastComponents.setResponse(0, east1Resp);
    eastComponents.setData(1, east2Rotated);
    eastComponents.setResponse(1, east2Resp);
    eastComponents.setData(2, east3Rotated);
    eastComponents.setResponse(2, east3Resp);
    
    DataStore vertComponents = new DataStore();
    vertComponents.setData(0, vert1Sensor);
    vertComponents.setResponse(0, vert1Resp);
    vertComponents.setData(1, vert2Sensor);
    vertComponents.setResponse(1, vert2Resp);
    vertComponents.setData(2, vert3Sensor);
    vertComponents.setResponse(2, vert3Resp);
    
    // get noise from each axis's data
    NoiseExperiment noiseExp = new NoiseExperiment();
    fireStateChange("Calculating north noise components...");
    noiseExp.setFreqSpace(freqSpace);
    noiseExp.runExperimentOnData(northComponents);
    XYSeriesCollection northXYS = noiseExp.getData().get(0);
    fireStateChange("Calculating east noise components...");
    noiseExp.runExperimentOnData(eastComponents);
    XYSeriesCollection eastXYS = noiseExp.getData().get(0);
    fireStateChange("Calculating vertical noise components...");
    noiseExp.runExperimentOnData(vertComponents);
    XYSeriesCollection vertXYS = noiseExp.getData().get(0);

    xySeriesData.add(northXYS);
    xySeriesData.add(eastXYS);
    xySeriesData.add(vertXYS);
    
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


}
