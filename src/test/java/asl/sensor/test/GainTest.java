package asl.sensor.test;

import static org.junit.Assert.*;

import java.io.FileNotFoundException;
import java.util.ArrayList;

import org.junit.Test;

import asl.sensor.experiment.GainExperiment;
import asl.sensor.input.DataStore;
import asl.sensor.utils.TimeSeriesUtils;

public class GainTest {

  @Test
  public void testGainCalculation() {
    
    DataStore ds = new DataStore();
    
    String currentDir = System.getProperty("user.dir");
    String folder = currentDir + "/data/relativeGain/";
    String[] prefixes = new String[2];
    prefixes[0] = "00_BHZ";
    prefixes[1] = "10_BHZ";
    String extension = ".512.seed";
    
    for (int i = 0; i < prefixes.length; ++i) {
      String fName = folder + prefixes[i] + extension;
      String seriesName = "";
      try {
        seriesName = 
            new ArrayList<String>( TimeSeriesUtils.getMplexNameSet(fName) ).
            get(0);
      } catch (FileNotFoundException e) {
        // TODO Auto-generated catch block
        fail();
        e.printStackTrace();
      }
      ds.setData(i, fName, seriesName);
    }
    
    folder = currentDir + "/responses/relativeGain/";
    String[] rnames = new String[2];
    rnames[0] = "RESP.IU.ANMO.00.BHZ_gainx100";
    rnames[1] = "RESP.IU.ANMO.10.BHZ";
    
    for (int i = 0; i < rnames.length; ++i) {
      String fName = folder + rnames[i];
      ds.setResponse(i, fName);
    }
    
    GainExperiment ge = new GainExperiment();
    ge.runExperimentOnData(ds);
    
    double[] stats = ge.getStatsFromPeak(0);
    double gain = stats[3];
    
    assertEquals(gain, 11714., 1.0);
    
  }
}
