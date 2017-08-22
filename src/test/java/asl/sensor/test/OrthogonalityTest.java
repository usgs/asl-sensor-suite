package asl.sensor.test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.io.FileNotFoundException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Calendar;
import java.util.TimeZone;

import org.junit.Test;

import asl.sensor.experiment.OrthogonalExperiment;
import asl.sensor.gui.InputPanel;
import asl.sensor.input.DataStore;
import asl.sensor.utils.TimeSeriesUtils;

public class OrthogonalityTest {

  @Test
  public void getsCorrectAngle() {
    
    DataStore ds = new DataStore();
    
    String currentDir = System.getProperty("user.dir");
    String folder = currentDir + "/data/orthogonality/";
    String[] prefixes = new String[4];
    prefixes[0] = "00_LH1";
    prefixes[1] = "00_LH2";
    prefixes[2] = "10_LH1";
    prefixes[3] = "10_LH2";
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
      ds.setBlock(i, fName, seriesName);
    }
    
    OrthogonalExperiment orth = new OrthogonalExperiment();
    
    assertTrue( orth.hasEnoughData(ds) );
    
    SimpleDateFormat sdf = InputPanel.SDF;
    sdf.setTimeZone( TimeZone.getTimeZone("UTC") );
    // sdf.setLenient(false);
    
    Calendar cCal = Calendar.getInstance( sdf.getTimeZone() );
    cCal.setTimeInMillis( ds.getBlock(0).getStartTime() / 1000 );
    cCal.set(Calendar.HOUR, 7);
    // cCal.set(Calendar.MINUTE, 30);
    System.out.println("start: " + sdf.format( cCal.getTime() ) );
    long start = cCal.getTime().getTime() * 1000L;
    cCal.set(Calendar.HOUR, 13);
    cCal.set(Calendar.MINUTE, 00);
    System.out.println("end: " + sdf.format( cCal.getTime() ) );
    long end = cCal.getTime().getTime() * 1000L;
    
    ds.trim(start, end);
    
    orth.runExperimentOnData(ds);
    
    System.out.println( orth.getFitAngle() );
    System.out.println( Arrays.toString( orth.getSolutionParams() ) );
    assertEquals( 94., orth.getFitAngle(), 1. );
    
  }
  
}
