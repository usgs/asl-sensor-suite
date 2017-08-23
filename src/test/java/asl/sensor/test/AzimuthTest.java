package asl.sensor.test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.io.FileNotFoundException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.TimeZone;

import org.junit.Test;

import asl.sensor.experiment.AzimuthExperiment;
import asl.sensor.gui.InputPanel;
import asl.sensor.input.DataStore;
import asl.sensor.utils.TimeSeriesUtils;

public class AzimuthTest {

  @Test
  public void getsCorrectAngle() {
    
    DataStore ds = new DataStore();
    
    String currentDir = System.getProperty("user.dir");
    String folder = currentDir + "/data/azi/";
    String[] prefixes = new String[3];
    prefixes[0] = "00_LH1";
    prefixes[1] = "00_LH2";
    prefixes[2] = "XX_LH1";
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
    
    AzimuthExperiment azi = new AzimuthExperiment();
    
    assertTrue( azi.hasEnoughData(ds) );
    
    SimpleDateFormat sdf = InputPanel.SDF;
    sdf.setTimeZone( TimeZone.getTimeZone("UTC") );
    // sdf.setLenient(false);
    
    Calendar cCal = Calendar.getInstance( sdf.getTimeZone() );
    cCal.setTimeInMillis( ds.getBlock(0).getStartTime() );
    cCal.set(Calendar.HOUR, 10);
    cCal.set(Calendar.MINUTE, 30);
    //System.out.println("start: " + sdf.format( cCal.getTime() ) );
    long start = cCal.getTime().getTime();
    cCal.set(Calendar.HOUR, 15);
    cCal.set(Calendar.MINUTE, 00);
    //System.out.println("end: " + sdf.format( cCal.getTime() ) );
    long end = cCal.getTime().getTime();
    
    ds.trim(start, end, 2);
    
    azi.runExperimentOnData(ds);
    
    System.out.println( azi.getFitAngle() );
    assertEquals( 16.0, azi.getFitAngle(), 0.5 );
    
  }
  
  @Test
  public void solvesNoRotation() {
    DataStore ds = new DataStore();
    
    String currentDir = System.getProperty("user.dir");
    String folder = currentDir + "/data/azimuth/";
    String[] prefixes = new String[3];
    prefixes[0] = "00_LH1";
    prefixes[1] = "00_LH2";
    prefixes[2] = "10_LH1";
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
    
    AzimuthExperiment azi = new AzimuthExperiment();
    
    assertTrue( azi.hasEnoughData(ds) );
    
    SimpleDateFormat sdf = InputPanel.SDF;
    sdf.setTimeZone( TimeZone.getTimeZone("UTC") );
    // sdf.setLenient(false);
    
    Calendar cCal = Calendar.getInstance( sdf.getTimeZone() );
    cCal.setTimeInMillis( ds.getBlock(0).getStartTime() );
    cCal.set(Calendar.HOUR_OF_DAY, 12);
    cCal.set(Calendar.MINUTE, 00);
    //System.out.println("start: " + sdf.format( cCal.getTime() ) );
    long start = cCal.getTime().getTime();
    cCal.set(Calendar.HOUR_OF_DAY, 14);
    cCal.set(Calendar.MINUTE, 00);
    //System.out.println("end: " + sdf.format( cCal.getTime() ) );
    long end = cCal.getTime().getTime();
    
    ds.trim(start, end, 2);
    
    azi.runExperimentOnData(ds);
    
    System.out.println( azi.getFitAngle() );
    assertEquals( 0., azi.getFitAngle(), 0.75 );
    
  }

}
  
