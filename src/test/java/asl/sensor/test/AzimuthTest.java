package asl.sensor.test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.io.FileNotFoundException;
import java.io.IOException;
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

  public String getNoisyData() {
    return "FUNA.00";
  }
  
  public String getCleanData() {
    return "ANMO.10";
  }
  
  @Test
  public void identifiesSprocketsAngle002Clean() {
    testsFromSprockets(2, getCleanData());
  }
  
  @Test
  public void identifiesSprocketsAngle002Noisy() {
    testsFromSprockets(2, getNoisyData());
  }
  
  @Test
  public void identifiesSprocketsAngle005Clean() {
    testsFromSprockets(5, getCleanData());
  }
  
  @Test
  public void identifiesSprocketsAngle005Noisy() {
    testsFromSprockets(5, getNoisyData());
  }
  
  @Test
  public void identifiesSprocketsAngle020Clean() {
    testsFromSprockets(20, getCleanData());
  }
  
  @Test
  public void identifiesSprocketsAngle020Noisy() {
    testsFromSprockets(20, getNoisyData());
  }
  
  @Test
  public void identifiesSprocketsAngle090Clean() {
    testsFromSprockets(90, getCleanData());
  }
  
  @Test
  public void identifiesSprocketsAngle090Noisy() {
    testsFromSprockets(90, getNoisyData());
  }
  
  @Test
  public void identifiesSprocketsAngle110Clean() {
    testsFromSprockets(110, getCleanData());
  }
  
  @Test
  public void identifiesSprocketsAngle110Noisy() {
    testsFromSprockets(110, getNoisyData());
  }
  
  @Test
  public void identifiesSprocketsAngle180Clean() {
    testsFromSprockets(180, getCleanData());
  }
  
  @Test
  public void identifiesSprocketsAngle180Noisy() {
    testsFromSprockets(180, getNoisyData());
  }
  
  @Test
  public void identifiesSprocketsAngle200Clean() {
    testsFromSprockets(200, getCleanData());
  }
  
  @Test
  public void identifiesSprocketsAngle200Noisy() {
    testsFromSprockets(200, getNoisyData());
  }
  
  @Test
  public void identifiesSprocketsAngle270Clean() {
    testsFromSprockets(270, getCleanData());
  }
  
  @Test
  public void identifiesSprocketsAngle270Noisy() {
    testsFromSprockets(270, getNoisyData());
  }
  
  @Test
  public void identifiesSprocketsAngle290Clean() {
    testsFromSprockets(290, getCleanData());
  }
  
  @Test
  public void identifiesSprocketsAngle290Noisy() {
    testsFromSprockets(290, getNoisyData());
  }
  
  @Test
  public void identifiesSprocketsAngle355Clean() {
    testsFromSprockets(355, getCleanData());
  }
  
  @Test
  public void identifiesSprocketsAngle355Noisy() {
    testsFromSprockets(355, getNoisyData());
  }
  
  @Test
  public void identifiesSprocketsAngle358Clean() {
    testsFromSprockets(358, getCleanData());
  }
  
  @Test
  public void identifiesSprocketsAngle358Noisy() {
    testsFromSprockets(358, getNoisyData());
  }
  
  public void testsFromSprockets(int angle, String staCha) {
    StringBuilder sb = new StringBuilder();
    sb.append(angle);
    // format for filenames is 002, 010, 358, etc.; prepend 0s if needed
    while(sb.length() < 3) {
      sb.insert(0, '0');
    }
    String data1 = "IU." + staCha + ".BH1";
    String data2 = "IU." + staCha + ".BH2";
    String refURL = "orientation/";
    // orientation/rotation/[002, etc.]/
    String testURL = refURL + "rotation/" + sb.toString() + "/";
    
    String refSubfolder = "orientation-reference/";
    String testSubfolder = "azimuth-" + sb.toString() + "/";
    
    try {
      // get reference data if needed
      TestUtils.downloadTestData(refURL, data1, refSubfolder, data1);
      // only need one reference point to get the data
      //TestUtils.downloadTestData(refURL, data2, refSubfolder, data2);
      // get test data if needed
      TestUtils.downloadTestData(testURL, data1, testSubfolder, data1);
      TestUtils.downloadTestData(testURL, data2, testSubfolder, data2);
    } catch (IOException e) {
      e.printStackTrace();
      fail();
    }
    
    DataStore ds = new DataStore();
    try {
      String root = "./test-data/sprockets/";
      testSubfolder = root + testSubfolder;
      refSubfolder = root + refSubfolder;
      ds.setBlock(0, TimeSeriesUtils.getFirstTimeSeries(testSubfolder + data1) );
      ds.setBlock(1, TimeSeriesUtils.getFirstTimeSeries(testSubfolder + data2) );
      ds.setBlock(2, TimeSeriesUtils.getFirstTimeSeries(refSubfolder + data1) );
      
      ds.trimToCommonTime();
      AzimuthExperiment ae = new AzimuthExperiment();
      ae.runExperimentOnData(ds);
      double fitAngle = ae.getFitAngle() - 180;
      if (fitAngle > 180) {
        fitAngle -= 360; // keep in range (-180, 180) for testing near 0 accurately
      } else if (angle > 180) {
        angle -= 360;
      }
      
      System.out.println(sb.toString() + " | " + ( (fitAngle % 360) + 360) % 360 );
      assertEquals((double) angle, fitAngle, 1.0);
      
    } catch (FileNotFoundException e) {
      e.printStackTrace();
      fail();
    }
    
  }
  
  @Test
  public void getsCorrectAngle() {
    
    DataStore ds = new DataStore();
    
    String currentDir = System.getProperty("user.dir");
    String folder = currentDir + "/test-data/azi-16off/";
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
    assertEquals( 16.0, azi.getFitAngle(), 2. );
    
  }
  
  @Test
  public void findsAntipolarCorrectly() {
    DataStore ds = new DataStore();
    
    String currentDir = System.getProperty("user.dir");
    String folder = currentDir + "/test-data/azi-polartest/";
    String[] prefixes = new String[3];
    prefixes[0] = "00_LH1";
    prefixes[1] = "00_LH2";
    prefixes[2] = "TST1.00_LH1";
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
    cCal.set(Calendar.HOUR_OF_DAY, 18);
    cCal.set(Calendar.MINUTE, 00);
    //System.out.println("start: " + sdf.format( cCal.getTime() ) );
    long start = cCal.getTime().getTime();
    cCal.set(Calendar.HOUR_OF_DAY, 20);
    cCal.set(Calendar.MINUTE, 30);
    //System.out.println("end: " + sdf.format( cCal.getTime() ) );
    long end = cCal.getTime().getTime();
    
    ds.trim(start, end, 2);
    
    System.out.println("FOR ANTIPOLAR TEST:");
    azi.runExperimentOnData(ds);
    
    System.out.println( azi.getFitAngle() );
    System.out.println("ANTIPOLAR TEST COMPLETED");
    assertEquals( 16., azi.getFitAngle(), 2. );
    
  }
  
  @Test
  public void solvesNoRotation() {
    DataStore ds = new DataStore();
    
    String currentDir = System.getProperty("user.dir");
    String folder = currentDir + "/test-data/azi-at0/";
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
    
    ds.trim(start, end);
    
    azi.runExperimentOnData(ds);
    double ang = azi.getFitAngle();
   
    if( ang > 180) {
      ang -= 360.;
    }
    System.out.println( ang );
    assertEquals( 0., ang, 2. );
    
  }

}
  
