package asl.sensor.test;

import static org.junit.Assert.*;
import java.io.FileNotFoundException;
import java.io.IOException;
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
  public void identifiesOrthogonalFromSprockets() {
    int[] days = new int[]{239, 243};
    double[] fitAngles = new double[days.length];
    // these are true for both test and ref data;
    String data1 = "00_LH1.512.seed";
    String data2 = "00_LH2.512.seed";
    String test1 = "00_LH1.512.seed";
    String test2 = "00_LH2.512.seed";
    String testPrepend = "TEST-";
    String refPrepend = "REF-";
    for (int i = 0; i < days.length; ++i) {
      int day = days[i];
      String testURL = "orientation/sensor/orthogonal/"+day+"/";
      String refURL = "orientation/sensor/reference/"+day+"/";
      String testSubfolder = "is-orthogonal/"+day+"/";

      try {
        // download data (if it doesn't already exist locally)
        // get reference data 1 and 2
        TestUtils.downloadTestData(refURL, data1, testSubfolder, refPrepend+data1);
        TestUtils.downloadTestData(refURL, data2, testSubfolder, refPrepend+data2);
        // get test data 1 and 2
        TestUtils.downloadTestData(testURL, test1, testSubfolder, testPrepend+data1);
        TestUtils.downloadTestData(testURL, test2, testSubfolder, testPrepend+data2);
      } catch (IOException e) {
        e.printStackTrace();
        fail();
      }
      
      DataStore ds = new DataStore();
      // TODO: check if working under windows too?
      String root = "./test-data/sprockets/" + testSubfolder;
      try {
        String filePre = root + refPrepend;
        ds.setBlock(0, TimeSeriesUtils.getFirstTimeSeries(filePre+data1) );
        ds.setBlock(1, TimeSeriesUtils.getFirstTimeSeries(filePre+data2) );
        filePre = root + testPrepend;
        ds.setBlock(2, TimeSeriesUtils.getFirstTimeSeries(filePre+data1) );
        ds.setBlock(3, TimeSeriesUtils.getFirstTimeSeries(filePre+data2) );
        Calendar cCal = TestUtils.getStartCalendar(ds);
        cCal.set(Calendar.HOUR_OF_DAY, 14);
        cCal.set(Calendar.MINUTE, 0);
        cCal.set(Calendar.SECOND, 0);
        cCal.set(Calendar.MILLISECOND, 0);
        long start = cCal.getTime().getTime();
        cCal.set(Calendar.HOUR_OF_DAY, 18);
        long end = cCal.getTime().getTime();
        ds.trim(start, end);
        OrthogonalExperiment oe = new OrthogonalExperiment();
        oe.runExperimentOnData(ds);
        fitAngles[i] = oe.getFitAngle();
        System.out.println(Arrays.toString(oe.getSolutionParams()));
      } catch (FileNotFoundException e) {
        e.printStackTrace();
        fail();
      }

    }
    boolean allEqual = true;
    for (double fitAngle : fitAngles) {
      System.out.println("SHOULD NOT BE 90!:" + fitAngle);
      allEqual &= (90. - 2. < fitAngle) && (90. + 2. > fitAngle);
    }
    assertTrue(allEqual);
  }
  
  @Test
  public void identifiesNotOrthogonalFromSprockets() {
    int[] days = new int[]{246, 249, 251, 259};
    double[] fitAngles = new double[days.length];
    // these are true for both test and ref data;
    String data1 = "00_LH1.512.seed";
    String data2 = "00_LH2.512.seed";
    String test1 = "00_LH1.512.seed";
    String test2 = "00_LH2.512.seed";
    String testPrepend = "TEST-";
    String refPrepend = "REF-";
    for (int i = 0; i < days.length; ++i) {
      int day = days[i];
      String testURL = "orientation/sensor/non-orthogonal/"+day+"/";
      String refURL = "orientation/sensor/reference/"+day+"/";
      String testSubfolder = "is-orthogonal/"+day+"/";

      try {
        // download data (if it doesn't already exist locally)
        // get reference data 1 and 2
        TestUtils.downloadTestData(refURL, data1, testSubfolder, refPrepend+data1);
        TestUtils.downloadTestData(refURL, data2, testSubfolder, refPrepend+data2);
        // get test data 1 and 2
        TestUtils.downloadTestData(testURL, test1, testSubfolder, testPrepend+data1);
        TestUtils.downloadTestData(testURL, test2, testSubfolder, testPrepend+data2);
      } catch (IOException e) {
        e.printStackTrace();
        fail();
      }
      
      DataStore ds = new DataStore();
      // TODO: check if working under windows too?
      String root = "./test-data/sprockets/" + testSubfolder;
      try {
        String filePre = root + refPrepend;
        ds.setBlock(0, TimeSeriesUtils.getFirstTimeSeries(filePre+data1) );
        ds.setBlock(1, TimeSeriesUtils.getFirstTimeSeries(filePre+data2) );
        filePre = root + testPrepend;
        ds.setBlock(2, TimeSeriesUtils.getFirstTimeSeries(filePre+data1) );
        ds.setBlock(3, TimeSeriesUtils.getFirstTimeSeries(filePre+data2) );
        Calendar cCal = TestUtils.getStartCalendar(ds);
        cCal.set(Calendar.HOUR_OF_DAY, 14);
        cCal.set(Calendar.MINUTE, 0);
        cCal.set(Calendar.SECOND, 0);
        cCal.set(Calendar.MILLISECOND, 0);
        long start = cCal.getTime().getTime();
        cCal.set(Calendar.HOUR_OF_DAY, 18);
        long end = cCal.getTime().getTime();
        ds.trim(start, end);
        OrthogonalExperiment oe = new OrthogonalExperiment();
        oe.runExperimentOnData(ds);
        fitAngles[i] = oe.getFitAngle();
        System.out.println(Arrays.toString(oe.getSolutionParams()));
      } catch (FileNotFoundException e) {
        e.printStackTrace();
        fail();
      }

    }
    
    // test should fail if any of the orthogonal within 2 degrees of 90
    boolean anyOrthogonal = false;
    for (double fitAngle : fitAngles) {
      System.out.println(fitAngle);
      anyOrthogonal |= (90. - 2. < fitAngle) && (90. + 2. > fitAngle);
    }
    assertFalse(anyOrthogonal);
  }
  

  @Test
  public void getsCorrectAngle() {
    
    DataStore ds = new DataStore();
    
    String currentDir = System.getProperty("user.dir");
    String folder = currentDir + "/test-data/orthog-94/";
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
    cCal.setTimeInMillis( ds.getBlock(0).getStartTime() );
    cCal.set(Calendar.HOUR, 7);
    // cCal.set(Calendar.MINUTE, 30);
    System.out.println("start: " + sdf.format( cCal.getTime() ) );
    long start = cCal.getTime().getTime();
    cCal.set(Calendar.HOUR, 13);
    cCal.set(Calendar.MINUTE, 00);
    System.out.println("end: " + sdf.format( cCal.getTime() ) );
    long end = cCal.getTime().getTime();
    
    ds.trim(start, end);
    
    orth.runExperimentOnData(ds);
    
    System.out.println( orth.getFitAngle() );
    System.out.println( Arrays.toString( orth.getSolutionParams() ) );
    assertEquals( 95., orth.getFitAngle(), 1. );
    
  }
  
}
