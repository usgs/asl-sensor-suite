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


  // @Test // commented out until I can better understand the issues here
  public void identifiesSprocketsDay239AsOrtho() {
    DataStore ds;
    try {
      ds = getOrthogonalSprocketsData(239);
      Calendar stCal = TestUtils.getStartCalendar(ds);
      stCal.set(Calendar.HOUR_OF_DAY, 7);
      stCal.set(Calendar.MINUTE, 0);
      stCal.set(Calendar.SECOND, 0);
      stCal.set(Calendar.MILLISECOND, 0);
      long start = stCal.getTime().getTime();
      Calendar edCal = stCal;
      edCal.set(Calendar.HOUR_OF_DAY, 11);
      //Calendar edCal = TestUtils.getEndCalendar(ds);
      long end = edCal.getTime().getTime();
      ds.trim(start, end);
      identifiesOrthogonalFromSprockets(ds);
    } catch (IOException e) {
      e.printStackTrace();
      fail();
    }
  }

  @Test
  public void identifiesSprocketsDay243AsOrtho() {
    DataStore ds;
    try {
      ds = getOrthogonalSprocketsData(243);
      Calendar stCal = TestUtils.getStartCalendar(ds);
      long start = stCal.getTime().getTime();
      Calendar edCal = stCal;
      edCal.set(Calendar.HOUR_OF_DAY, 14);
      edCal.set(Calendar.MINUTE, 0);
      edCal.set(Calendar.SECOND, 0);
      edCal.set(Calendar.MILLISECOND, 0);
      long end = edCal.getTime().getTime();
      ds.trim(start, end);
      identifiesOrthogonalFromSprockets(ds);
    } catch (IOException e) {
      e.printStackTrace();
      fail();
    }
  }

  public DataStore getOrthogonalSprocketsData(int day) throws IOException {
    // these are true for both test and ref data;
    String data1 = "00_LH1.512.seed";
    String data2 = "00_LH2.512.seed";
    String test1 = "00_LH1.512.seed";
    String test2 = "00_LH2.512.seed";
    String testPrepend = "TEST-";
    String refPrepend = "REF-";

    String testURL = "orientation/sensor/orthogonal/"+day+"/";
    String refURL = "orientation/sensor/reference/"+day+"/";
    String testSubfolder = "is-orthogonal/"+day+"/";

    // download data (if it doesn't already exist locally)
    // get reference data 1 and 2
    TestUtils.downloadTestData(refURL, data1, testSubfolder, refPrepend+data1);
    TestUtils.downloadTestData(refURL, data2, testSubfolder, refPrepend+data2);
    // get test data 1 and 2
    TestUtils.downloadTestData(testURL, test1, testSubfolder, testPrepend+data1);
    TestUtils.downloadTestData(testURL, test2, testSubfolder, testPrepend+data2);


    DataStore ds = new DataStore();
    // TODO: check if working under windows too?
    String root = "./test-data/sprockets/" + testSubfolder;

    String filePre = root + refPrepend;
    ds.setBlock(0, TimeSeriesUtils.getFirstTimeSeries(filePre+data1) );
    ds.setBlock(1, TimeSeriesUtils.getFirstTimeSeries(filePre+data2) );
    filePre = root + testPrepend;
    ds.setBlock(2, TimeSeriesUtils.getFirstTimeSeries(filePre+data1) );
    ds.setBlock(3, TimeSeriesUtils.getFirstTimeSeries(filePre+data2) );
    
    System.out.println(ds.getBlock(2).getName());
    System.out.println(ds.getBlock(3).getName());
    /*
      Calendar cCal = TestUtils.getStartCalendar(ds);
      cCal.set(Calendar.HOUR_OF_DAY, 13);
      cCal.set(Calendar.MINUTE, 0);
      cCal.set(Calendar.SECOND, 0);
      cCal.set(Calendar.MILLISECOND, 0);
      long start = cCal.getTime().getTime();
      cCal.set(Calendar.HOUR_OF_DAY, 18);
      long end = TestUtils.getEndCalendar(ds).getTime().getTime();
     */
    return ds;

  }

  public void identifiesOrthogonalFromSprockets(DataStore ds) {
    double fitAngle = 0;
    OrthogonalExperiment oe = new OrthogonalExperiment();
    oe.runExperimentOnData(ds);
    fitAngle = oe.getFitAngle();
    System.out.println(Arrays.toString(oe.getSolutionParams()));

    System.out.println("SHOULD BE 90! - " + fitAngle);
    assertEquals(90., fitAngle, 2.0);
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
      System.out.println("SHOULD NOT BE 90! - " + fitAngle);
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
    assertEquals( 94., orth.getFitAngle(), 1. );

  }

}
