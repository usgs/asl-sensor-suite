package asl.sensor.experiment;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import asl.sensor.gui.ExperimentPanel;
import asl.sensor.input.DataStore;
import asl.sensor.test.TestUtils;
import asl.sensor.utils.TimeSeriesUtils;
import edu.iris.dmc.seedcodec.CodecException;
import edu.sc.seis.seisFile.mseed.SeedFormatException;
import java.io.IOException;
import java.util.Arrays;
import java.util.Calendar;
import java.util.Random;
import org.junit.Test;

public class AzimuthExperimentTest {

  public static String folder = TestUtils.TEST_DATA_LOCATION + TestUtils.SUBPAGE;

  @Test
  public void matchArrayLengths_offByOne() {
    double[] arr1 = new double[41753];
    double[] arr2 = new double[41754];
    double[] arr3 = new double[41754];
    double[][] results = AzimuthExperiment.matchArrayLengths(arr1, arr2, arr3);

    assertEquals(3, results.length);
    double[] testArr1 = results[0];
    double[] testArr2 = results[1];
    double[] testArr3 = results[2];

    assertEquals(testArr1.length, 41753);
    assertEquals(testArr2.length, 41753);
    assertEquals(testArr3.length, 41753);
  }

  @Test
  public void matchArrayLengths_multipleLow() {
    double[] arr1 = new double[41855];
    double[] arr2 = new double[41754];
    double[] arr3 = new double[41751];
    double[][] results = AzimuthExperiment.matchArrayLengths(arr1, arr2, arr3);

    assertEquals(3, results.length);
    double[] testArr1 = results[0];
    double[] testArr2 = results[1];
    double[] testArr3 = results[2];

    assertEquals(testArr1.length, 41751);
    assertEquals(testArr2.length, 41751);
    assertEquals(testArr3.length, 41751);
  }

  @Test
  public void matchArrayLengths_offByMany() {
    double[] arr1 = new double[3000];
    double[] arr2 = new double[41753];
    double[] arr3 = new double[41754];
    double[][] results = AzimuthExperiment.matchArrayLengths(arr1, arr2, arr3);

    assertEquals(3, results.length);
    double[] testArr1 = results[0];
    double[] testArr2 = results[1];
    double[] testArr3 = results[2];

    assertEquals(testArr1.length, 3000);
    assertEquals(testArr2.length, 3000);
    assertEquals(testArr3.length, 3000);
  }

  @Test
  public void matchArrayLengths_equalLengths() {
    double[] arr1 = new double[41753];
    double[] arr2 = new double[41753];
    double[] arr3 = new double[41753];
    double[][] results = AzimuthExperiment.matchArrayLengths(arr1, arr2, arr3);

    assertEquals(3, results.length);
    double[] testArr1 = results[0];
    double[] testArr2 = results[1];
    double[] testArr3 = results[2];

    assertEquals(testArr1.length, 41753);
    assertEquals(testArr2.length, 41753);
    assertEquals(testArr3.length, 41753);
  }

  @Test
  public void matchArrayLengths_singleArray() {
    double[] arr1 = new double[41753];
    double[][] results = AzimuthExperiment.matchArrayLengths(arr1);

    assertEquals(1, results.length);
    double[] testArr1 = results[0];

    assertEquals(testArr1.length, 41753);
  }

  @Test
  public void findsAntipolarCorrectly() {
    DataStore ds = new DataStore();

    String dataFolder = folder + "azi-polartest/";
    String[] prefixes = new String[3];
    prefixes[0] = "00_LH1";
    prefixes[1] = "00_LH2";
    prefixes[2] = "TST1.00_LH1";
    String extension = ".512.seed";

    for (int i = 0; i < prefixes.length; ++i) {
      String fName = dataFolder + prefixes[i] + extension;
      try {
        ds.setBlock(i, fName);
      } catch (SeedFormatException | CodecException | IOException e) {
        e.printStackTrace();
        fail();
      }
    }

    AzimuthExperiment azi = new AzimuthExperiment();

    assertTrue(azi.hasEnoughData(ds));

    Calendar cCal = Calendar.getInstance(ExperimentPanel.DATE_TIME_FORMAT.get().getTimeZone());
    cCal.setTimeInMillis(ds.getBlock(0).getStartTime());
    cCal.set(Calendar.HOUR_OF_DAY, 18);
    cCal.set(Calendar.MINUTE, 0);
    long start = cCal.getTime().getTime();
    cCal.set(Calendar.HOUR_OF_DAY, 20);
    cCal.set(Calendar.MINUTE, 30);
    long end = cCal.getTime().getTime();

    ds.trim(start, end, 2);

    azi.runExperimentOnData(ds);
    assertEquals(16., azi.getFitAngle(), 2.);

  }

  private String getCleanData() {
    return "ANMO.10";
  }

  private String getNoisyData() {
    return "FUNA.00";
  }

  @Test
  public void getsCorrectAngleANMO() {

    DataStore ds = new DataStore();

    String dataFolder = folder + "azi-ANMO-test/";
    String[] prefixes = new String[3];
    prefixes[0] = "ANMO.00_LH1";
    prefixes[1] = "ANMO.00_LH2";
    prefixes[2] = "TST.00_LH1";
    String extension = ".512.seed";

    for (int i = 0; i < prefixes.length; ++i) {
      String fName = dataFolder + prefixes[i] + extension;
      try {
        ds.setBlock(i, fName);
      } catch (SeedFormatException | CodecException | IOException e) {
        e.printStackTrace();
        fail();
      }
    }

    AzimuthExperiment azi = new AzimuthExperiment();

    assertTrue(azi.hasEnoughData(ds));

    azi.runExperimentOnData(ds);
    assertEquals(15.0, azi.getFitAngle(), 1.);
    assertEquals(0.4, azi.getUncertainty(), 0.5);

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

  @Test
  public void solvesNoRotation() {
    DataStore ds = new DataStore();

    String dataFolder = folder + "azi-at0/";
    String[] prefixes = new String[3];
    prefixes[0] = "00_LH1";
    prefixes[1] = "00_LH2";
    prefixes[2] = "10_LH1";
    String extension = ".512.seed";

    for (int i = 0; i < prefixes.length; ++i) {
      String fName = dataFolder + prefixes[i] + extension;
      try {
        ds.setBlock(i, fName);
      } catch (SeedFormatException | CodecException | IOException e) {
        e.printStackTrace();
        fail();
      }
    }

    AzimuthExperiment azi = new AzimuthExperiment();

    assertTrue(azi.hasEnoughData(ds));

    Calendar cCal = Calendar.getInstance(ExperimentPanel.DATE_TIME_FORMAT.get().getTimeZone());
    cCal.setTimeInMillis(ds.getBlock(0).getStartTime());
    cCal.set(Calendar.HOUR_OF_DAY, 12);
    cCal.set(Calendar.MINUTE, 0);
    long start = cCal.getTime().getTime();
    cCal.set(Calendar.HOUR_OF_DAY, 14);
    cCal.set(Calendar.MINUTE, 0);
    long end = cCal.getTime().getTime();

    ds.trim(start, end);

    azi.runExperimentOnData(ds);
    double ang = azi.getFitAngle();

    if (ang > 180) {
      ang -= 360.;
    }
    assertEquals(0., ang, 2.);

  }

  @Test
  public void testGetsCorrectWindowAnglesCOWI() {
    String filename = folder + "cowi-multitests/C100823215422_COWI.LHx";
    String dataname1 = "US_COWI_  _LHN";
    String dataname2 = "US_COWI_  _LHE";
    String filename2 = folder + "cowi-multitests/DT000110.LH1";
    try {
      DataStore ds = new DataStore();
      ds.setBlock(0, filename, dataname1);
      ds.setBlock(1, filename, dataname2);
      ds.setBlock(2, filename2);

      String startString = "2010-236T02:00:00.0";
      String endString = "2010-236T13:00:00.0";
      long st = TestUtils.timeStringToEpochMilli(startString);
      long ed = TestUtils.timeStringToEpochMilli(endString);
      ds.trim(st, ed);

      AzimuthExperiment az = new AzimuthExperiment();
      az.runExperimentOnData(ds);

      assertEquals(3.2, az.getFitAngle(), 0.5);
      assertEquals(0.5, az.getUncertainty(), 0.5);

    } catch (SeedFormatException | CodecException | IOException e) {
      e.printStackTrace();
      fail();
    }
  }

  private void testsFromSprockets(int angle, String staCha) {
    StringBuilder sb = new StringBuilder();
    sb.append(angle);
    // format for filenames is 002, 010, 358, etc.; prepend 0s if needed
    while (sb.length() < 3) {
      sb.insert(0, '0');
    }
    String data1 = "IU." + staCha + ".BH1";
    String data2 = "IU." + staCha + ".BH2";

    String refSubfolder = "mock_data/orientation/";
    String testSubfolder = refSubfolder + "rotation/" + sb.toString() + "/";

    DataStore ds = new DataStore();
    String root = TestUtils.TEST_DATA_LOCATION;
    testSubfolder = root + testSubfolder;
    refSubfolder = root + refSubfolder;
    try {
      ds.setBlock(0, TimeSeriesUtils.getFirstTimeSeries(testSubfolder + data1));
      ds.setBlock(1, TimeSeriesUtils.getFirstTimeSeries(testSubfolder + data2));
      ds.setBlock(2, TimeSeriesUtils.getFirstTimeSeries(refSubfolder + data1));
    } catch (SeedFormatException | CodecException | IOException e) {
      e.printStackTrace();
      fail();
    }

    ds.trimToCommonTime();
    AzimuthExperiment ae = new AzimuthExperiment();
    ae.runExperimentOnData(ds);
    double fitAngle = ae.getFitAngle();
    assertEquals(angle, fitAngle, 1.0);

  }

  @Test
  public void testSimpleCOWI() {
    String filename = folder + "cowi-multitests/C100823215422_COWI.LHx";
    String dataname1 = "US_COWI_  _LHN";
    String dataname2 = "US_COWI_  _LHE";
    String filename2 = folder + "cowi-multitests/DT000110.LH1";
    try {
      DataStore ds = new DataStore();
      ds.setBlock(0, filename, dataname1);
      ds.setBlock(1, filename, dataname2);
      ds.setBlock(2, filename2);

      String startString = "2010-236T02:00:00.0";
      String endString = "2010-236T13:00:00.0";
      long st = TestUtils.timeStringToEpochMilli(startString);
      long ed = TestUtils.timeStringToEpochMilli(endString);
      ds.trim(st, ed);

      AzimuthExperiment az = new AzimuthExperiment();
      az.setSimple(true);
      az.runExperimentOnData(ds);
      assertEquals(3.4, az.getFitAngle(), 0.5);

    } catch (SeedFormatException | CodecException | IOException e) {
      e.printStackTrace();
      fail();
    }
  }

  @Test
  public void alternateEntryPoint_testInitializationUnchanged_DefaultFalseSimpleCalc() {
    double[] north = new double[10000];
    Arrays.fill(north, 1);
    double[] east = new double[10000];
    Arrays.fill(east, 2);
    double[] referenceNorth = new double[10000];
    Arrays.fill(referenceNorth, 1);
    long interval = 1;
    long start = 0;
    long end = 9;
    AzimuthExperiment experiment = new AzimuthExperiment();
    experiment.alternateEntryPoint(north, east, referenceNorth, interval, start, end);

    assertFalse(experiment.getSimpleCalc());
    assertTrue(experiment.dataNames.contains("N"));
    assertTrue(experiment.dataNames.contains("E"));
    assertTrue(experiment.dataNames.contains("R"));

    assertNotNull(experiment.xySeriesData);
  }

  @Test
  public void alternateEntryPoint_testInitializationUnchangedTrue() {
    double[] north = new double[10000];
    Arrays.fill(north, 1);
    double[] east = new double[10000];
    Arrays.fill(east, 2);
    double[] referenceNorth = new double[10000];
    Arrays.fill(referenceNorth, 1);
    long interval = 1;
    long start = 0;
    long end = 9;
    AzimuthExperiment experiment = new AzimuthExperiment();
    experiment.setSimple(true);
    experiment.alternateEntryPoint(north, east, referenceNorth, interval, start, end);

    assertTrue(experiment.getSimpleCalc());
    assertTrue(experiment.dataNames.contains("N"));
    assertTrue(experiment.dataNames.contains("E"));
    assertTrue(experiment.dataNames.contains("R"));

    assertNotNull(experiment.xySeriesData);
  }

  @Test
  public void getAzimuth_simpleRotatedTest() {
    //create data 90 degrees off.
    double[] north = new double[10000];
    double[] east = new double[10000];
    double[] referenceNorth = new double[10000];
    Random rand = new Random();
    rand.setSeed(42);
    for (int i = 0; i < north.length; i++) {
      north[i] = rand.nextDouble();
      east[i] = rand.nextDouble();
      referenceNorth[i] = east[i];
    }
    long interval = 40;
    long start = 0;
    long end = 10000;
    assertEquals(
        4.71250648,
        AzimuthExperiment.getAzimuth(north, east, referenceNorth, interval, start, end),
        10E-7);
  }

  @Test
  public void getAzimuth_simpleUnrotatedTest() {
    //create data 90 degrees off.
    double[] north = new double[10000];
    double[] east = new double[10000];
    double[] referenceNorth = new double[10000];
    Random rand = new Random();
    rand.setSeed(43);
    for (int i = 0; i < north.length; i++) {
      north[i] = rand.nextDouble();
      east[i] = rand.nextDouble();
      referenceNorth[i] = north[i];
    }
    long interval = 40;
    long start = 0;
    long end = 10000;
    assertEquals(
        0.0,
        AzimuthExperiment.getAzimuth(north, east, referenceNorth, interval, start, end),
        10E-7);
  }
}

