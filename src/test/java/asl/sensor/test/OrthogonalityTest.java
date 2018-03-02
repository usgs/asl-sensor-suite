package asl.sensor.test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.time.OffsetDateTime;
import java.util.Arrays;
import java.util.Calendar;
import java.util.TimeZone;
import org.junit.Before;
import org.junit.Test;
import asl.sensor.experiment.OrthogonalExperiment;
import asl.sensor.gui.InputPanel;
import asl.sensor.input.DataStore;
import edu.iris.dmc.seedcodec.CodecException;
import edu.sc.seis.seisFile.mseed.SeedFormatException;

public class OrthogonalityTest {

  public static String folder = TestUtils.DL_DEST_LOCATION + TestUtils.SUBPAGE;

  @Before
  public void getReferencedData() {

    // place in sprockets folder under 'from-sensor-test/[test-name]'
    String refSubfolder = TestUtils.SUBPAGE + "orthog-94/";
    String[] prefixes = new String[4];
    prefixes[0] = "00_LH1";
    prefixes[1] = "00_LH2";
    prefixes[2] = "10_LH1";
    prefixes[3] = "10_LH2";
    String extension = ".512.seed";

    for (String prefix : prefixes) {
      String fileID = prefix + extension;
      try {
        TestUtils.downloadTestData(refSubfolder, fileID, refSubfolder, fileID);
      } catch (IOException e) {
        e.printStackTrace();
      }
    }
  }

  public String getNoisyData() {
    return "FUNA.00";
  }

  public String getCleanData() {
    return "ANMO.10";
  }

  @Test
  public void identifiesSprocketsAngle060Clean() {
    testsFromSprockets(60, getCleanData());
  }

  @Test
  public void identifiesSprocketsAngle060Noisy() {
    testsFromSprockets(60, getNoisyData());
  }

  @Test
  public void identifiesSprocketsAngle080Clean() {
    testsFromSprockets(80, getCleanData());
  }

  @Test
  public void identifiesSprocketsAngle080Noisy() {
    testsFromSprockets(80, getNoisyData());
  }

  @Test
  public void identifiesSprocketsAngle085Clean() {
    testsFromSprockets(85, getCleanData());
  }

  @Test
  public void identifiesSprocketsAngle085Noisy() {
    testsFromSprockets(85, getNoisyData());
  }

  @Test
  public void identifiesSprocketsAngle087Clean() {
    testsFromSprockets(87, getCleanData());
  }

  @Test
  public void identifiesSprocketsAngle087Noisy() {
    testsFromSprockets(87, getNoisyData());
  }

  @Test
  public void identifiesSprocketsAngle088Clean() {
    testsFromSprockets(88, getCleanData());
  }

  @Test
  public void identifiesSprocketsAngle088Noisy() {
    testsFromSprockets(88, getNoisyData());
  }

  @Test
  public void identifiesSprocketsAngle089Clean() {
    testsFromSprockets(89, getCleanData());
  }

  @Test
  public void identifiesSprocketsAngle089Noisy() {
    testsFromSprockets(89, getNoisyData());
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
  public void identifiesSprocketsAngle091Clean() {
    testsFromSprockets(91, getCleanData());
  }

  @Test
  public void identifiesSprocketsAngle091Noisy() {
    testsFromSprockets(91, getNoisyData());
  }

  @Test
  public void identifiesSprocketsAngle092Clean() {
    testsFromSprockets(92, getCleanData());
  }

  @Test
  public void identifiesSprocketsAngle092Noisy() {
    testsFromSprockets(92, getNoisyData());
  }

  @Test
  public void identifiesSprocketsAngle100Clean() {
    testsFromSprockets(100, getCleanData());
  }

  @Test
  public void identifiesSprocketsAngle100Noisy() {
    testsFromSprockets(100, getNoisyData());
  }

  @Test
  public void identifiesSprocketsAngle120Clean() {
    testsFromSprockets(120, getCleanData());
  }

  @Test
  public void identifiesSprocketsAngle120Noisy() {
    testsFromSprockets(120, getNoisyData());
  }

  @Test
  public void identifiesSprocketsAngle270Clean() {
    testsFromSprockets(270, getCleanData());
  }

  @Test
  public void identifiesSprocketsAngle270Noisy() {
    testsFromSprockets(270, getNoisyData());
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
    String testURL = refURL + "orthogonality/" + sb.toString() + "/";

    String refSubfolder = "orientation-reference/";
    String testSubfolder = "orthogonality-" + sb.toString() + "/";
    try {
      // get reference data if needed
      TestUtils.downloadTestData(refURL, data1, refSubfolder, data1);
      TestUtils.downloadTestData(refURL, data2, refSubfolder, data2);
      // get test data if needed
      TestUtils.downloadTestData(testURL, data1, testSubfolder, data1);
      TestUtils.downloadTestData(testURL, data2, testSubfolder, data2);
    } catch (IOException e) {
      e.printStackTrace();
      fail();
    }

    DataStore ds = new DataStore();
    try {
      String root = TestUtils.DL_DEST_LOCATION;
      testSubfolder = root + testSubfolder;
      refSubfolder = root + refSubfolder;

      ds.setBlock(0, refSubfolder + data1);
      ds.setBlock(1, refSubfolder + data2);

      ds.setBlock(2, testSubfolder + data1);
      ds.setBlock(3, testSubfolder + data2);

      OffsetDateTime cCal = TestUtils.getStartCalendar(ds);
      cCal = cCal.withHour(10);
      cCal = cCal.withMinute(0);
      cCal = cCal.withSecond(0);
      cCal = cCal.withNano(0);
      long start = cCal.toInstant().toEpochMilli();
      cCal = cCal.withHour(14);
      long end = cCal.toInstant().toEpochMilli();

      ds.trim(start, end);

      // ds.trimToCommonTime();
      OrthogonalExperiment oe = new OrthogonalExperiment();
      oe.runExperimentOnData(ds);
      double fitAngle = oe.getFitAngle();
      System.out.println( Arrays.toString(oe.getSolutionParams()) );

      double expectedAngle = angle;
      if (expectedAngle > 180) {
        expectedAngle = 360 - expectedAngle;
      }
      System.out.println(expectedAngle + " | " + fitAngle + " [" + sb.toString() + "]");

      assertEquals(expectedAngle, fitAngle, 1.0);

    } catch (SeedFormatException | CodecException e) {
      e.printStackTrace();
      fail();
    }

  }

  @Test
  public void getsCorrectAngle() {

    DataStore ds = new DataStore();
    String testFolder = folder + "orthog-94/";
    String[] prefixes = new String[4];
    prefixes[0] = "00_LH1";
    prefixes[1] = "00_LH2";
    prefixes[2] = "10_LH1";
    prefixes[3] = "10_LH2";
    String extension = ".512.seed";

    for (int i = 0; i < prefixes.length; ++i) {
      String fName = testFolder + prefixes[i] + extension;
      try {
        ds.setBlock(i, fName);
      } catch (SeedFormatException | CodecException e) {
        e.printStackTrace();
        fail();
      }
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
