package asl.sensor.experiment;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;

import asl.sensor.input.DataStore;
import asl.sensor.test.TestUtils;
import edu.iris.dmc.seedcodec.CodecException;
import edu.sc.seis.seisFile.mseed.SeedFormatException;
import java.io.IOException;
import java.time.OffsetDateTime;
import java.time.ZoneOffset;
import org.junit.Test;

public class GainSixExperimentTest {

  private static final String folder = TestUtils.TEST_DATA_LOCATION + TestUtils.SUBPAGE;
  private final double EXPECTED_GAIN = 1143;

  @Test
  public void testGainCalculationFirstReference() {

    DataStore ds = new DataStore();

    String testFolder = folder + "relative-gain-sixin/";
    String[] prefixes = new String[6];
    prefixes[0] = "00_BH1";
    prefixes[1] = "00_BH2";
    prefixes[2] = "00_BHZ";
    prefixes[3] = "10_BH1";
    prefixes[4] = "10_BH2";
    prefixes[5] = "10_BHZ";
    String extension = ".512.seed";

    for (int i = 0; i < prefixes.length; ++i) {
      String fName = testFolder + prefixes[i] + extension;
      try {
        ds.setBlock(i, fName);
      } catch (IOException | SeedFormatException | CodecException e) {
        e.printStackTrace();
        fail();
      }
    }

    String[] rnames = new String[6];
    rnames[0] = "RESP.IU.ANMO.00.BH1";
    rnames[1] = "RESP.IU.ANMO.00.BH2";
    rnames[2] = "RESP.IU.ANMO.00.BHZ";
    rnames[3] = "RESP.IU.ANMO.10.BH1";
    rnames[4] = "RESP.IU.ANMO.10.BH2";
    rnames[5] = "RESP.IU.ANMO.10.BHZ";

    for (int i = 0; i < rnames.length; ++i) {
      String fName = testFolder + rnames[i];
      try {
        ds.setResponse(i, fName);
      } catch (IOException e) {
        fail(e.getMessage());
      }
    }

    OffsetDateTime start =
        OffsetDateTime.ofInstant(ds.getBlock(0).getStartInstant(), ZoneOffset.UTC);
    start = start.withHour(10);
    OffsetDateTime end = start.withHour(14);

    ds.trim(start.toInstant(), end.toInstant());

    GainSixExperiment ge = new GainSixExperiment();
    ge.runExperimentOnData(ds);

    ge.setRangeForStatistics(
        GainExperiment.DEFAULT_LOW_BOUND, GainExperiment.DEFAULT_UP_BOUND);
    double[][] allStats = ge.getStatistics();
    double[] stats = allStats[0];
    double gain = stats[3];
    // System.out.println( Arrays.toString(stats) );
    assertEquals("ACTUAL GAIN VALUE WAS " + gain, EXPECTED_GAIN, gain, 10.);
  }

  @Test
  public void testGainCalculationSecondReference() {

    DataStore ds = new DataStore();

    String testFolder = folder + "relative-gain-sixin/";
    String[] prefixes = new String[6];
    prefixes[0] = "10_BH1";
    prefixes[1] = "10_BH2";
    prefixes[2] = "10_BHZ";
    prefixes[3] = "00_BH1";
    prefixes[4] = "00_BH2";
    prefixes[5] = "00_BHZ";
    String extension = ".512.seed";

    for (int i = 0; i < prefixes.length; ++i) {
      String fName = testFolder + prefixes[i] + extension;
      try {
        ds.setBlock(i, fName);
      } catch (IOException | SeedFormatException | CodecException e) {
        e.printStackTrace();
        fail();
      }
    }

    String[] rnames = new String[6];
    rnames[0] = "RESP.IU.ANMO.10.BH1";
    rnames[1] = "RESP.IU.ANMO.10.BH2";
    rnames[2] = "RESP.IU.ANMO.10.BHZ";
    rnames[3] = "RESP.IU.ANMO.00.BH1";
    rnames[4] = "RESP.IU.ANMO.00.BH2";
    rnames[5] = "RESP.IU.ANMO.00.BHZ";

    for (int i = 0; i < rnames.length; ++i) {
      String fName = testFolder + rnames[i];
      try {
        ds.setResponse(i, fName);
      } catch (IOException e) {
        fail(e.getMessage());
      }
    }

    OffsetDateTime start =
        OffsetDateTime.ofInstant(ds.getBlock(3).getStartInstant(), ZoneOffset.UTC);
    start = start.withHour(10);
    OffsetDateTime end = start.withHour(14);

    ds.trim(start.toInstant(), end.toInstant());

    GainSixExperiment ge = new GainSixExperiment();
    ge.setSecondDataAsAngleReference();
    ge.runExperimentOnData(ds);
    ge.setReferenceIndex(1);

    ge.setRangeForStatistics(
        GainExperiment.DEFAULT_LOW_BOUND, GainExperiment.DEFAULT_UP_BOUND);
    double[][] allStats = ge.getStatistics();
    double[] stats = allStats[0];
    double gain = stats[3];
    assertEquals("ACTUAL GAIN VALUE WAS " + gain, EXPECTED_GAIN, gain, 10.);
  }

}
