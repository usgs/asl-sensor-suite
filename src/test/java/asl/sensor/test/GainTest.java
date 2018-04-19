package asl.sensor.test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;
import java.io.IOException;
import java.time.OffsetDateTime;
import java.time.ZoneOffset;
import org.junit.Test;
import asl.sensor.experiment.GainExperiment;
import asl.sensor.input.DataStore;
import edu.iris.dmc.seedcodec.CodecException;
import edu.sc.seis.seisFile.mseed.SeedFormatException;

public class GainTest {

  public static String folder = TestUtils.TEST_DATA_LOCATION + TestUtils.SUBPAGE;

  @Test
  public void testGainCalculation() {

    DataStore ds = new DataStore();

    String testFolder = folder + "relative-gain-100/";
    String[] prefixes = new String[2];
    prefixes[0] = "00_BHZ";
    prefixes[1] = "10_BHZ";
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

    String[] rnames = new String[2];
    rnames[0] = "RESP.IU.ANMO.00.BHZ_gainx100";
    //rnames[0] = "RESP.IU.ANMO.00.BHZ";
    rnames[1] = "RESP.IU.ANMO.10.BHZ";

    for (int i = 0; i < rnames.length; ++i) {
      String fName = testFolder + rnames[i];
      try {
        ds.setResponse(i, fName);
      } catch (IOException e) {
        // TODO Auto-generated catch block
        e.printStackTrace();
        fail();
      }
    }

    OffsetDateTime start =
        OffsetDateTime.ofInstant(ds.getBlock(0).getStartInstant(), ZoneOffset.UTC);
    start = start.withHour(10);
    OffsetDateTime end = start.withHour(14);

    ds.trim(start.toInstant(), end.toInstant());

    GainExperiment ge = new GainExperiment();
    ge.runExperimentOnData(ds);

    double[] stats = ge.getStatsFromPeak(0);
    double gain = stats[3];
    // System.out.println( Arrays.toString(stats) );
    assertEquals(11714., gain, 2.0);

  }
}
