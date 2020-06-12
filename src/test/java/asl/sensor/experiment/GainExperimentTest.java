package asl.sensor.experiment;

import static java.lang.Math.abs;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import asl.sensor.input.DataStore;
import asl.sensor.test.TestUtils;
import asl.utils.input.InstrumentResponse;
import edu.iris.dmc.seedcodec.CodecException;
import edu.sc.seis.seisFile.mseed.SeedFormatException;
import java.io.IOException;
import java.time.OffsetDateTime;
import java.time.ZoneOffset;
import org.apache.commons.math3.complex.Complex;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.junit.Test;

public class GainExperimentTest {

  private static final String folder = TestUtils.TEST_DATA_LOCATION + TestUtils.SUBPAGE;

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
    rnames[1] = "RESP.IU.ANMO.10.BHZ";

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

    GainExperiment ge = new GainExperiment();
    ge.runExperimentOnData(ds);

    double[] stats = ge.getStatsFromPeak(0);
    double gain = stats[3];
    // System.out.println( Arrays.toString(stats) );
    assertEquals(11719., gain, 2.0);
  }

  @Test
  public void testA0Calculation() throws IOException, SeedFormatException, CodecException {
    // attempts to replicate the logic of the A0 calculation in the backend
    String testFolder = folder + "relative-gain-100/";
    String seedName = testFolder + "00_BHZ.512.seed";
    DataStore ds = new DataStore();
    ds.setBlock(0, seedName);
    ds.setBlock(1, seedName);

    String filename = "RESP.IU.ANMO.00.BHZ_gainx100";
    InstrumentResponse ir = new InstrumentResponse(testFolder + filename);
    double initialA0 = ir.getNormalization();
    ds.setResponse(0, ir);
    ir = new InstrumentResponse(testFolder + filename);
    ir.setNormalization(1);
    ds.setResponse(1, ir);

    OffsetDateTime start =
        OffsetDateTime.ofInstant(ds.getBlock(0).getStartInstant(), ZoneOffset.UTC);
    start = start.withHour(10);
    OffsetDateTime end = start.withHour(12);
    ds.trim(start.toInstant(), end.toInstant());

    GainExperiment ge = new GainExperiment();
    ge.runExperimentOnData(ds);

    XYSeriesCollection xys = ge.getData().get(0);
    XYSeries ref = xys.getSeries(0);
    XYSeries test = xys.getSeries(1);
    for (int i = 0; i < ref.getItemCount(); ++i) {
      assertEquals((double) ref.getY(i), (double) test.getY(i), 1E-3);
    }
    double[] resultValues = ge.getStatsFromFreqs();
    double referenceGain = resultValues[2];
    double calculatedGain = resultValues[3];
    assertEquals(referenceGain, calculatedGain, 1E-1);
    double testA0 = resultValues[9];
    assertEquals(initialA0, testA0, 1E-1);
  }
}
