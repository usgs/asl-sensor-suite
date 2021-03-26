package asl.sensor.experiment;

import static asl.utils.response.ResponseBuilders.deepCopyWithNewPolesZeros;
import static asl.utils.response.ResponseParser.loadEmbeddedResponse;
import static asl.utils.response.ResponseParser.parseResponse;
import static junit.framework.TestCase.assertFalse;
import static junit.framework.TestCase.assertTrue;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;

import asl.sensor.input.DataStore;
import asl.sensor.test.TestUtils;
import asl.utils.response.ChannelMetadata;
import asl.utils.response.ChannelMetadata.ResponseStageException;
import asl.utils.response.PolesZeros;
import asl.utils.response.ResponseBuilders.PolesZerosBuilder;
import asl.utils.response.ResponseUnits.ResolutionType;
import asl.utils.response.ResponseUnits.SensorType;
import edu.iris.dmc.seedcodec.CodecException;
import edu.sc.seis.seisFile.mseed.SeedFormatException;
import java.io.IOException;
import java.time.OffsetDateTime;
import java.time.ZoneOffset;
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
    ge.setRangeForStatistics(GainExperiment.DEFAULT_LOW_BOUND, GainExperiment.DEFAULT_UP_BOUND);
    double[] stats = ge.getStatsFromFreqs();
    double gain = stats[3];
    assertEquals(11720., gain, 2.0);
  }

  @Test
  public void testEmbedCheck_trueForOneEmbed() throws IOException {

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

    String fName = testFolder + "RESP.IU.ANMO.00.BHZ_gainx100";
    ds.setResponse(0, fName);
    ds.setResponse(1, loadEmbeddedResponse(SensorType.STS1T5, ResolutionType.HIGH));

    OffsetDateTime start =
        OffsetDateTime.ofInstant(ds.getBlock(0).getStartInstant(), ZoneOffset.UTC);
    start = start.withHour(10);
    OffsetDateTime end = start.withHour(14);

    ds.trim(start.toInstant(), end.toInstant());

    GainExperiment ge = new GainExperiment();
    ge.runExperimentOnData(ds);
    assertTrue(ge.skipsFIRStages());
  }

  @Test
  public void testEmbedCheck_falseForTwoEmbeds() throws IOException {

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

    ds.setResponse(0, loadEmbeddedResponse(SensorType.STS1T5, ResolutionType.HIGH));
    ds.setResponse(1, loadEmbeddedResponse(SensorType.STS1T5, ResolutionType.HIGH));

    OffsetDateTime start =
        OffsetDateTime.ofInstant(ds.getBlock(0).getStartInstant(), ZoneOffset.UTC);
    start = start.withHour(10);
    OffsetDateTime end = start.withHour(14);

    ds.trim(start.toInstant(), end.toInstant());

    GainExperiment ge = new GainExperiment();
    ge.runExperimentOnData(ds);
    assertFalse(ge.skipsFIRStages());
  }

  @Test
  public void testA0Calculation()
      throws IOException, SeedFormatException, CodecException, ResponseStageException {
    // attempts to replicate the logic of the A0 calculation in the backend
    String testFolder = folder + "relative-gain-100/";
    String seedName = testFolder + "00_BHZ.512.seed";
    DataStore ds = new DataStore();
    ds.setBlock(0, seedName);
    ds.setBlock(1, seedName);

    String filename = "RESP.IU.ANMO.00.BHZ_gainx100";
    ChannelMetadata ir = parseResponse(testFolder + filename);
    double initialA0 = ir.getPoleZeroStage().getNormalizationFactor();
    ds.setResponse(0, ir);
    ir = parseResponse(testFolder + filename);
    PolesZeros original = ir.getPoleZeroStage();
    PolesZerosBuilder builder = new PolesZerosBuilder()
        .withPoles(original.getPoleDoubleList())
        .withZeros(original.getZeroDoubleList())
        .withInputUnits(original.getInputMotionUnit())
        .withOutputUnits(original.getOutputUnits())
        .withTransferFunction(original.getTransferFunction())
        .withNormalization(1)
        .withNormalizationFrequency(original.getNormalizationFreq());
    ChannelMetadata replacedNormalization = deepCopyWithNewPolesZeros(ir, builder.build());
    ds.setResponse(1, replacedNormalization);

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
