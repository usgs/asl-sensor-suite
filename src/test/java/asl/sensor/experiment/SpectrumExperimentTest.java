package asl.sensor.experiment;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;
import java.io.IOException;
import org.junit.Test;
import asl.sensor.input.DataStore;
import asl.sensor.test.TestUtils;
import edu.iris.dmc.seedcodec.CodecException;
import edu.iris.dmc.seedcodec.UnsupportedCompressionType;
import edu.sc.seis.seisFile.mseed.SeedFormatException;

public class SpectrumExperimentTest {

  public static String folder = TestUtils.TEST_DATA_LOCATION + TestUtils.SUBPAGE;

  @Test
  public void spectrumExperiment_constructorInitializes() {
    SpectrumExperiment experiment = new SpectrumExperiment();

    assertFalse(experiment.getFreqSpace());
    assertArrayEquals(new int[3], experiment.listActiveResponseIndices());

    //This checks if super() was called as well
    assertTrue(experiment.dataNames.isEmpty());
  }

  @Test
  public void backend_standardData() {
    //There doesn't appear to be any "good results" to compare against since it primarily plots.
    //We can at least check the below

    //assert datanames was correctly populated

    //assert respIndices was correctly populated

    //assert xySeriesData was populated

    fail();
  }

  @Test
  public void blocksNeeded() {
    SpectrumExperiment experiment = new SpectrumExperiment();
    assertEquals(3, experiment.blocksNeeded());
  }

  @Test
  public void hasEnoughData_notEnoughData() throws SeedFormatException, UnsupportedCompressionType,
  CodecException, IOException {
    SpectrumExperiment experiment = new SpectrumExperiment();
    DataStore ds = new DataStore();
    String resp = TestUtils.RESP_LOCATION + "T-compact_Q330HR_BH_40";
    String testData1 = folder + "noise-neg159db/" + "00_BH0.512.seed";
    String testData2 = folder + "noise-neg159db/" + "TST6.00_BH0.512.seed";
    ds.setBlock(0, testData1);
    ds.setResponse(1, resp);
    ds.setBlock(2, testData2);
    // fails when all data within range (i.e., indices 0-2) has only one of either resp or seed data
    // (has enough data when at least 1 has both resp and seed)
    assertFalse(experiment.hasEnoughData(ds));
  }
  @Test
  public void hasEnoughData_allNeededDataPresent() throws SeedFormatException,
  UnsupportedCompressionType, CodecException, IOException {
    SpectrumExperiment experiment = new SpectrumExperiment();
    DataStore ds = new DataStore();
    String resp = TestUtils.RESP_LOCATION + "T-compact_Q330HR_BH_40";
    String testData1 = folder + "noise-neg159db/" + "00_BH0.512.seed";
    String testData2 = folder + "noise-neg159db/" + "TST6.00_BH0.512.seed";
    ds.setBlock(0, testData1);
    ds.setResponse(0, resp);
    ds.setResponse(1, resp);
    ds.setBlock(2, testData2);
    //Verify it returns true when it should
    assertTrue(experiment.hasEnoughData(ds));
  }

  @Test
  public void setFreqSpace() {
    SpectrumExperiment experiment = new SpectrumExperiment();
    assertFalse(experiment.getFreqSpace());

    //Verify not toggling
    experiment.setFreqSpace(false);
    assertFalse(experiment.getFreqSpace());

    //Verify set actually sets
    experiment.setFreqSpace(true);
    assertTrue(experiment.getFreqSpace());

    //Verify isn't a one way operation
    experiment.setFreqSpace(false);
    assertFalse(experiment.getFreqSpace());
  }
}