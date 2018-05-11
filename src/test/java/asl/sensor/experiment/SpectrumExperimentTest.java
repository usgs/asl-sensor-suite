package asl.sensor.experiment;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
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
  public void backend_standardData() throws SeedFormatException, UnsupportedCompressionType,
  CodecException, IOException {
    //There doesn't appear to be any "good results" to compare against since it primarily plots.
    //We can at least check the below
    SpectrumExperiment experiment = new SpectrumExperiment();
    DataStore ds = new DataStore();
    String resp = TestUtils.RESP_LOCATION + "T-compact_Q330HR_BH_40";
    String testData1 = folder + "noise-neg159db/" + "00_BH0.512.seed";
    ds.setBlock(0, testData1);
    long start = ds.getCommonTime().getFirst();
    long end = ds.getCommonTime().getSecond();
    end = (end - start) / 4 + start; // trim data to 1/4 its current length to speed up cal
    ds.trim(start, end);
    ds.setResponse(0, resp);
    experiment.runExperimentOnData(ds);
    //assert datanames was correctly populated
    assertEquals(2, experiment.dataNames.size());
    assertEquals("XX_TST5_00_BH0", experiment.dataNames.get(0));
    assertEquals("T-compact_Q330HR_BH_40", experiment.dataNames.get(1));
    //assert respIndices was correctly populated
    assertEquals(0, experiment.listActiveResponseIndices()[0]);
    //assert xySeriesData was populated
    assertEquals(1, experiment.getData().size());
    String[] keycheck = {"PSD XX_TST5_00_BH0 [0]", "NLNM", "NHNM"};
    for (int i = 0; i < experiment.getData().get(0).getSeriesCount(); ++i) {
      String key = (String) experiment.getData().get(0).getSeriesKey(i);
      assertEquals(keycheck[i], key);
    }
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