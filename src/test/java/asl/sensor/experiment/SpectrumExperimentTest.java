package asl.sensor.experiment;

import static asl.sensor.test.TestUtils.RESP_LOCATION;
import static asl.sensor.test.TestUtils.getSeedFolder;
import static org.junit.Assert.*;

import asl.sensor.input.DataStore;
import asl.sensor.input.DataStoreUtils;
import org.junit.Test;

public class SpectrumExperimentTest {

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
  public void hasEnoughData_notEnoughData() {
    //Verify this fails when it should fail.
    //Does it need a response and data for all 3 indices?
    //What if 1 index is good and the others have data, but are missing resps?
    //What if 1 index is good and the others have resps, but are missing data?
    assertFalse(true);
  }
  @Test
  public void hasEnoughData_allNeededDataPresent() {
    //Verify it returns true when it should
    assertTrue(false);
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