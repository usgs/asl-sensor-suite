package asl.sensor.experiment;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertTrue;
import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertNull;

import asl.sensor.input.DataStore;
import java.util.ArrayList;
import org.junit.Test;

public class ExperimentTest {

  @Test
  public void experiment_constructorInitializes(){
    Experiment experiment = new MockExperiment();
    assertEquals(0, experiment.getStart());
    assertEquals(0, experiment.getEnd());
    assertTrue(experiment.getInputNames().isEmpty());
    assertEquals("", experiment.getStatus());
    assertNull(experiment.getGapRegions());
  }

  @Test
  public void fireStateChange_updatesStatus(){
    MockExperiment experiment = new MockExperiment();
    experiment.fireStateChange("Fired Status Change");
    assertEquals("Fired Status Change", experiment.getStatus());

    assertEquals(1, experiment.numberOfChangesFired);
  }

  @Test
  public void listActiveResponseIndices_defaultImplementation(){
    Experiment experiment = new MockExperiment();
    assertArrayEquals(new int[]{}, experiment.listActiveResponseIndices());
  }

  @Test
  public void runExperimentOnData_doesItCheckForDataAndBlocksThenShortCircuit_doesItCallBackend() {
    MockExperiment experiment = new MockExperiment();
    experiment.runExperimentOnData(null);

    assertTrue(experiment.hasEnoughDataCalled);
    assertTrue(experiment.backendCalled);
    assertTrue(experiment.blocksNeededCalled);

    assertEquals(1, experiment.numberOfChangesFired);
  }

  @Test
  public void runExperimentOnData_doesItReinitializeFields_shortCircuited() {
    Experiment experiment = new MockExperiment();
    //Dirty everything
    experiment.start = 1L;
    experiment.end = 2L;
    experiment.dataNames.add("Not Empty");
    //experiment.xySeriesData defaults to null
    //experiment.getGapRegions defaults to null;

    experiment.runExperimentOnData(null);

    //These should have been reinitialized
    assertEquals(0L, experiment.getStart());
    assertEquals(0L, experiment.getEnd());
    assertTrue(experiment.getInputNames().isEmpty());
    assertTrue(experiment.getGapRegions().isEmpty());
  }

  @Test(expected = IndexOutOfBoundsException.class)
  public void runExperimentOnData_notEnoughData_throwsException() {
    MockExperiment experiment = new MockExperiment();
    //Dirty everything
    experiment.start = 1L;
    experiment.end = 2L;
    experiment.dataNames.add("Not Empty");
    //experiment.xySeriesData defaults to null
    //experiment.getGapRegions defaults to null;

    //Set it to not shortCircuit, but still not try loading
    experiment.setHasEnoughData = false;
    try {
      experiment.runExperimentOnData(new DataStore());
    } catch (IndexOutOfBoundsException e) {
      //Should be called still
      assertTrue(experiment.hasEnoughDataCalled);

      assertEquals(1, experiment.numberOfChangesFired);

      //These should have been untouched
      assertEquals(1L, experiment.getStart());
      assertEquals(2L, experiment.getEnd());

      //These should have been reinitialized still
      assertTrue(experiment.getInputNames().isEmpty());
      assertTrue(experiment.getGapRegions().isEmpty());
      throw e;
    }
  }


}
