package asl.sensor.experiment;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertTrue;
import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertNull;

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
  public void runExperimentOnData_doesItReinitializeFields() {
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


}
