package asl.sensor.experiment;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertTrue;
import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertNull;

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


}
