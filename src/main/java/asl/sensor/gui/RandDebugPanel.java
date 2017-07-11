package asl.sensor.gui;

import asl.sensor.experiment.ExperimentEnum;
import asl.sensor.experiment.RandDebugExperiment;

public class RandDebugPanel extends RandomizedPanel {

  /**
   * 
   */
  private static final long serialVersionUID = -4314038309895842615L;

  public RandDebugPanel(ExperimentEnum exp) {
    super(exp);
    expResult = new RandDebugExperiment();
  }

}
