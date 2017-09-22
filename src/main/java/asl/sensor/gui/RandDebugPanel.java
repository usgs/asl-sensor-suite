package asl.sensor.gui;

import asl.sensor.experiment.ExperimentEnum;
import asl.sensor.experiment.RandDebugExperiment;

/**
 * Used to determine possible flaws in response choice or data deconvolution
 * that may prevent the solver for the response parameters from converging in
 * a randomized calibration experiment.
 * @author akearns
 *
 */
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
