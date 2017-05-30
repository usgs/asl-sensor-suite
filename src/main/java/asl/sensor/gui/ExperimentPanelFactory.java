package asl.sensor.gui;

import asl.sensor.experiment.ExperimentEnum;

/**
 * Simple factory method for creating an experiment panel
 * @author akearns
 *
 */
public class ExperimentPanelFactory {

  /**
   * Instantiate an ExperimentPanel based on the enumerated type passed in
   * This enumerated type is then used to instantiate a backend using
   * the ExperimentFactory class.
   * @param exp Type of experiment (see the enum for details on kinds)
   * @return A concrete implementation of the ExperimentPanel abstract class
   */
  public static ExperimentPanel createPanel(
      ExperimentEnum exp){
    switch(exp){
    case ORTHO:
      return new OrthogonalPanel(exp);
    case NOIS9:
      return new NoiseNinePanel(exp);
    case NOISE:
      return new NoisePanel(exp);
    case RGAIN:
      return new GainPanel(exp);
    case GAIN6:
      return new GainSixPanel(exp);
    case STCAL:
      return new StepPanel(exp);
    case AZMTH:
      return new AzimuthPanel(exp);
    case RANDM:
      return new RandomizedPanel(exp);
    case RESPN:
      return new ResponsePanel(exp);
    case RNDBG:
      return new RandDebugPanel(exp);
    default:
      // this shouldn't happen unless someone added to the enum
      // and forgot to follow-through on implementation
      throw new IllegalArgumentException("Invalid enum type specified");
    }
  }
  
}
