package asl.sensor.gui;

import asl.sensor.ExperimentFactory;

/**
 * Simple factory method for creating an experiment panel
 *
 * @author akearns
 */
public class ExperimentPanelFactory {

  /**
   * Instantiate an ExperimentPanel based on the enumerated type passed in
   * This enumerated type is then used to instantiate a backend.
   *
   * @param exp Type of experiment (see the enum for details on kinds)
   * @return A concrete implementation of the ExperimentPanel abstract class
   */
  public static ExperimentPanel createPanel(
      ExperimentFactory exp) {
    switch (exp) {
      case ORTHOGONALITY:
        return new OrthogonalPanel(exp);
      case NOISE9:
        return new NoiseNinePanel(exp);
      case NOISE:
        return new NoisePanel(exp);
      case GAIN:
        return new GainPanel(exp);
      case GAIN6:
        return new GainSixPanel(exp);
      case STEPCAL:
        return new StepPanel(exp);
      case AZIMUTH:
        return new AzimuthPanel(exp);
      case RANDOMCAL:
        return new RandomizedPanel(exp);
      case RESPONSE:
        return new ResponsePanel(exp);
      case SINECAL:
        return new SinePanel(exp);
      case SPECTRUM:
        return new SpectrumPanel(exp);
      default:
        // this shouldn't happen unless someone added to the enum
        // and forgot to follow-through on implementation
        throw new IllegalArgumentException("Invalid enum type specified");
    }
  }

}
