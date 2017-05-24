package asl.sensor.experiment;

/**
 * Enumerated type defining each kind of test, for GUI and test factory
 * @author akearns
 *
 */
public enum ExperimentEnum {

  // if adding a new test to this, make sure to also create a new extender for
  // experiment and a corresponding constructor call in experimentfactory
  
  NOISE ("Self-noise"),
  RGAIN ("Relative gain"),
  STCAL ("Step calibration"),
  RANDM ("Randomized calibration"),
  RNDBG ("R-Cal. (no solver)"),
  AZMTH ("Azimuth"),
  ORTHO ("Orthogonality"),
  RESPN ("Response"),
  NOIS9 ("Self-noise (9-Component)");

  private String name;
  
  ExperimentEnum(String name) {
    this.name = name;  }
  
  
  /**
   * Get the full name of this experiment (used for plot & tab names)
   * @return Name of experiment, as String
   */
  public String getName() {
    return name;
  }
  
}
