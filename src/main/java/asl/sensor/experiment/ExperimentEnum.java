package asl.sensor.experiment;

/**
 * Enumerated type defining each kind of test, for GUI and test factory
 * @author akearns
 *
 */
public enum ExperimentEnum {

  // if adding a new test to this, make sure to also create a new extender for
  // experiment and a corresponding constructor call in experimentfactory
  // due to how iterating through an enum works, the order in which panels'
  // tabs appear in the GUI should match up with the 
  NOISE ("Self-noise"),
  NOIS9 ("Self-noise (3-component)"),
  RGAIN ("Relative gain"),
  GAIN6 ("Relative gain (3-component)"),
  STCAL ("Step calibration"),
  RANDM ("Randomized calibration"),
  RNDBG ("Rdm. cal. verification"),
  AZMTH ("Azimuth"),
  ORTHO ("Orthogonality"),
  RESPN ("Response");

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
