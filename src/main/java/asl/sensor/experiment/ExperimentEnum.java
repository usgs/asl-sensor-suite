package asl.sensor.experiment;

/**
 * Enumerated type defining each kind of test, for GUI and test factory
 * @author akearns
 *
 */
public enum ExperimentEnum {

  // if adding a new test to this, make sure to also create a new extender for
  // experiment and a corresponding constructor call in experimentfactory
  
  NOISE ("Self-noise", 3, 3),
  RGAIN ("Relative gain", 2, 2),
  STCAL ("Step calibration", 2, 1),
  AZMTH ("Azimuth", 3, 0),
  ORTHO ("Orthogonality", 4, 0);

  private String name;
  private int blocksNeeded;
  private int fullDataNeeded;
    // not my favorite way of doing this but it's an easy modification to make
  
  ExperimentEnum(String name, int blocksNeeded, int fullDataNeeded) {
    this.name = name;
    this.blocksNeeded = blocksNeeded;
    this.fullDataNeeded = fullDataNeeded;
  }
  
  /**
   * Get the amount of data blocks needed for the calculation
   * @return
   */
  public int blocksNeeded() {
    return blocksNeeded;
  }
  
  /**
   * Get the amount of data and response pairs needed for the calculation
   * @return
   */
  public int fullDataNeeded() {
    return fullDataNeeded;
  }
  
  /**
   * Get the full name of this experiment (used for plot & tab names)
   * @return
   */
  public String getName() {
    return name;
  }
  
}
