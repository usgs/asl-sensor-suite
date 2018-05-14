package asl.sensor.experiment;

/**
 * Enumerated type defining each kind of test, done so GUI has list of all experiments available
 * and for creating the associated Experiment class.
 *
 * @author akearns - KBRWyle
 */
public enum ExperimentEnum {

  /**
   * If adding a new test, make sure to also create a new extension for Experiment, associated
   * ExperimentPanel extension, and entry in ExperimentPanelFactory.
   * Due to how iterating through an enum works, the order in which panel tabs appear in the
   * GUI should match up with the order they are listed here.
   **/
  NOISE("Self-noise") {
    @Override
    public Experiment createExperiment() {
      return new NoiseExperiment();
    }
  },
  NOIS9("Self-noise (3-component)") {
    @Override
    public Experiment createExperiment() {
      return new NoiseNineExperiment();
    }
  },
  RGAIN("Relative gain") {
    @Override
    public Experiment createExperiment() {
      return new GainExperiment();
    }
  },
  GAIN6("Relative gain (3-component)") {
    @Override
    public Experiment createExperiment() {
      return new GainSixExperiment();
    }
  },
  STCAL("Step calibration") {
    @Override
    public Experiment createExperiment() {
      return new StepExperiment();
    }
  },
  RANDM("Randomized calibration") {
    @Override
    public Experiment createExperiment() {
      return new RandomizedExperiment();
    }
  },
  SINCL("Sine calibration") {
    @Override
    public Experiment createExperiment() {
      return new SineExperiment();
    }
  },
  AZMTH("Azimuth") {
    @Override
    public Experiment createExperiment() {
      return new AzimuthExperiment();
    }
  },
  ORTHO("Orthogonality") {
    @Override
    public Experiment createExperiment() {
      return new OrthogonalExperiment();
    }
  },
  SPECT("Power-spectrum") {
    @Override
    public Experiment createExperiment() {
      return new SpectrumExperiment();
    }
  },
  RESPN("Response") {
    @Override
    public Experiment createExperiment() {
      return new ResponseExperiment();
    }
  };
  private String name;

  ExperimentEnum(String name) {
    this.name = name;
  }

  public abstract Experiment createExperiment();

  /**
   * Get the full name of this experiment (used for plot & tab names)
   *
   * @return Name of experiment, as String
   */
  public String getName() {
    return name;
  }

}
