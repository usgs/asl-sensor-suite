package asl.sensor;

import asl.sensor.experiment.AzimuthExperiment;
import asl.sensor.experiment.Experiment;
import asl.sensor.experiment.GainExperiment;
import asl.sensor.experiment.GainSixExperiment;
import asl.sensor.experiment.NoiseExperiment;
import asl.sensor.experiment.NoiseNineExperiment;
import asl.sensor.experiment.OrthogonalExperiment;
import asl.sensor.experiment.RandomizedExperiment;
import asl.sensor.experiment.ResponseExperiment;
import asl.sensor.experiment.SineExperiment;
import asl.sensor.experiment.SpectrumExperiment;
import asl.sensor.experiment.StepExperiment;

/**
 * Enumerated type defining each kind of test, done so GUI has list of all experiments available
 * and for creating the associated Experiment class.
 *
 * If adding a new test, make sure to also create a new extension for Experiment, associated
 * ExperimentPanel extension.
 * Due to how iterating through an enum works, the order in which panel tabs appear in the
 * GUI should match up with the order they are listed here.
 *
 * @author akearns - KBRWyle
 * @author jholland - USGS
 */
public enum ExperimentFactory {

  NOISE("Self-noise") {
    @Override
    public Experiment createExperiment() {
      return new NoiseExperiment();
    }
  },
  NOISE9("Self-noise (3-component)") {
    @Override
    public Experiment createExperiment() {
      return new NoiseNineExperiment();
    }
  },
  GAIN("Relative gain") {
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
  STEPCAL("Step calibration") {
    @Override
    public Experiment createExperiment() {
      return new StepExperiment();
    }
  },
  RANDOMCAL("Randomized calibration") {
    @Override
    public Experiment createExperiment() {
      return new RandomizedExperiment();
    }
  },
  SINECAL("Sine calibration") {
    @Override
    public Experiment createExperiment() {
      return new SineExperiment();
    }
  },
  AZIMUTH("Azimuth") {
    @Override
    public Experiment createExperiment() {
      return new AzimuthExperiment();
    }
  },
  ORTHOGONALITY("Orthogonality") {
    @Override
    public Experiment createExperiment() {
      return new OrthogonalExperiment();
    }
  },
  SPECTRUM("Power-spectrum") {
    @Override
    public Experiment createExperiment() {
      return new SpectrumExperiment();
    }
  },
  RESPONSE("Response") {
    @Override
    public Experiment createExperiment() {
      return new ResponseExperiment();
    }
  };
  private String name;

  ExperimentFactory(String name) {
    this.name = name;
  }

  /**
   * Creates the associated Experiment for the enum
   * @return new Experiment
   */
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
