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
import asl.sensor.gui.AzimuthPanel;
import asl.sensor.gui.ExperimentPanel;
import asl.sensor.gui.GainPanel;
import asl.sensor.gui.GainSixPanel;
import asl.sensor.gui.NoiseNinePanel;
import asl.sensor.gui.NoisePanel;
import asl.sensor.gui.OrthogonalPanel;
import asl.sensor.gui.RandomizedPanel;
import asl.sensor.gui.ResponsePanel;
import asl.sensor.gui.SinePanel;
import asl.sensor.gui.SpectrumPanel;
import asl.sensor.gui.StepPanel;

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

    @Override
    public ExperimentPanel createPanel() {
      return new NoisePanel(this);
    }
  },
  NOISE9("Self-noise (3-component)") {
    @Override
    public Experiment createExperiment() {
      return new NoiseNineExperiment();
    }

    @Override
    public ExperimentPanel createPanel() {
      return new NoiseNinePanel(this);
    }
  },
  GAIN("Relative gain") {
    @Override
    public Experiment createExperiment() {
      return new GainExperiment();
    }

    @Override
    public ExperimentPanel createPanel() {
      return new GainPanel(this);
    }
  },
  GAIN6("Relative gain (3-component)") {
    @Override
    public Experiment createExperiment() {
      return new GainSixExperiment();
    }

    @Override
    public ExperimentPanel createPanel() {
      return new GainSixPanel(this);
    }
  },
  STEPCAL("Step calibration") {
    @Override
    public Experiment createExperiment() {
      return new StepExperiment();
    }

    @Override
    public ExperimentPanel createPanel() {
      return new StepPanel(this);
    }
  },
  RANDOMCAL("Randomized calibration") {
    @Override
    public Experiment createExperiment() {
      return new RandomizedExperiment();
    }

    @Override
    public ExperimentPanel createPanel() {
      return new RandomizedPanel(this);
    }
  },
  SINECAL("Sine calibration") {
    @Override
    public Experiment createExperiment() {
      return new SineExperiment();
    }

    @Override
    public ExperimentPanel createPanel() {
      return new SinePanel(this);
    }
  },
  AZIMUTH("Azimuth") {
    @Override
    public Experiment createExperiment() {
      return new AzimuthExperiment();
    }

    @Override
    public ExperimentPanel createPanel() {
      return new AzimuthPanel(this);
    }
  },
  ORTHOGONALITY("Orthogonality") {
    @Override
    public Experiment createExperiment() {
      return new OrthogonalExperiment();
    }

    @Override
    public ExperimentPanel createPanel() {
      return new OrthogonalPanel(this);
    }
  },
  SPECTRUM("Power-spectrum") {
    @Override
    public Experiment createExperiment() {
      return new SpectrumExperiment();
    }

    @Override
    public ExperimentPanel createPanel() {
      return new SpectrumPanel(this);
    }
  },
  RESPONSE("Response") {
    @Override
    public Experiment createExperiment() {
      return new ResponseExperiment();
    }

    @Override
    public ExperimentPanel createPanel() {
      return new ResponsePanel(this);
    }
  };
  private final String name;

  ExperimentFactory(String name) {
    this.name = name;
  }

  /**
   * Creates the associated Experiment for the enum
   *
   * @return new Experiment
   */
  public abstract Experiment createExperiment();

  /**
   * Create a new panel for the GUI.
   *
   * @return an ExperimentPanel for associated class.
   */
  public abstract ExperimentPanel createPanel();

  /**
   * Get the full name of this experiment (used for plot & tab names)
   *
   * @return Name of experiment, as String
   */
  public String getName() {
    return name;
  }

}
