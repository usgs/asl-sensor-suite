package asl.sensor;

import static org.junit.Assert.*;
import static org.hamcrest.CoreMatchers.instanceOf;


import asl.sensor.experiment.AzimuthExperiment;
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
import org.junit.Test;

/**
 * These tests verify that each enum matches its expected classes for createExperiment() and createPanel().
 */
public class ExperimentFactoryTest {

  @Test
  public void noise_createExperiment() {
    assertThat(ExperimentFactory.NOISE.createExperiment(), instanceOf(NoiseExperiment.class));
  }

  @Test
  public void noise_createPanel() {
    assertThat(ExperimentFactory.NOISE.createPanel(), instanceOf(NoisePanel.class));
  }

  @Test
  public void noise9_createExperiment() {
    assertThat(ExperimentFactory.NOISE9.createExperiment(), instanceOf(NoiseNineExperiment.class));
  }

  @Test
  public void noise9_createPanel() {
    assertThat(ExperimentFactory.NOISE9.createPanel(), instanceOf(NoiseNinePanel.class));
  }

  @Test
  public void gain_createExperiment() {
    assertThat(ExperimentFactory.GAIN.createExperiment(), instanceOf(GainExperiment.class));
  }

  @Test
  public void gain_createPanel() {
    assertThat(ExperimentFactory.GAIN.createPanel(), instanceOf(GainPanel.class));
  }

  @Test
  public void gain6_createExperiment() {
    assertThat(ExperimentFactory.GAIN6.createExperiment(), instanceOf(GainSixExperiment.class));
  }

  @Test
  public void gain6_createPanel() {
    assertThat(ExperimentFactory.GAIN6.createPanel(), instanceOf(GainSixPanel.class));
  }

  @Test
  public void stepCal_createExperiment() {
    assertThat(ExperimentFactory.STEPCAL.createExperiment(), instanceOf(StepExperiment.class));
  }

  @Test
  public void stepCal_createPanel() {
    assertThat(ExperimentFactory.STEPCAL.createPanel(), instanceOf(StepPanel.class));
  }

  @Test
  public void randomCal_createExperiment() {
    assertThat(ExperimentFactory.RANDOMCAL.createExperiment(), instanceOf(RandomizedExperiment.class));
  }

  @Test
  public void randomCal_createPanel() {
    assertThat(ExperimentFactory.RANDOMCAL.createPanel(), instanceOf(RandomizedPanel.class));
  }

  @Test
  public void sineCal_createExperiment() {
    assertThat(ExperimentFactory.SINECAL.createExperiment(), instanceOf(SineExperiment.class));
  }

  @Test
  public void sineCal_createPanel() {
    assertThat(ExperimentFactory.SINECAL.createPanel(), instanceOf(SinePanel.class));
  }

  @Test
  public void azimuth_createExperiment() {
    assertThat(ExperimentFactory.AZIMUTH.createExperiment(), instanceOf(AzimuthExperiment.class));
  }

  @Test
  public void azimuth_createPanel() {
    assertThat(ExperimentFactory.AZIMUTH.createPanel(), instanceOf(AzimuthPanel.class));
  }

  @Test
  public void orthogonality_createExperiment() {
    assertThat(ExperimentFactory.ORTHOGONALITY.createExperiment(), instanceOf(OrthogonalExperiment.class));
  }

  @Test
  public void orthogonality_createPanel() {
    assertThat(ExperimentFactory.ORTHOGONALITY.createPanel(), instanceOf(OrthogonalPanel.class));
  }

  @Test
  public void spectrum_createExperiment() {
    assertThat(ExperimentFactory.SPECTRUM.createExperiment(), instanceOf(SpectrumExperiment.class));
  }

  @Test
  public void spectrum_createPanel() {
    assertThat(ExperimentFactory.SPECTRUM.createPanel(), instanceOf(SpectrumPanel.class));
  }

  @Test
  public void response_createExperiment() {
    assertThat(ExperimentFactory.RESPONSE.createExperiment(), instanceOf(ResponseExperiment.class));
  }

  @Test
  public void response_createPanel() {
    assertThat(ExperimentFactory.RESPONSE.createPanel(), instanceOf(ResponsePanel.class));
  }
}