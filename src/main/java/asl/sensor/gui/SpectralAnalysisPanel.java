package asl.sensor.gui;

import asl.sensor.ExperimentFactory;
import asl.sensor.experiment.SpectralAnalysisExperiment;
import asl.sensor.input.Configuration;
import asl.sensor.input.NoiseModel.NoiseModelFormatException;
import java.awt.event.ActionEvent;
import java.io.File;
import java.io.IOException;
import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JOptionPane;

public abstract class SpectralAnalysisPanel extends ExperimentPanel {

  final JButton noiseModelButton;

  /**
   * Construct a new panel, using a backend defined by the passed-in enum
   *
   * @param experiment Experiment enum with corresponding backend for factory instantiation
   */
  SpectralAnalysisPanel(ExperimentFactory experiment) {
    super(experiment);

    noiseModelButton = new JButton("Set noise model");
    noiseModelButton.setEnabled(true);
    noiseModelButton.addActionListener(this);
  }

  @Override
  public void actionPerformed(ActionEvent event) {
    super.actionPerformed(event);

    if (event.getSource() == noiseModelButton) {
      SpectralAnalysisExperiment exp = (SpectralAnalysisExperiment) expResult;
      if (exp.noiseModelLoaded()) {
        Object[] options = {"Cancel", "Remove", "Load new"};
        String message = "Current model: " + exp.getNoiseModel().getName()
            + "\nLoad or remove current model?";
        int selected = JOptionPane.showOptionDialog(this,
            message, "Model load check",
            JOptionPane.DEFAULT_OPTION, JOptionPane.QUESTION_MESSAGE, null, options, options[1]);
        if (selected == 0) {
          return; // cancel option selected, do nothing
        } else if (selected == 1) {
          exp.clearNoiseModel();
          if (exp.getData() != null) {
            drawCharts();
          }
          return;
        }
      }
      // if we're here, load new was selected or no noise model was loaded yet
      JFileChooser chooser = new JFileChooser(
          Configuration.getInstance().getDefaultNoiseModelFolder());
      int returnVal = chooser.showOpenDialog(this);
      if (returnVal == JFileChooser.APPROVE_OPTION) {
        File file = chooser.getSelectedFile().getAbsoluteFile();
        try {
          exp.loadNoiseModel(file);
          if (exp.getData() != null) {
            drawCharts();
          }
        } catch (IOException | NoiseModelFormatException e) {
          String errorMessage = "Could not load " + file.getAbsolutePath() + "\n"
              + "Check that file is valid noise model format.";
          JOptionPane.showMessageDialog(this, errorMessage, "NOISE MODEL ERROR",
              JOptionPane.ERROR_MESSAGE);
          e.printStackTrace();
        }
      }
    }
  }
}
