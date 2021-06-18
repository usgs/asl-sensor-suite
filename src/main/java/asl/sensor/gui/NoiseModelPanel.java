package asl.sensor.gui;

import asl.sensor.ExperimentFactory;
import asl.sensor.experiment.SpectralAnalysisExperiment;
import asl.sensor.input.DataStore;
import java.awt.Color;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import javax.swing.JCheckBox;
import org.jfree.chart.axis.LogarithmicAxis;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

public class NoiseModelPanel extends SpectralAnalysisPanel {

  final JCheckBox freqSpaceBox;

  private final NumberAxis freqAxis;

  /**
   * Construct a new panel, using a backend defined by the passed-in enum
   *
   * @param experiment Experiment enum with corresponding backend for factory instantiation
   */
  public NoiseModelPanel(ExperimentFactory experiment) {
    super(experiment);

    for (int i = 0; i < 3; ++i) {
      channelType[i] = "Input data (RESP required)";
    }

    plotTheseInBold = new String[]{"NLNM", "NHNM"};

    xAxis = new LogarithmicAxis("Period (s)");
    freqAxis = new LogarithmicAxis("Frequency (Hz)");
    yAxis = new NumberAxis("Power (rel. 1 (m/s^2)^2/Hz)");
    yAxis.setAutoRange(true);
    ((NumberAxis) yAxis).setAutoRangeIncludesZero(false);
    Font bold = xAxis.getLabelFont();
    bold = bold.deriveFont(Font.BOLD, bold.getSize() + 2);
    xAxis.setLabelFont(bold);
    yAxis.setLabelFont(bold);
    freqAxis.setLabelFont(bold);

    freqSpaceBox = new JCheckBox("Use Hz units");
    freqSpaceBox.setSelected(false);
    freqSpaceBox.addActionListener(this);

    applyAxesToChart(); // now that we've got axes defined

    // set the GUI components
    this.setLayout(new GridBagLayout());
    GridBagConstraints constraints = new GridBagConstraints();

    constraints.fill = GridBagConstraints.BOTH;
    constraints.gridx = 0;
    constraints.gridy = 0;
    constraints.weightx = 1.0;
    constraints.weighty = 1.0;
    constraints.gridwidth = 3;
    constraints.anchor = GridBagConstraints.CENTER;
    this.add(chartPanel, constraints);

    // place the other UI elements in a single row below the chart
    constraints.gridwidth = 1;
    constraints.weighty = 0.0;
    constraints.weightx = 0.0;
    constraints.anchor = GridBagConstraints.WEST;
    constraints.fill = GridBagConstraints.HORIZONTAL;
    constraints.gridy += 1;
    constraints.gridx = 0;
    freqSpaceBox.setPreferredSize(noiseModelButton.getPreferredSize());
    freqSpaceBox.setMaximumSize(noiseModelButton.getMaximumSize());
    freqSpaceBox.setMinimumSize(noiseModelButton.getMinimumSize());
    freqSpaceBox.setSize(noiseModelButton.getSize());
    this.add(freqSpaceBox, constraints);

    constraints.gridx += 1;
    constraints.weightx = 1.0;
    constraints.fill = GridBagConstraints.NONE;
    constraints.anchor = GridBagConstraints.CENTER;
    this.add(save, constraints);

    // add an empty panel as a spacer to keep the save button in the center
    constraints.fill = GridBagConstraints.HORIZONTAL;
    constraints.gridx += 1;
    constraints.weightx = 0;
    constraints.anchor = GridBagConstraints.EAST;
    this.add(noiseModelButton, constraints);
  }

  @Override
  protected void drawCharts() {
    XYSeriesCollection data = new XYSeriesCollection();
    {
      SpectralAnalysisExperiment spectExp = (SpectralAnalysisExperiment) expResult;
      if (spectExp.noiseModelLoaded()) {
        XYSeries series = spectExp.getPlottableNoiseModelData(freqSpaceBox.isSelected());
        for (int i = 0; i < series.getItemCount(); ++i) {
          System.out.println(series.getDataItem(i));
        }
        seriesColorMap.put((String) series.getKey(), Color.BLACK);
        data.addSeries(spectExp.getPlottableNoiseModelData(freqSpaceBox.isSelected()));
      }
    }
    setChart(data);
    chartPanel.setChart(chart);
    chartPanel.setMouseZoomable(true);
  }

  @Override
  public int panelsNeeded() {
    return 0;
  }

  @Override
  protected void updateData(DataStore dataStore) {
    // this also does nothing, as no data beyond the noise model is needed
  }

  /**
   * Gets the x-axis for this panel based on whether or not the
   * selection box to plot in units of Hz is selected. If it is, this
   * plot will have frequency units of Hz in the x-axis, otherwise it will have
   * interval units of seconds in it
   */
  @Override
  public ValueAxis getXAxis() {
    // true if using Hz units
    if (freqSpaceBox.isSelected()) {
      return freqAxis;
    }
    return xAxis;

  }

  @Override
  public void actionPerformed(ActionEvent event) {
    super.actionPerformed(event);
    // override because this way we draw the charts anyway
    if (event.getSource() == noiseModelButton || event.getSource() == freqSpaceBox) {
      drawCharts();
    }
  }

}
