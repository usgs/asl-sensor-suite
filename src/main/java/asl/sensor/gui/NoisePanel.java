package asl.sensor.gui;

import asl.sensor.ExperimentFactory;
import asl.sensor.experiment.NoiseExperiment;
import asl.sensor.experiment.SpectralAnalysisExperiment;
import asl.sensor.input.DataStore;
import java.awt.Color;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import javax.swing.JCheckBox;
import org.jfree.chart.annotations.XYTitleAnnotation;
import org.jfree.chart.axis.LogarithmicAxis;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.title.TextTitle;
import org.jfree.chart.ui.RectangleAnchor;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

/**
 * Panel for displaying the results of the self-noise experiment (3-input).
 * In addition to general requirements of output panels, also includes
 * a checkbox to choose between frequency and interval x-axis and
 * the variant axes for when that box is checked.
 *
 * @author akearns - KBRWyle
 */
public class NoisePanel extends SpectralAnalysisPanel {

  // TODO: add buttons for specifying and clearing a station noise model
  //   and logic to combine the noise model with other plotted data
  //   (will likely need to add this to 3D version of panel as well)

  private static final long serialVersionUID = 9018553361096758354L;

  final JCheckBox freqSpaceBox;

  final int NOISE_PLOT_COUNT = 6;

  private final NumberAxis freqAxis;

  /**
   * Constructs a new panel and lays out all the components in it
   */
  public NoisePanel(ExperimentFactory experiment) {

    // create chart, chartPanel, save button & file chooser,
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

    freqSpaceBox = new JCheckBox("Use Hz units (requires regen)");
    freqSpaceBox.setSelected(false);

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
    XYSeriesCollection data = expResult.getData().get(0);
    {
      SpectralAnalysisExperiment spectExp = (SpectralAnalysisExperiment) expResult;
      if (spectExp.noiseModelLoaded()) {
        XYSeries series = spectExp.getPlottableNoiseModelData(freqSpaceBox.isSelected());
        seriesDashedSet.add((String) series.getKey());
        seriesColorMap.put((String) series.getKey(), Color.BLACK);
        data.addSeries(spectExp.getPlottableNoiseModelData(freqSpaceBox.isSelected()));
      }
    }
    setChart(data);
    setTitle(chart.getXYPlot(), expResult.getInsetStrings()[0]);
    chartPanel.setChart(chart);
    chartPanel.setMouseZoomable(true);
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
  public int panelsNeeded() {
    return 3;
  }

  /**
   * Initially called function to calculate self-noise when data is passed in
   */
  @Override
  protected void updateData(final DataStore dataStore) {
    set = true;

    ((NoiseExperiment) expResult).setFreqSpace(freqSpaceBox.isSelected());
    expResult.runExperimentOnData(dataStore);

    XYSeriesCollection timeseries = expResult.getData().get(0);

    for (int i = 0; i < NOISE_PLOT_COUNT; ++i) {
      String name = (String) timeseries.getSeriesKey(i);
      Color plotColor = getColor(i);
      seriesColorMap.put(name, plotColor);
      if (i >= 3) {
        seriesDashedSet.add(name);
      }
    }
  }

  void setTitle(XYPlot plot, String insetString) {
    TextTitle angle = getDefaultTextTitle();
    angle.setText(insetString);
    XYTitleAnnotation title = new XYTitleAnnotation(0.98, 0.98, angle,
        RectangleAnchor.TOP_RIGHT);
    plot.clearAnnotations();
    plot.addAnnotation(title);
  }

}
