package asl.sensor.gui;

import asl.sensor.ExperimentFactory;
import asl.sensor.experiment.SineExperiment;
import asl.sensor.input.DataStore;
import java.awt.Color;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import javax.swing.JComboBox;
import javax.swing.JPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.annotations.XYTitleAnnotation;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.title.TextTitle;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.RectangleAnchor;

public class SinePanel extends ExperimentPanel {

  private static final long serialVersionUID = 2453757553804095685L;
  private final NumberAxis calAxis;
  private final NumberAxis outAxis;
  private final JComboBox<String> plotSelection;
  private JFreeChart sinesChart, linearChart;

  public SinePanel(ExperimentFactory experiment) {
    super(experiment);
    channelType[0] = "Calibration input";
    channelType[1] = "Calibration output (no RESP needed)";

    xAxis = new NumberAxis("Time (s)");
    xAxis.setAutoRange(true);
    yAxis = new NumberAxis("Normalized sine wave signal");
    yAxis.setAutoRange(true);
    Font bold = xAxis.getLabelFont().deriveFont(Font.BOLD);
    xAxis.setLabelFont(bold);
    yAxis.setLabelFont(bold);
    applyAxesToChart();

    calAxis = new NumberAxis("Calibration input sample");
    outAxis = new NumberAxis("Output sample");
    calAxis.setAutoRange(true);
    outAxis.setAutoRange(true);
    calAxis.setLabelFont(bold);
    outAxis.setLabelFont(bold);

    plotSelection = new JComboBox<>();
    plotSelection.addItem("Sine plot");
    plotSelection.addItem("Linearity plot");
    plotSelection.addActionListener(this);
    plotSelection.setEnabled(false);

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
    constraints.fill = GridBagConstraints.NONE;
    constraints.gridy += 1;
    constraints.gridx = 0;
    JPanel space = new JPanel();
    space.setPreferredSize(plotSelection.getPreferredSize());
    this.add(space, constraints);

    constraints.gridx += 1;
    constraints.weightx = 1.0;
    constraints.fill = GridBagConstraints.NONE;
    constraints.anchor = GridBagConstraints.CENTER;
    this.add(save, constraints);

    // plot selection combo box
    constraints.fill = GridBagConstraints.HORIZONTAL;
    constraints.gridx += 1;
    constraints.weightx = 0;
    constraints.anchor = GridBagConstraints.WEST;
    this.add(plotSelection, constraints);
  }

  @Override
  public void actionPerformed(ActionEvent event) {
    super.actionPerformed(event);

    if (event.getSource() == plotSelection) {

      if (!set) {
        XYPlot plot = chart.getXYPlot();
        String label = getXAxis().getLabel();
        plot.getDomainAxis().setLabel(label);
        label = getYAxis().getLabel();
        plot.getRangeAxis().setLabel(label);
        return;
      }

      JFreeChart[] charts =
          new JFreeChart[]{sinesChart, linearChart};
      chart = charts[plotSelection.getSelectedIndex()];
      chartPanel.setChart(chart);
    }
  }

  @Override
  protected void drawCharts() {
    plotSelection.setSelectedIndex(0);
    plotSelection.setEnabled(true);
    sinesChart = buildChart(expResult.getData().get(0));
    XYPlot plot = (XYPlot) sinesChart.getPlot();

    TextTitle result = new TextTitle();
    result.setText(expResult.getInsetStrings()[0]);
    result.setBackgroundPaint(Color.white);
    XYTitleAnnotation title = new XYTitleAnnotation(0.98, 0.98, result,
        RectangleAnchor.TOP_RIGHT);
    plot.clearAnnotations();
    plot.addAnnotation(title);

    linearChart = buildChart(expResult.getData().get(1), calAxis, outAxis);
    appendChartTitle(linearChart, " (Linearity plot)");
    chart = sinesChart;
    chartPanel.setChart(chart);
    chartPanel.setMouseZoomable(true);
  }

  @Override
  public int panelsNeeded() {
    return expResult.blocksNeeded();
  }

  @Override
  protected void updateData(DataStore dataStore) {
    expResult.runExperimentOnData(dataStore);

    XYSeriesCollection timeseries = expResult.getData().get(0);
    for (int i = 0; i < timeseries.getSeriesCount(); ++i) {
      Color toColor = COLORS[i % COLORS.length];
      String curve = (String) timeseries.getSeriesKey(i);
      seriesColorMap.put(curve, toColor);
    }
    set = true;
  }

  @Override
  public JFreeChart[] getCharts() {
    return new JFreeChart[]{sinesChart, linearChart};
  }

  /**
   * Get the index of the data holding the sensor output.
   * Note that the input data list is listed as CAL, OUT, RESP, so the
   * relevant index is the second one
   */
  @Override
  protected int getIndexOfMainData() {
    return 1;
  }
}
