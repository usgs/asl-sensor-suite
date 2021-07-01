package asl.sensor.gui;

import asl.sensor.ExperimentFactory;
import asl.sensor.input.DataStore;
import java.awt.Color;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import javax.swing.JComboBox;
import javax.swing.JPanel;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.annotations.XYTitleAnnotation;
import org.jfree.chart.axis.DateAxis;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.title.TextTitle;
import org.jfree.chart.ui.RectangleAnchor;
import org.jfree.data.xy.XYSeriesCollection;

public class LagTimePanel extends ExperimentPanel {

  private final JComboBox<String> plotSelection; // which chart to display in window?
  private JFreeChart correlationChart; // plot of correlation function
  private JFreeChart correctedChart; // plot of time-corrected data
  private ValueAxis dateAxis, signalAxis;

  private static final String[] TITLES =
      {" (Correlation)", " (Time Correction)"};

  /**
   * Construct a new panel, using a backend defined by the passed-in enum
   *
   * @param experiment Experiment enum with corresponding backend for factory instantiation
   */
  public LagTimePanel(ExperimentFactory experiment) {
    super(experiment);

    channelType[0] = "Sensor under test (RESP required)";
    channelType[1] = "Sensor for lag reference (RESP required)";

    String xAxisTitle = "Time (UTC), corrected";
    dateAxis = new DateAxis(xAxisTitle);
    ((DateAxis) dateAxis).setDateFormatOverride(ExperimentPanel.DATE_TIME_FORMAT.get());
    Font bold = dateAxis.getLabelFont();
    bold = bold.deriveFont(Font.BOLD, bold.getSize() + 2);
    dateAxis.setLabelFont(bold);
    xAxis = new NumberAxis("Time delay (ms)");
    yAxis = new NumberAxis("Correlation value between data");
    xAxis.setAutoRange(true);
    yAxis.setAutoRange(true);
    ((NumberAxis) yAxis).setAutoRangeIncludesZero(false);
    xAxis.setLabelFont(bold);
    yAxis.setLabelFont(bold);
    signalAxis = new NumberAxis("Counts");
    signalAxis.setAutoRange(true);
    ((NumberAxis) signalAxis).setAutoRangeIncludesZero(false);
    signalAxis.setLabelFont(bold);

    plotSelection = new JComboBox<>();
    plotSelection.addItem("Correlation function");
    plotSelection.addItem("Delay-adjusted data");
    plotSelection.addActionListener(this);

    correlationChart = ChartFactory.createXYLineChart(expType.getName() + TITLES[0],
        "", "", null);
    correctedChart = ChartFactory.createXYLineChart(expType.getName() + TITLES[1],
        "", "", null);
    chart = getCharts()[plotSelection.getSelectedIndex()];
    applyAxesToChart();
    chartPanel.setChart(chart);

    removeAll();
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

    constraints.gridx = 0;
    constraints.gridy += 1;
    constraints.weighty = 0;
    constraints.fill = GridBagConstraints.NONE;
    constraints.anchor = GridBagConstraints.WEST;
    JPanel spacer = new JPanel();
    spacer.setMaximumSize(plotSelection.getMaximumSize());
    spacer.setMinimumSize(plotSelection.getMinimumSize());
    spacer.setPreferredSize(plotSelection.getPreferredSize());
    spacer.setSize(plotSelection.getSize());
    this.add(spacer, constraints);
    constraints.weightx = 0;
    constraints.gridx += 1;
    constraints.anchor = GridBagConstraints.CENTER;
    this.add(save, constraints);

    constraints.gridx += 1;
    constraints.weightx = 0;
    constraints.anchor = GridBagConstraints.EAST;
    add(plotSelection, constraints);
  }

  @Override
  public ValueAxis getXAxis() {
    if (null == plotSelection) {
      return xAxis;
    }

    ValueAxis[] array = new ValueAxis[]{xAxis, dateAxis};
    return array[plotSelection.getSelectedIndex()];
  }

  @Override
  public ValueAxis getYAxis() {
    if (null == plotSelection) {
      return yAxis;
    }

    ValueAxis[] array = new ValueAxis[]{yAxis, signalAxis};
    return array[plotSelection.getSelectedIndex()];
  }

  @Override
  protected void drawCharts() {
    {
      XYSeriesCollection timeSeries = expResult.getData().get(0);
      seriesColorMap.put((String) timeSeries.getSeriesKey(0), Color.BLACK);
      correlationChart = buildChart(timeSeries);
      correlationChart.setTitle(expType.getName() + TITLES[0]);
    }
    {
      XYSeriesCollection timeSeries = expResult.getData().get(1);
      for (int i = 0; i < timeSeries.getSeriesCount(); ++i) {
        seriesColorMap.put((String) timeSeries.getSeriesKey(i), getColor(i));
      }
      correctedChart = buildChart(timeSeries);
      correctedChart.setTitle(expType.getName() + TITLES[1]);
    }
    setTitle();
    plotSelection.setSelectedIndex(0);
    chartPanel.setMouseZoomable(true);
  }

  @Override
  public JFreeChart[] getCharts() {
    return new JFreeChart[]{correlationChart, correctedChart};
  }

  @Override
  public int panelsNeeded() {
    return expResult.blocksNeeded();
  }

  @Override
  protected void updateData(DataStore dataStore) {
    set = true;
    expResult.runExperimentOnData(dataStore);
  }

  @Override
  public void actionPerformed(ActionEvent event) {
    super.actionPerformed(event); // saving?

    if (event.getSource() == plotSelection) {
      int index = plotSelection.getSelectedIndex();
      JFreeChart[] charts = getCharts();
      chart = charts[index];
      applyAxesToChart();
      chartPanel.setChart(chart);
    }
  }

  /**
   * Displays the statistic results when the calculate button is hit
   * in an inset box on the chart, also used as text in report generation
   */
  private void setTitle() {
    JFreeChart[] charts = getCharts();
    String[] results = expResult.getInsetStrings();
    for (int i = 0; i < results.length; ++i) {
      XYPlot plot = charts[i].getXYPlot();
      TextTitle result = getDefaultTextTitle();
      result.setText(results[i]);
      XYTitleAnnotation title = new XYTitleAnnotation(0.98, 0.98, result,
          RectangleAnchor.TOP_RIGHT);
      plot.clearAnnotations();
      plot.addAnnotation(title);
    }
  }

}
