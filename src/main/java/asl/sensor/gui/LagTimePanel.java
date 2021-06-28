package asl.sensor.gui;

import asl.sensor.ExperimentFactory;
import asl.sensor.input.DataStore;
import java.awt.Color;
import java.awt.Font;
import org.jfree.chart.annotations.XYTitleAnnotation;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.title.TextTitle;
import org.jfree.chart.ui.RectangleAnchor;
import org.jfree.data.xy.XYSeriesCollection;

public class LagTimePanel extends ExperimentPanel {

  /**
   * Construct a new panel, using a backend defined by the passed-in enum
   *
   * @param experiment Experiment enum with corresponding backend for factory instantiation
   */
  public LagTimePanel(ExperimentFactory experiment) {
    super(experiment);

    channelType[0] = "Sensor under test (RESP required)";
    channelType[1] = "Sensor for lag reference (RESP required)";

    xAxis = new NumberAxis("Lag time (ms)");
    yAxis = new NumberAxis("Correlation value between data");
    xAxis.setAutoRange(true);
    yAxis.setAutoRange(true);
    ((NumberAxis) yAxis).setAutoRangeIncludesZero(false);
    Font bold = xAxis.getLabelFont();
    bold = bold.deriveFont(Font.BOLD, bold.getSize() + 2);
    xAxis.setLabelFont(bold);
    yAxis.setLabelFont(bold);
    applyAxesToChart();
  }

  @Override
  protected void drawCharts() {
    XYSeriesCollection timeSeries = expResult.getData().get(0);
    setChart(timeSeries);
    chartPanel.setChart(chart);
    // setTitle();
    chartPanel.setMouseZoomable(true);
  }

  @Override
  public int panelsNeeded() {
    return expResult.blocksNeeded();
  }

  @Override
  protected void updateData(DataStore dataStore) {
    set = true;
    expResult.runExperimentOnData(dataStore);

    XYSeriesCollection timeseries = expResult.getData().get(0);
    seriesColorMap.put((String) timeseries.getSeriesKey(0), Color.BLACK);
  }

  /**
   * Displays the statistic results when the calculate button is hit
   * in an inset box on the chart, also used as text in report generation
   */
  private void setTitle() {
    XYPlot plot = (XYPlot) chartPanel.getChart().getPlot();
    TextTitle result = getDefaultTextTitle();
    result.setText(expResult.getInsetStrings()[0]);
    XYTitleAnnotation xyt = new XYTitleAnnotation(0.98, 0.98, result,
        RectangleAnchor.TOP_RIGHT);
    plot.clearAnnotations();
    plot.addAnnotation(xyt);
  }

}
