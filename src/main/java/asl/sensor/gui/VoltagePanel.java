package asl.sensor.gui;

import asl.sensor.ExperimentFactory;
import asl.sensor.experiment.VoltageExperiment;
import asl.sensor.input.DataStore;
import java.awt.Color;
import java.awt.Font;
import java.util.List;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.annotations.XYTitleAnnotation;
import org.jfree.chart.axis.DateAxis;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.title.TextTitle;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.RectangleAnchor;

public class VoltagePanel extends ExperimentPanel {

 private int plotCount;

  public VoltagePanel(ExperimentFactory experiment) {
    super(experiment);

    plotCount = 0;

    for (int i = 0; i < panelsNeeded(); ++i) {
      channelType[i] = "Sensor under calibration (RESP required)";

    }

    String xAxisTitle = "Time (UTC)";
    String yAxisTitle = "Digital counts";
    xAxis = new DateAxis(xAxisTitle);
    ((DateAxis) xAxis).setDateFormatOverride(ExperimentPanel.DATE_TIME_FORMAT.get());
    Font bold = xAxis.getLabelFont();
    bold = bold.deriveFont(Font.BOLD, bold.getSize() + 2);
    xAxis.setLabelFont(bold);
    yAxis = new NumberAxis(yAxisTitle);

    applyAxesToChart();
  }

  @Override
  protected void drawCharts() {
    JFreeChart chart = ChartFactory.createScatterPlot(getName(),
        "", "", expResult.getData().get(0));
    chart.getXYPlot().setRangeAxis(yAxis);
    chart.getXYPlot().setDomainAxis(xAxis);
    this.chart = chart;
    chartPanel.setChart(chart);
    setTitle();
    chartPanel.setMouseZoomable(true);
  }

  @Override
  public int panelsNeeded() {
    return 3;
  }

  /**
   * Displays the statistic results when the calculate button is hit
   * in an inset box on the chart, also used as text in report generation
   */
  private void setTitle() {
    // chart = chartPanel.getChart();
    String results = expResult.getInsetStrings()[0];
    XYPlot plot = chart.getXYPlot();
    TextTitle result = getDefaultTextTitle();
    result.setText(results);
    XYTitleAnnotation title = new XYTitleAnnotation(0.5, 0.5, result,
        RectangleAnchor.CENTER);
    plot.clearAnnotations();
    plot.addAnnotation(title);

  }

  @Override
  protected void updateData(DataStore dataStore) {
    expResult.runExperimentOnData(dataStore);

    plotCount = 0;
    for (int i = 0; i < panelsNeeded(); ++i) {
      if (dataStore.bothComponentsSet(i)) {
        ++plotCount;
      }
    }

    expResult.runExperimentOnData(dataStore);

    XYSeriesCollection timeseries = expResult.getData().get(0);

    for (int i = 0; i < plotCount; ++i) {
      String name = (String) timeseries.getSeriesKey(i);
      if (null == name) {
        continue;
      }
      Color plotColor = COLORS[i % 3];
      seriesColorMap.put(name, plotColor);
    }

    set = true;
    /*
    List<XYSeriesCollection> xyscList = expResult.getData();
    chart = ChartFactory.createXYLineChart(expType.getName(),
        "Data timing (ms)", "Digital counts", xyscList.get(0));
    setTitle();
    */
  }
}
