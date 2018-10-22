package asl.sensor.gui;

import asl.sensor.ExperimentFactory;
import asl.sensor.experiment.VoltageExperiment;
import asl.sensor.input.DataStore;
import java.util.List;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.annotations.XYTitleAnnotation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.title.TextTitle;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.RectangleAnchor;

public class VoltagePanel extends ExperimentPanel {


  public VoltagePanel(ExperimentFactory experiment) {
    super(experiment);

    for (int i = 0; i < panelsNeeded(); ++i) {
      channelType[i] = "Sensor under calibration (RESP required)";

    }

    chart = ChartFactory.createXYLineChart(expType.getName(),
        "Data timing (ms)", "Digital counts", null);
  }

  @Override
  public JFreeChart[] getCharts() {
    return new JFreeChart[]{chart};
  }

  @Override
  protected void drawCharts() {
    chartPanel.setChart(chart);
  }

  @Override
  public int panelsNeeded() {
    return 1;
  }

  /**
   * Displays the statistic results when the calculate button is hit
   * in an inset box on the chart, also used as text in report generation
   */
  private void setTitle() {
    JFreeChart[] charts = getCharts();
    String[] results = expResult.getInsetStrings();
    for (int i = 0; i < charts.length; ++i) {
      XYPlot plot = charts[i].getXYPlot();
      TextTitle result = getDefaultTextTitle();
      result.setText(results[i]);
      XYTitleAnnotation title = new XYTitleAnnotation(0.98, 0.98, result,
          RectangleAnchor.TOP_RIGHT);
      plot.clearAnnotations();
      plot.addAnnotation(title);
    }

  }

  @Override
  protected void updateData(DataStore dataStore) {
    VoltageExperiment ve = (VoltageExperiment) expResult;
    ve.runExperimentOnData(dataStore);
    List<XYSeriesCollection> xyscList = ve.getData();
    chart = ChartFactory.createXYLineChart(expType.getName(),
        "Data timing (ms)", "Digital counts", xyscList.get(0));
    setTitle();

  }
}
