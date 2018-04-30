package asl.sensor.gui;

import asl.sensor.experiment.ExperimentEnum;
import asl.sensor.experiment.OrthogonalExperiment;
import asl.sensor.input.DataStore;
import java.awt.Color;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.util.Arrays;
import org.jfree.chart.annotations.XYTitleAnnotation;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.title.TextTitle;
import org.jfree.ui.RectangleAnchor;

/**
 * Panel to display results of Orthogonal Experiment.
 * This plots the difference between reference and rotated sensor signals
 * and also displays the azimuth angles in an inset box.
 * No additional interface components are created for this panel.
 *
 * @author akearns - KBRWyle
 */
public class OrthogonalPanel extends ExperimentPanel {

  private static final long serialVersionUID = -2749224338484110043L;

  private static String getInsetString(OrthogonalExperiment ort) {
    double[] fit = ort.getSolutionParams();
    double angle = ort.getFitAngle();

    return "Calculated angle between non-reference sensors:\n"
        + angle
        + "\nRough est. orientation angles for (non-ref) LH1, LH2 respectively:\n"
        + Arrays.toString(fit);
  }

  OrthogonalPanel(ExperimentEnum experiment) {
    super(experiment);

    channelType[0] = "North reference sensor";
    channelType[1] = "East reference sensor";
    channelType[2] = "Assumed-north test sensor";
    channelType[3] = "Assumed-east test sensor";

    String xAxisTitle = "Time (s)";
    String yAxisTitle = "Amplitude difference (counts)";
    xAxis = new NumberAxis(xAxisTitle);
    xAxis.setAutoRange(true);

    yAxis = new NumberAxis(yAxisTitle);
    Font bold = xAxis.getLabelFont().deriveFont(Font.BOLD);
    xAxis.setLabelFont(bold);
    yAxis.setLabelFont(bold);

    applyAxesToChart();

    this.setLayout(new GridBagLayout());
    GridBagConstraints constraints = new GridBagConstraints();
    constraints.gridx = 0;
    constraints.gridy = 0;
    constraints.weightx = 1.0;
    constraints.weighty = 1.0;
    constraints.fill = GridBagConstraints.BOTH;
    constraints.anchor = GridBagConstraints.CENTER;
    this.add(chartPanel, constraints);
    constraints.weighty = 0.0;
    constraints.fill = GridBagConstraints.NONE;
    constraints.gridy += 1;
    this.add(save, constraints);

  }

  @Override
  protected void drawCharts() {
    setChart(expResult.getData().get(0));
    XYPlot plot = (XYPlot) chart.getPlot();

    TextTitle result = new TextTitle();
    result.setText(getInsetStrings());
    result.setBackgroundPaint(Color.white);
    XYTitleAnnotation xyt = new XYTitleAnnotation(0.98, 0.98, result,
        RectangleAnchor.TOP_RIGHT);
    plot.clearAnnotations();
    plot.addAnnotation(xyt);
    chartPanel.setChart(chart);
  }

  @Override
  String getInsetStrings() {
    return getInsetString((OrthogonalExperiment) expResult);
  }


  @Override
  public int panelsNeeded() {
    return 4;
  }

  @Override
  protected void updateData(final DataStore dataStore) {
    set = true;
    expResult.runExperimentOnData(dataStore);
  }
}
