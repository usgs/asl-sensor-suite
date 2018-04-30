package asl.sensor.gui;

import asl.sensor.experiment.AzimuthExperiment;
import asl.sensor.experiment.ExperimentEnum;
import asl.sensor.input.DataStore;
import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.util.List;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JSpinner;
import javax.swing.SpinnerModel;
import javax.swing.SpinnerNumberModel;
import javax.swing.SwingConstants;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.annotations.XYTitleAnnotation;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.Marker;
import org.jfree.chart.plot.PolarPlot;
import org.jfree.chart.plot.ValueMarker;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.DefaultXYItemRenderer;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.chart.title.TextTitle;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.Layer;
import org.jfree.ui.RectangleAnchor;

/**
 * Wrapper class to display result from Azimuth. Overrides some parent
 * functions because the main plot uses polar orientation rather than typical
 * x-y plotting.
 *
 * @author akearns - KBRWyle
 */
public class AzimuthPanel extends ExperimentPanel {

  private static final long serialVersionUID = 4088024342809622854L;

  private final JSpinner offsetSpinner; // select how far from north to set reference data
  private JFreeChart angleChart, estimationChart; // plot angle, plot windowed estimation angle and correlation
  // note that some overrides are necessary because angle chart is a polar plot, not xy plot
  // so things like progress updates are called in a different manner

  private final JComboBox<String> chartSelector;

  AzimuthPanel(ExperimentEnum experiment) {
    super(experiment);

    SpinnerModel spinModel = new SpinnerNumberModel(0, -360, 360, 0.1);
    offsetSpinner = new JSpinner(spinModel);

    JLabel offsetSpinnerLabel = new JLabel("Offset angle (deg.):");
    offsetSpinnerLabel.setLabelFor(offsetSpinner);
    offsetSpinnerLabel.setHorizontalTextPosition(SwingConstants.RIGHT);
    offsetSpinnerLabel.setHorizontalAlignment(SwingConstants.RIGHT);
    JPanel labelPanel = new JPanel();
    labelPanel.add(offsetSpinnerLabel);

    chartSelector = new JComboBox<>();
    chartSelector.addItem("Azimuth angle");
    chartSelector.addItem("Estimation");
    chartSelector.setSelectedItem(0);
    chartSelector.addActionListener(this);

    plotTheseInBold = new String[]{}; // shouldn't be used anyway

    channelType[0] = "North test sensor";
    channelType[1] = "East test sensor";
    channelType[2] = "Reference sensor (use offset to specify degrees from north)";

    // don't bother instantiating axes, we need to build a custom polar plot
    // and so will just use the ChartFactory methods to do our building anyway

    angleChart = ChartFactory.createPolarChart(expType.getName(),
        null, false, false, false);
    chart = angleChart;
    chartPanel.setChart(chart);

    estimationChart =
        ChartFactory.createXYLineChart(expType.getName() + " Windowing",
            "Window start", "Correlation of aligned data, Angle of rotation", null);

    this.setLayout(new GridBagLayout());

    GridBagConstraints constraints = new GridBagConstraints();
    constraints.gridx = 0;
    constraints.gridy = 0;
    constraints.weightx = 1;
    constraints.weighty = 0;
    constraints.fill = GridBagConstraints.NONE;
    constraints.anchor = GridBagConstraints.EAST;

    constraints.anchor = GridBagConstraints.CENTER;
    constraints.gridx = 0;
    constraints.gridy = 0;
    constraints.gridwidth = 3;
    constraints.weightx = 1.0;
    constraints.weighty = 1.0;
    constraints.fill = GridBagConstraints.BOTH;
    this.add(chartPanel, constraints);

    constraints.weighty = 0.0;
    constraints.gridy += 1;
    constraints.gridwidth = 1;
    constraints.fill = GridBagConstraints.NONE;
    constraints.anchor = GridBagConstraints.EAST;
    this.add(offsetSpinnerLabel, constraints);

    constraints.gridx += 1;
    constraints.anchor = GridBagConstraints.WEST;
    this.add(offsetSpinner, constraints);

    constraints.gridx += 1;
    constraints.fill = GridBagConstraints.HORIZONTAL;
    constraints.anchor = GridBagConstraints.CENTER;
    this.add(chartSelector, constraints);

    constraints.gridx = 0;
    constraints.gridy += 1;
    constraints.gridwidth = 3;
    constraints.fill = GridBagConstraints.NONE;
    constraints.anchor = GridBagConstraints.CENTER;
    this.add(save, constraints);
  }

  @Override
  public void actionPerformed(ActionEvent event) {

    if (event.getSource() == chartSelector) {
      JFreeChart[] charts = getCharts();
      chart = charts[chartSelector.getSelectedIndex()];
      chartPanel.setChart(chart);
      return;
    }

    super.actionPerformed(event);
  }

  @Override
  protected void clearChartAndSetProgressData() {
    chartSelector.setSelectedIndex(0);
    angleChart = ChartFactory.createPolarChart(expType.getName(),
        null, false, false, false);
    chart = angleChart;
    chartPanel.setChart(chart);
    displayInfoMessage("Running calculation...");
  }

  @Override
  public void displayErrorMessage(String errorMsg) {

    if (chartSelector.getSelectedIndex() == 0) {
      PolarPlot plot = (PolarPlot) angleChart.getPlot();
      plot.clearCornerTextItems();
      plot.addCornerTextItem(errorMsg);
    } else {
      super.displayInfoMessage(errorMsg);
    }

  }

  @Override
  public void displayInfoMessage(String infoMsg) {

    if (chartSelector.getSelectedIndex() == 0) {
      PolarPlot plot = (PolarPlot) angleChart.getPlot();
      plot.clearCornerTextItems();
      plot.addCornerTextItem(infoMsg);
    } else {
      super.displayInfoMessage(infoMsg);
    }

  }

  @Override
  protected void drawCharts() {
    chartSelector.setEnabled(true);
    chartSelector.setSelectedIndex(0);
    chart = angleChart;
    chartPanel.setChart(chart);
  }

  @Override
  public String[] getAdditionalReportPages() {
    AzimuthExperiment azimuthExperiment = (AzimuthExperiment) expResult;
    double[] correlations = azimuthExperiment.getCorrelations();
    StringBuilder sb = new StringBuilder("Best-fit correlation value per-window:\n");
    for (double correlation : correlations) {
      sb.append(AzimuthPanel.DECIMAL_FORMAT.get().format(correlation)).append("  ");
    }
    sb.append("\n");

    return new String[]{sb.toString()};
  }

  @Override
  public JFreeChart[] getCharts() {
    return new JFreeChart[]{angleChart, estimationChart};
  }

  @Override
  public String getInsetStrings() {
    AzimuthExperiment experiment = (AzimuthExperiment) expResult;
    double value = experiment.getOffset();
    double angle = experiment.getFitAngle();
    StringBuilder angleStr = new StringBuilder();
    angleStr.append("FIT ANGLE: ").append(DECIMAL_FORMAT.get().format(angle));
    double result = ((value + angle) % 360 + 360) % 360;

    angleStr.append(" + ").append(DECIMAL_FORMAT.get().format(value)).append(" = ");
    angleStr.append(DECIMAL_FORMAT.get().format(result)).append(" (+/- ");
    angleStr.append(DECIMAL_FORMAT.get().format(experiment.getUncertainty())).append(")");
    if (!experiment.hadEnoughPoints()) {
      angleStr.append(" | WARNING: SMALL RANGE");
    }
    return angleStr.toString();
  }

  @Override
  public int panelsNeeded() {
    return 3;
  }

  @Override
  protected void updateData(DataStore dataStore) {

    set = true;

    double value = (double) offsetSpinner.getValue();

    if (value < 0) {
      value += 360;
    }

    AzimuthExperiment experiment = (AzimuthExperiment) expResult;
    experiment.setOffset(value);

    XYPlot estimationPlot;

    expResult.runExperimentOnData(dataStore);
    List<XYSeriesCollection> allData = experiment.getData();
    XYSeriesCollection polars = allData.get(0);

    angleChart = ChartFactory.createPolarChart(expType.getName(),
        polars, true, true, false);

    String angleStr = getInsetStrings();

    PolarPlot plot = (PolarPlot) angleChart.getPlot();
    plot.clearCornerTextItems();
    plot.addCornerTextItem(angleStr);

    XYSeriesCollection angleEstimation = allData.get(1);
    XYSeriesCollection coherenceEstimation = allData.get(2);
    String titleEst = expType.getName() + " Accuracy Estimation";
    estimationChart = ChartFactory.createXYLineChart(titleEst,
        "xAxis", "yAxis", angleEstimation);
    estimationPlot = estimationChart.getXYPlot();
    estimationPlot.setDataset(0, angleEstimation);
    estimationPlot.setDataset(1, coherenceEstimation);
    estimationPlot.setRenderer(0, new DefaultXYItemRenderer());

    // set color of second dataset to be blue
    XYItemRenderer renderer = new DefaultXYItemRenderer();
    renderer.setSeriesPaint(0, Color.BLUE);
    estimationPlot.setRenderer(1, renderer);

    NumberAxis angleEstimationAxis = new NumberAxis("Angle est. (deg)");
    estimationPlot.setRangeAxis(0, angleEstimationAxis);
    NumberAxis correlationAxis = new NumberAxis("Correlation est. of best fit angle");
    estimationPlot.setRangeAxis(1, correlationAxis);
    NumberAxis xAxis = new NumberAxis("Time from data start of (2000s) window (s)");
    xAxis.setAutoRangeIncludesZero(false);
    estimationPlot.setDomainAxis(xAxis);
    estimationPlot.mapDatasetToRangeAxis(0, 0);
    estimationPlot.mapDatasetToRangeAxis(1, 1);

    if (!experiment.hadEnoughPoints()) {
      estimationPlot = estimationChart.getXYPlot();
      TextTitle result = new TextTitle();
      result.setText("WARNING: NOT ENOUGH DATA FOR WINDOWED COHERENCE ESTIMATION");
      result.setBackgroundPaint(Color.red);
      result.setPaint(Color.white);
      XYTitleAnnotation xyt = new XYTitleAnnotation(0.5, 0.5, result,
          RectangleAnchor.CENTER);
      estimationPlot.clearAnnotations();
      estimationPlot.addAnnotation(xyt);
    } else {
      double cutOff = experiment.getMinCorr();
      Marker highWater = new ValueMarker(cutOff);
      highWater.setStroke(new BasicStroke((float) 1.5));
      highWater.setPaint(Color.BLACK);
      estimationPlot.addRangeMarker(1, highWater, Layer.BACKGROUND);
    }

    chartSelector.setSelectedIndex(0);
  }
}
