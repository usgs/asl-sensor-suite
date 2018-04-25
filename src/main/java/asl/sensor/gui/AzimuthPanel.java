package asl.sensor.gui;

import asl.sensor.experiment.AzimuthExperiment;
import asl.sensor.experiment.ExperimentEnum;
import asl.sensor.input.DataStore;
import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.text.DecimalFormat;
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
  /**
   * Thread safe reference to a shared DecimalFormat object.
   */
  private static final ThreadLocal<DecimalFormat> DECIMAL_FORMAT =
      ThreadLocal.withInitial(() ->  new DecimalFormat("#.###"));
  private final JSpinner offsetSpinner; // select how far from north to set reference data
  private JFreeChart angleChart, estimChart; // plot angle, plot windowed estimation angle and correlation
  // note that some overrides are necessary because angle chart is a polar plot, not xy plot
  // so things like progress updates are called in a different manner

  private final JComboBox<String> chartSelector;

  AzimuthPanel(ExperimentEnum exp) {
    super(exp);

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
    // and so will just use the chartfactory methods to do our building anyway

    angleChart = ChartFactory.createPolarChart(expType.getName(),
        null, false, false, false);
    chart = angleChart;
    chartPanel.setChart(chart);

    estimChart =
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
  public void actionPerformed(ActionEvent e) {

    if (e.getSource() == chartSelector) {
      JFreeChart[] charts = getCharts();
      chart = charts[chartSelector.getSelectedIndex()];
      chartPanel.setChart(chart);
      return;
    }

    super.actionPerformed(e);
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
  public void displayErrorMessage(String errMsg) {

    if (chartSelector.getSelectedIndex() == 0) {
      PolarPlot plot = (PolarPlot) angleChart.getPlot();
      plot.clearCornerTextItems();
      plot.addCornerTextItem(errMsg);
    } else {
      super.displayInfoMessage(errMsg);
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
    double[] corr = azimuthExperiment.getCorrelations();
    StringBuilder sb = new StringBuilder("Best-fit correlation value per-window:\n");
    for (double aCorr : corr) {
      sb.append(AzimuthPanel.DECIMAL_FORMAT.get().format(aCorr)).append("  ");
    }
    sb.append("\n");

    return new String[]{sb.toString()};
  }

  @Override
  public JFreeChart[] getCharts() {
    return new JFreeChart[]{angleChart, /*coherenceChart,*/ estimChart};
  }

  @Override
  public String getInsetStrings() {
    AzimuthExperiment az = (AzimuthExperiment) expResult;
    double value = az.getOffset();
    double angle = az.getFitAngle();
    StringBuilder angleStr = new StringBuilder();
    angleStr.append("FIT ANGLE: ").append(DECIMAL_FORMAT.get().format(angle));
    double result = ((value + angle) % 360 + 360) % 360;

    angleStr.append(" + ").append(DECIMAL_FORMAT.get().format(value)).append(" = ")
        .append(DECIMAL_FORMAT.get().format(result));
    angleStr.append(" (+/- ").append(DECIMAL_FORMAT.get().format(az.getUncertainty())).append(")");
    if (!az.hadEnoughPoints()) {
      angleStr.append(" | WARNING: SMALL RANGE");
    }
    return angleStr.toString();
  }

  @Override
  public int panelsNeeded() {
    return 3;
  }

  @Override
  protected void updateData(DataStore ds) {

    set = true;

    double value = (double) offsetSpinner.getValue();

    if (value < 0) {
      value += 360;
    }

    AzimuthExperiment az = (AzimuthExperiment) expResult;
    az.setOffset(value);

    XYPlot xyp;

    expResult.runExperimentOnData(ds);
    List<XYSeriesCollection> allData = az.getData();
    XYSeriesCollection polars = allData.get(0);

    /*
    XYSeriesCollection xysc = allData.get(3); // coherence per-frequency
    coherenceChart = ChartFactory.createXYLineChart(
        expType.getName() + " Coherence", "Frequency (Hz)", "Coherence", xysc);
    */

    angleChart = ChartFactory.createPolarChart(expType.getName(),
        polars, true, true, false);

    String angleStr = getInsetStrings();

    /*
    XYPlot xyp = (XYPlot) coherenceChart.getPlot();
    TextTitle title = new TextTitle(angleStr);
    title.setBackgroundPaint(Color.white);
    XYTitleAnnotation xyt = new XYTitleAnnotation(0.98, 0.02, title,
        RectangleAnchor.BOTTOM_RIGHT);
    xyp.clearAnnotations();
    xyp.addAnnotation(xyt);
    // plot.addCornerTextItem(angleStr);
    */

    PolarPlot plot = (PolarPlot) angleChart.getPlot();
    plot.clearCornerTextItems();
    plot.addCornerTextItem(angleStr);

    XYSeriesCollection angleEstim = allData.get(1);
    XYSeriesCollection coherEstim = allData.get(2);
    String titleEst = expType.getName() + " Accuracy Estimation";
    estimChart = ChartFactory.createXYLineChart(titleEst,
        "xaxis", "yaxis", angleEstim);
    xyp = estimChart.getXYPlot();
    xyp.setDataset(0, angleEstim);
    xyp.setDataset(1, coherEstim);
    xyp.setRenderer(0, new DefaultXYItemRenderer());
    // set color of second dataset to be blue
    XYItemRenderer renderer = new DefaultXYItemRenderer();
    renderer.setSeriesPaint(0, Color.BLUE);
    xyp.setRenderer(1, renderer);
    NumberAxis angleEstimationAxis = new NumberAxis("Angle est. (deg)");
    // angleEstimationAxis.setAutoRangeIncludesZero(false);
    xyp.setRangeAxis(0, angleEstimationAxis);
    NumberAxis correlationAxis = new NumberAxis("Correlation est. of best fit angle");
    // correlationAxis.setAutoRangeIncludesZero(false);
    xyp.setRangeAxis(1, correlationAxis);
    NumberAxis xAx = new NumberAxis("Time from data start of (2000s) window (s)");
    xAx.setAutoRangeIncludesZero(false);
    xyp.setDomainAxis(xAx);
    xyp.mapDatasetToRangeAxis(0, 0);
    xyp.mapDatasetToRangeAxis(1, 1);

    if (!az.hadEnoughPoints()) {
      xyp = estimChart.getXYPlot();
      TextTitle result = new TextTitle();
      result.setText("WARNING: NOT ENOUGH DATA FOR WINDOWED COHERENCE ESTIMATION");
      result.setBackgroundPaint(Color.red);
      result.setPaint(Color.white);
      XYTitleAnnotation xyt = new XYTitleAnnotation(0.5, 0.5, result,
          RectangleAnchor.CENTER);
      xyp.clearAnnotations();
      xyp.addAnnotation(xyt);
    } else {
      double cutOff = az.getMinCorr();
      Marker highWater = new ValueMarker(cutOff);
      highWater.setStroke(new BasicStroke((float) 1.5));
      highWater.setPaint(Color.BLACK);
      xyp.addRangeMarker(1, highWater, Layer.BACKGROUND);
    }

    chartSelector.setSelectedIndex(0);
  }

}
