package asl.sensor.gui;

import static asl.utils.ReportingUtils.setNonNumericPrintable;

import asl.sensor.ExperimentFactory;
import asl.sensor.experiment.RandomizedExperiment;
import asl.sensor.experiment.ResponseExperiment;
import asl.sensor.input.DataStore;
import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.text.DecimalFormat;
import java.util.HashMap;
import java.util.List;
import javax.swing.BoxLayout;
import javax.swing.ButtonGroup;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JSpinner;
import javax.swing.SpinnerModel;
import javax.swing.SpinnerNumberModel;
import javax.swing.SwingConstants;
import org.apache.commons.math3.complex.Complex;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.LogarithmicAxis;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.block.BlockContainer;
import org.jfree.chart.block.FlowArrangement;
import org.jfree.chart.plot.Marker;
import org.jfree.chart.plot.ValueMarker;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.chart.title.CompositeTitle;
import org.jfree.chart.title.LegendTitle;
import org.jfree.chart.title.TextTitle;
import org.jfree.chart.ui.RectangleEdge;
import org.jfree.chart.ui.VerticalAlignment;
import org.jfree.data.xy.XYSeriesCollection;

/**
 * Panel to display results from a randomized calibration experiment.
 * This includes plots of response magnitude and argument, selectable from a
 * drop-down combo box on the panel.
 * The inclusion of two selectable plots means that overrides are necessary
 * to produce output of both plots when creating a report of the results,
 * and that the typical means of assigning the visible chart cannot be used.
 *
 * @author akearns - KBRWyle
 */
public class RandomizedPanel extends ExperimentPanel {

  private static final long serialVersionUID = -1791709117080520178L;
  private final JComboBox<String> plotSelection;
  private final JSpinner nyquistMultiplier;
  private ValueAxis degreeAxis, residualPhaseAxis, residualAmplitudeAxis, periodAxis,
      residualXAxis, residualPeriodAxis;
  private JRadioButton lowFrequency, autoFrequency;
  // high frequency button created in constructor but not checked explicitly for status
  private JCheckBox showParams, frequencySpace, capacitiveCal;
  private JFreeChart magnitudeChart, argumentChart, residualAmplitudeChart, residualPhaseChart;

  public RandomizedPanel(ExperimentFactory experiment) {
    super(experiment);

    SpinnerModel spinModel = new SpinnerNumberModel(50., 30., 90., 1.);
    nyquistMultiplier = new JSpinner(spinModel);
    JLabel nyquistMultiplierLabel = new JLabel("% bound of nyquist for HF cals");
    nyquistMultiplierLabel.setLabelFor(nyquistMultiplier);
    nyquistMultiplierLabel.setHorizontalTextPosition(SwingConstants.RIGHT);
    nyquistMultiplierLabel.setHorizontalAlignment(SwingConstants.RIGHT);
    JPanel labelPanel = new JPanel();
    labelPanel.add(nyquistMultiplierLabel);

    JLabel introText = new JLabel("Calibration type: ");
    lowFrequency = new JRadioButton("low freq.");
    JRadioButton highFrequency = new JRadioButton("high freq.");
    autoFrequency = new JRadioButton("auto");

    ButtonGroup group = new ButtonGroup();
    group.add(lowFrequency);
    group.add(highFrequency);
    group.add(autoFrequency);
    autoFrequency.setSelected(true);

    channelType[0] = "Calibration input";
    channelType[1] = "Calibration output from sensor (RESP required)";

    initAxes();

    applyAxesToChart(); // now that we've got axes defined

    magnitudeChart = buildChart(null, xAxis, yAxis);
    argumentChart = buildChart(null, xAxis, degreeAxis);
    residualPhaseChart = buildChart(null, xAxis, residualPhaseAxis);
    residualAmplitudeChart = buildChart(null, xAxis, residualAmplitudeAxis);

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
    constraints.gridheight = 2;
    constraints.weighty = 0.0;
    constraints.weightx = 0.0;
    constraints.anchor = GridBagConstraints.WEST;
    constraints.fill = GridBagConstraints.NONE;
    constraints.gridy += 1;
    constraints.gridx = 0;
    JPanel optionsPanel = new JPanel();
    optionsPanel.setLayout(new BoxLayout(optionsPanel, BoxLayout.Y_AXIS));

    JPanel radioButtonPanel = new JPanel();
    radioButtonPanel.setAlignmentX(LEFT_ALIGNMENT);
    radioButtonPanel.setLayout(new BoxLayout(radioButtonPanel, BoxLayout.X_AXIS));
    radioButtonPanel.add(introText);
    radioButtonPanel.add(lowFrequency);
    radioButtonPanel.add(highFrequency);
    radioButtonPanel.add(autoFrequency);

    optionsPanel.add(radioButtonPanel);
    showParams.setAlignmentX(LEFT_ALIGNMENT);
    optionsPanel.add(showParams);
    frequencySpace.setAlignmentX(LEFT_ALIGNMENT);
    optionsPanel.add(frequencySpace);
    capacitiveCal.setAlignmentX(LEFT_ALIGNMENT);
    optionsPanel.add(capacitiveCal);
    this.add(optionsPanel, constraints);

    constraints.gridx += 1;
    constraints.weightx = 1.0;
    constraints.fill = GridBagConstraints.NONE;
    constraints.anchor = GridBagConstraints.SOUTH;
    this.add(save, constraints);

    // plot selection combo box
    constraints.fill = GridBagConstraints.BOTH;
    constraints.gridheight = 1;
    constraints.gridy = 1;
    constraints.gridx += 1;
    constraints.weightx = 0;
    constraints.anchor = GridBagConstraints.CENTER;
    JPanel nyquistLimitPanel = new JPanel();
    nyquistLimitPanel.setLayout(new BoxLayout(nyquistLimitPanel, BoxLayout.X_AXIS));
    nyquistLimitPanel.add(nyquistMultiplier);
    nyquistLimitPanel.add(nyquistMultiplierLabel);
    this.add(nyquistLimitPanel, constraints);

    constraints.fill = GridBagConstraints.NONE;
    constraints.gridy += 1;
    plotSelection = new JComboBox<>();
    plotSelection.addItem(ResponseExperiment.MAGNITUDE);
    plotSelection.addItem(ResponseExperiment.ARGUMENT);
    plotSelection.addItem("Residual amplitude plot");
    plotSelection.addItem("Residual phase plot");
    plotSelection.addActionListener(this);
    this.add(plotSelection, constraints);
  }

  /**
   * Utility function for formatting additional report pages from the
   * underlying experiment backend; can be called without constructing a
   * panel. Called by a non-static function in order to implement overrides, as
   * static functions do not get overridden by inheritance.
   *
   * @param experiment RandomizedExperiment to pull data from (i.e., from a panel instance)
   * @return List of strings, each one representing a new page's worth of data
   */
  public static String[] getAdditionalReportPages(RandomizedExperiment experiment) {

    // TODO: refactor this now that period values are included in
    // the inset portion of the report text instead of merely in the extra data
    StringBuilder resultString = new StringBuilder();

    StringBuilder csvPoles = new StringBuilder();
    StringBuilder csvZeros = new StringBuilder();
    StringBuilder csvTitle = new StringBuilder();
    DecimalFormat csvFormat = new DecimalFormat("+#.###;-#.###");
    setNonNumericPrintable(csvFormat);

    final int COL_WIDTH = 9;
    String[] columns = new String[]{"Init", "Fit", "Diff", "Mean", "PctDiff"};
    for (String column : columns) {
      StringBuilder paddedColumn = new StringBuilder(column);
      while (paddedColumn.length() < COL_WIDTH) {
        paddedColumn.append(" ");
      }
      csvTitle.append(paddedColumn);
    }

    List<Complex> fitP = experiment.getFitPoles();
    List<Complex> initP = experiment.getInitialPoles();
    List<Complex> fitZ = experiment.getFitZeros();
    List<Complex> initZ = experiment.getInitialZeros();

    // get statistics for differences between initial and solved parameters
    csvPoles.append("POLE VARIABLES, AS CSV:\n").append(csvTitle).append("\n");

    for (int i = 0; i < fitP.size(); ++i) {
      double realPartFit = fitP.get(i).getReal();
      double imagPartFit = fitP.get(i).getImaginary();

      double realPartInit = initP.get(i).getReal();
      double imagPartInit = initP.get(i).getImaginary();

      // make sure sign of the imaginary parts are the same
      if (Math.signum(imagPartFit) != Math.signum(imagPartInit)) {
        imagPartFit *= -1;
      }

      double realDiff = realPartInit - realPartFit;
      double imagDiff = imagPartInit - imagPartFit;

      double realAvg = (realPartInit + realPartFit) / 2.;
      double imagAvg = (imagPartInit + imagPartFit) / 2.;

      double realPct = realDiff * 100 / realPartFit;
      if (realPartFit == 0.) {
        realPct = 0.;
      }
      double imagPct = imagDiff * 100 / imagPartFit;
      if (imagPartFit == 0.) {
        imagPct = 0.;
      }

      // INIT, FIT, DIFF, AVG, PCT

      double[] realRow = new double[]
          {realPartInit, realPartFit, realDiff, realAvg, realPct};

      double[] imagRow = new double[]
          {imagPartInit, imagPartFit, imagDiff, imagAvg, imagPct};

      for (double colNumber : realRow) {
        String column = csvFormat.format(colNumber);
        StringBuilder paddedColumn = new StringBuilder(column);
        while (paddedColumn.length() < COL_WIDTH) {
          paddedColumn.append(" "); // add a space
        }
        csvPoles.append(paddedColumn);
      }
      csvPoles.append("\n");

      for (double colNumber : imagRow) {
        String column = csvFormat.format(colNumber);
        StringBuilder paddedColumn = new StringBuilder(column);
        while (paddedColumn.length() < COL_WIDTH) {
          paddedColumn.append(" "); // add a space
        }
        csvPoles.append(paddedColumn);
      }
      csvPoles.append("\n");

      if (imagPartFit != 0.) {
        ++i; // skip complex conjugate
      }
    }

    // get statistics for differences between initial and solved parameters
    if (fitZ.size() > 0) {
      csvZeros = new StringBuilder("ZERO VARIABLES, AS CSV:\n");
      csvZeros.append(csvTitle);
      csvZeros.append("\n");
    }

    for (int i = 0; i < fitZ.size(); ++i) {
      double realPartFit = fitZ.get(i).getReal();
      double imagPartFit = fitZ.get(i).getImaginary();

      double realPartInit = initZ.get(i).getReal();
      double imagPartInit = initZ.get(i).getImaginary();

      // make sure sign of the imaginary parts are the same
      if (Math.signum(imagPartFit) != Math.signum(imagPartInit)) {
        imagPartFit *= -1;
      }

      double realDiff = realPartInit - realPartFit;
      double imagDiff = imagPartInit - imagPartFit;

      double realAvg = (realPartInit + realPartFit) / 2.;
      double imagAvg = (imagPartInit + imagPartFit) / 2.;

      double realPct = realDiff * 100 / realPartFit;
      if (realPartFit == 0.) {
        realPct = 0.;
      }
      double imagPct = imagDiff * 100 / imagPartFit;
      if (imagPartFit == 0.) {
        imagPct = 0.;
      }

      double[] realRow = new double[]
          {realPartInit, realPartFit, realDiff, realAvg, realPct};

      double[] imagRow = new double[]
          {imagPartInit, imagPartFit, imagDiff, imagAvg, imagPct};

      for (double colNumber : realRow) {
        String column = csvFormat.format(colNumber);
        StringBuilder paddedColumn = new StringBuilder(column);
        while (paddedColumn.length() < COL_WIDTH) {
          paddedColumn.append(" "); // add a space
        }
        csvZeros.append(paddedColumn);
      }
      csvZeros.append("\n");

      for (double colNumber : imagRow) {
        String column = csvFormat.format(colNumber);
        StringBuilder paddedColumn = new StringBuilder(column);
        while (paddedColumn.length() < COL_WIDTH) {
          paddedColumn.append(" "); // add a space
        }
        csvZeros.append(paddedColumn);
      }
      csvZeros.append("\n");

      if (imagPartFit != 0.) {
        ++i; // skip complex conjugate
      }
    }

    resultString.append(csvPoles);
    resultString.append(csvZeros);

    return new String[]{resultString.toString()};
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

      int index = plotSelection.getSelectedIndex();
      JFreeChart[] charts =
          new JFreeChart[]{magnitudeChart, argumentChart, residualAmplitudeChart,
              residualPhaseChart};
      chart = charts[index];
      chartPanel.setChart(chart);
      return;
    }

    if (event.getSource() == showParams) {
      if (!showParams.isSelected()) {
        for (JFreeChart chart : getCharts()) {
          LegendTitle legend = chart.getLegend();
          chart.clearSubtitles();
          chart.addLegend(legend);
        }
      }

      if (showParams.isSelected()) {
        setSubtitles();
      }
    }
  }

  @Override
  protected void drawCharts() {
    // disable the listener here so that this selection override doesn't cause a race condition
    // since it's not being completed on the dispatch thread
    showParams.removeActionListener(this);
    showParams.setSelected(true);
    showParams.setEnabled(true);
    setSubtitles();
    showParams.addActionListener(this);
    // just force the active plot at the start to be the amplitude plot
    chart = magnitudeChart;
    chartPanel.setChart(chart);
    plotSelection.setSelectedIndex(0);
    chartPanel.setMouseZoomable(true);
  }

  @Override
  public String[] getAdditionalReportPages() {
    // produce output of poles and zeros as period values in new report page
    return getAdditionalReportPages((RandomizedExperiment) expResult);
  }

  @Override
  public JFreeChart[] getCharts() {
    return new JFreeChart[]{magnitudeChart, argumentChart};
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

  @Override
  public String getMetadataString() {
    RandomizedExperiment experiment = (RandomizedExperiment) expResult;
    StringBuilder result = new StringBuilder();

    int iters = experiment.getIterations();
    result.append("Iteration count from solver: ");
    result.append(iters);
    result.append("\n");

    result.append(super.getMetadataString());

    double[] weights = experiment.getWeights();
    result.append("Residuals weighting:\n");
    result.append("    Amplitude: ");
    result.append(weights[0]);
    result.append("\n");
    result.append("    Phase: ");
    result.append(weights[1]);
    return result.toString();
  }

  /**
   * Produce the filename of the report generated from this experiment.
   * Since response data is not directly associated with data at a given
   * time, rather than a sensor as a whole, we merely use the current date
   * and the first response used in the experiment.
   *
   * @return String that will be default filename of PDF generated from data
   */
  @Override
  public String getPDFFilename() {
    RandomizedExperiment rExp = (RandomizedExperiment) expResult;
    if (rExp.isLowFrequencyCalibration()) {
      return "Low_Frq_" + super.getPDFFilename();
    } else {
      return "High_Frq_" + super.getPDFFilename();
    }
  }

  private ValueAxis getResidAxis() {
    if (null == plotSelection || frequencySpace.isSelected()) {
      return residualXAxis;
    } else {
      return residualPeriodAxis;
    }
  }

  @Override
  public JFreeChart[] getSecondPageCharts() {
    return new JFreeChart[]{residualAmplitudeChart, residualPhaseChart};
  }

  @Override
  public ValueAxis getXAxis() {
    if (null == plotSelection || frequencySpace.isSelected()) {
      return xAxis;
    } else {
      return periodAxis;
    }
  }

  @Override
  public ValueAxis getYAxis() {
    if (null == plotSelection) {
      return yAxis;
    }

    int index = plotSelection.getSelectedIndex();
    ValueAxis[] out =
        new ValueAxis[]{yAxis, degreeAxis, residualAmplitudeAxis, residualPhaseAxis};
    return out[index];
  }

  private void initAxes() {
    String yAxisTitle = "20 * log10( RESP(f) )";
    String xAxisTitle = "Frequency (Hz)";
    String prdAxisTitle = "Period (s)";
    String degreeAxisTitle = "phi(RESP(f))";

    xAxis = new LogarithmicAxis(xAxisTitle);
    periodAxis = new LogarithmicAxis(prdAxisTitle);
    residualXAxis = new LogarithmicAxis(xAxisTitle);
    residualPeriodAxis = new LogarithmicAxis(prdAxisTitle);

    yAxis = new NumberAxis(yAxisTitle);
    yAxis.setAutoRange(true);

    degreeAxis = new NumberAxis(degreeAxisTitle);
    degreeAxis.setAutoRange(true);

    residualPhaseAxis = new NumberAxis("Phase error (percentage)");
    residualAmplitudeAxis = new NumberAxis("Amplitude error (percentage)");

    ((NumberAxis) yAxis).setAutoRangeIncludesZero(false);
    Font bold = xAxis.getLabelFont();
    bold = bold.deriveFont(Font.BOLD, bold.getSize() + 2);
    xAxis.setLabelFont(bold);
    yAxis.setLabelFont(bold);
    degreeAxis.setLabelFont(bold);
    residualPhaseAxis.setLabelFont(bold);
    residualAmplitudeAxis.setLabelFont(bold);

    showParams = new JCheckBox("Show params");
    showParams.setEnabled(false);
    showParams.addActionListener(this);

    frequencySpace = new JCheckBox("Use Hz units (req. regen)");
    frequencySpace.setSelected(true);

    capacitiveCal = new JCheckBox("Capacitive calibration");
    capacitiveCal.setSelected(false);
  }

  @Override
  public int panelsNeeded() {
    return 2;
  }

  private void setSubtitles() {
    BlockContainer bc = new BlockContainer(new FlowArrangement());
    CompositeTitle ct = new CompositeTitle(bc);
    String[] insets = expResult.getInsetStrings();
    for (String inset : insets) {
      TextTitle result = getDefaultTextTitle();
      result.setText(inset);
      bc.add(result);
    }

    TextTitle result = getDefaultTextTitle();
    RandomizedExperiment re = (RandomizedExperiment) expResult;
    int numIters = re.getIterations();
    result.setText("NUMBER OF ITERATIONS: " + numIters);

    ct.setVerticalAlignment(VerticalAlignment.BOTTOM);
    ct.setPosition(RectangleEdge.BOTTOM);
    result.setVerticalAlignment(VerticalAlignment.BOTTOM);
    result.setPosition(RectangleEdge.BOTTOM);
    for (JFreeChart chart : getCharts()) {
      chart.addSubtitle(ct);
      chart.addSubtitle(result);
    }
  }

  @Override
  protected void updateData(DataStore dataStore) {
    set = true;
    showParams.setSelected(false);

    final boolean isLowFreq = lowFrequency.isSelected();
    seriesColorMap = new HashMap<>();

    // we display as % but backend expects this to be used
    double multiplier = (double) nyquistMultiplier.getValue() / 100.;

    RandomizedExperiment rndExp = (RandomizedExperiment) expResult;
    if (autoFrequency.isSelected()) {
      rndExp.autoDetermineCalibrationStatus(dataStore);
    } else {
      rndExp.setLowFrequencyCalibration(isLowFreq);
    }
    rndExp.setPlotUsingHz(frequencySpace.isSelected());
    rndExp.setNyquistMultiplier(multiplier);
    rndExp.setCapactiveCalibration(capacitiveCal.isSelected());
    expResult.runExperimentOnData(dataStore);

    String appendFreqTitle;

    if (isLowFreq) {
      appendFreqTitle = " (LOW FREQ.)";
    } else {
      appendFreqTitle = " (HIGH FREQ.)";
    }

    List<XYSeriesCollection> xysc = expResult.getData();

    XYSeriesCollection magSeries = xysc.get(0);
    XYSeriesCollection argSeries = xysc.get(1);

    for (int i = 0; i < magSeries.getSeriesCount(); ++i) {

      Color toColor = getColor(i);
      String magName = (String) magSeries.getSeriesKey(i);
      seriesColorMap.put(magName, toColor);

      String argName = (String) argSeries.getSeriesKey(i);
      seriesColorMap.put(argName, toColor);

    }

    getXAxis().setAutoRange(true);

    argumentChart = buildChart(argSeries, getXAxis(), degreeAxis);
    argumentChart.getXYPlot().getRangeAxis().setAutoRange(true);
    invertSeriesRenderingOrder(argumentChart);

    magnitudeChart = buildChart(magSeries, getXAxis(), yAxis);
    magnitudeChart.getXYPlot().getRangeAxis().setAutoRange(true);
    invertSeriesRenderingOrder(magnitudeChart);

    if (!isLowFreq) {
      Marker maxFitMarker = new ValueMarker(rndExp.getMaxFitFrequency());
      maxFitMarker.setStroke(new BasicStroke((float) 1.5));
      maxFitMarker.setPaint(Color.BLACK);
      magnitudeChart.getXYPlot().addDomainMarker(maxFitMarker);
      argumentChart.getXYPlot().addDomainMarker(maxFitMarker);
    }

    appendChartTitle(argumentChart, appendFreqTitle);
    appendChartTitle(magnitudeChart, appendFreqTitle);

    // get residuals plots
    XYItemRenderer renderer; // use this to set series paints to set the last color to use 3rd color
    residualAmplitudeChart = buildChart(xysc.get(2), getResidAxis(), residualAmplitudeAxis);
    renderer = residualAmplitudeChart.getXYPlot().getRenderer();
    renderer.setSeriesPaint(1, getColor(2));
    residualPhaseChart = buildChart(xysc.get(3), getResidAxis(), residualPhaseAxis);
    renderer = residualPhaseChart.getXYPlot().getRenderer();
    renderer.setSeriesPaint(1, getColor(2));

  }

}
