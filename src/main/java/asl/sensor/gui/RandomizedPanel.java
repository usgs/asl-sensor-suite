package asl.sensor.gui;

import asl.sensor.experiment.ExperimentEnum;
import asl.sensor.experiment.RandomizedExperiment;
import asl.sensor.experiment.ResponseExperiment;
import asl.sensor.input.DataStore;
import asl.sensor.utils.NumericUtils;
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
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JPanel;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.complex.ComplexFormat;
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
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.RectangleEdge;
import org.jfree.ui.VerticalAlignment;

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

  private static String complexListToString(List<Complex> complexList) {
    final int MAX_LINE = 2; // maximum number of entries per line

    ComplexFormat complexFormat = new ComplexFormat(DECIMAL_FORMAT.get());
    StringBuilder stringBuilder = new StringBuilder();
    int numInLine = 0;

    for (Complex number : complexList) {

      double initPrd = NumericUtils.TAU / number.abs();

      stringBuilder.append(complexFormat.format(number));
      stringBuilder.append(" (");
      stringBuilder.append(DECIMAL_FORMAT.get().format(initPrd));
      stringBuilder.append(")");
      ++numInLine;
      // want to fit two to a line for paired values
      if (numInLine >= MAX_LINE) {
        stringBuilder.append("\n");
        numInLine = 0;
      } else {
        stringBuilder.append(", ");
      }
    }

    return stringBuilder.toString();
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
    NumericUtils.setInfinityPrintable(csvFormat);

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

    boolean solverNotRun = experiment.getSolverState();

    if (solverNotRun) {
      return new String[]{};
    }

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


  /**
   * Static helper method for getting the formatted inset string directly
   * from a RandomizedExperiment
   *
   * @param experiment RandomizedExperiment with data to be extracted
   * @return String format representation of data from the experiment
   */
  public static String[] getInsetString(RandomizedExperiment experiment) {

    List<Complex> fitPoles = experiment.getFitPoles();
    List<Complex> initialPoles = experiment.getInitialPoles();
    List<Complex> fitZeros = experiment.getFitZeros();
    List<Complex> initialZeros = experiment.getInitialZeros();

    if (fitPoles == null) {
      return new String[]{""};
    }

    boolean solverNotRun = experiment.getSolverState();

    double initialResidual = experiment.getInitResidual();
    double fitResidual = experiment.getFitResidual();

    StringBuilder sbInitialPoles = new StringBuilder();
    StringBuilder sbFitPoles = new StringBuilder();
    // add poles, initial then fit (single loop, append the two builders)
    sbInitialPoles.append("Initial poles: \n");
    sbFitPoles.append("Fit poles: \n");

    sbInitialPoles.append(complexListToString(initialPoles));
    sbFitPoles.append(complexListToString(fitPoles));

    sbInitialPoles.append("\n");
    sbFitPoles.append("\n");

    StringBuilder sbInitZ = new StringBuilder();
    StringBuilder sbFitZ = new StringBuilder();

    if (fitZeros.size() > 0) {
      sbInitZ.append("Initial zeros: \n");
      sbFitZ.append("Fit zeros: \n");
    }

    sbInitZ.append(complexListToString(initialZeros));
    sbFitZ.append(complexListToString(fitZeros));

    sbFitPoles.append("\n");
    sbInitialPoles.append("\n");
    sbInitZ.append("\n");
    sbFitZ.append("\n");

    if (!solverNotRun) {
      sbInitialPoles.append(sbFitPoles);
    }

    if (!solverNotRun) {
      sbInitZ.append(sbFitZ);
    }

    StringBuilder sbR = new StringBuilder();
    sbR.append("Residuals:");
    sbR.append('\n');
    sbR.append("Initial (nom. resp curve): ");
    sbR.append(initialResidual);
    sbR.append('\n');
    if (!solverNotRun) {
      sbR.append("Best fit: ");
      sbR.append(fitResidual);
    }

    return new String[]{sbInitialPoles.toString(), sbInitZ.toString(), sbR.toString()};
  }

  private ValueAxis degreeAxis, residualPhaseAxis, residualAmplitudeAxis, periodAxis,
      residualXAxis, residualPeriodAxis;
  private final JComboBox<String> plotSelection;
  private JCheckBox lowFrequencyBox, showParams, frequencySpace;
  private JFreeChart magnitudeChart, argumentChart, residualAmplitudeChart, residualPhaseChart;

  public RandomizedPanel(ExperimentEnum experiment) {
    super(experiment);

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
    constraints.weighty = 0.0;
    constraints.weightx = 0.0;
    constraints.anchor = GridBagConstraints.WEST;
    constraints.fill = GridBagConstraints.NONE;
    constraints.gridy += 1;
    constraints.gridx = 0;
    JPanel checkBoxPanel = new JPanel();
    checkBoxPanel.setLayout(new BoxLayout(checkBoxPanel, BoxLayout.Y_AXIS));
    checkBoxPanel.add(lowFrequencyBox);
    checkBoxPanel.add(showParams);
    checkBoxPanel.add(frequencySpace);
    this.add(checkBoxPanel, constraints);

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
    plotSelection = new JComboBox<>();
    plotSelection.addItem(ResponseExperiment.MAGNITUDE);
    plotSelection.addItem(ResponseExperiment.ARGUMENT);
    plotSelection.addItem("Residual amplitude plot");
    plotSelection.addItem("Residual phase plot");
    plotSelection.addActionListener(this);
    this.add(plotSelection, constraints);
  }

  @Override
  public void actionPerformed(ActionEvent e) {
    super.actionPerformed(e);

    if (e.getSource() == plotSelection) {
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

    if (e.getSource() == showParams) {
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
    // just force the active plot at the start to be the amplitude plot
    showParams.setSelected(true);
    showParams.setEnabled(true);
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

  /**
   * Used to get the text that will represent the title text in the PDF result
   */
  @Override
  public String getInsetStrings() {
    StringBuilder sb = new StringBuilder();
    for (String str : getInsetStringsAsList()) {
      sb.append(str);
      sb.append("\n");
    }
    return sb.toString();
  }

  /**
   * Produce arrays of pole, zero, and residual data for text titles
   *
   * @return Array of strings
   */
  private String[] getInsetStringsAsList() {
    return getInsetString((RandomizedExperiment) expResult);
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
    if (lowFrequencyBox.isSelected()) {
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

    residualPhaseAxis = new NumberAxis("Phase error (degrees)");
    residualAmplitudeAxis = new NumberAxis("Amplitude error (percentage)");

    ((NumberAxis) yAxis).setAutoRangeIncludesZero(false);
    Font bold = xAxis.getLabelFont().deriveFont(Font.BOLD);
    xAxis.setLabelFont(bold);
    yAxis.setLabelFont(bold);
    degreeAxis.setLabelFont(bold);
    residualPhaseAxis.setLabelFont(bold);
    residualAmplitudeAxis.setLabelFont(bold);

    lowFrequencyBox = new JCheckBox("Low frequency calibration");
    lowFrequencyBox.setSelected(true);

    showParams = new JCheckBox("Show params");
    showParams.setEnabled(false);
    showParams.addActionListener(this);

    frequencySpace = new JCheckBox("Use Hz units (req. regen)");
    frequencySpace.setSelected(true);
  }

  @Override
  public int panelsNeeded() {
    return 2;
  }

  private void setSubtitles() {
    BlockContainer bc = new BlockContainer(new FlowArrangement());
    CompositeTitle ct = new CompositeTitle(bc);
    String[] insets = getInsetStringsAsList();
    for (String inset : insets) {
      TextTitle result = new TextTitle();
      result.setText(inset);
      result.setBackgroundPaint(Color.white);
      bc.add(result);
    }

    TextTitle result = new TextTitle();
    RandomizedExperiment re = (RandomizedExperiment) expResult;
    int numIters = re.getIterations();
    result.setText("NUMBER OF ITERATIONS: " + numIters);
    result.setBackgroundPaint(Color.white);

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

    final boolean isLowFreq = lowFrequencyBox.isSelected();
    seriesColorMap = new HashMap<>();

    RandomizedExperiment rndExp = (RandomizedExperiment) expResult;
    rndExp.setLowFreq(isLowFreq);
    rndExp.useFreqUnits(frequencySpace.isSelected());
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

      Color toColor = COLORS[i % COLORS.length];
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

    String inset = getInsetStrings();
    TextTitle result = new TextTitle();
    result.setText(inset);
    result.setBackgroundPaint(Color.white);

    appendChartTitle(argumentChart, appendFreqTitle);
    appendChartTitle(magnitudeChart, appendFreqTitle);

    // get residuals plots
    XYItemRenderer renderer; // use this to set series paints to set the last color to use 3rd color
    residualAmplitudeChart = buildChart(xysc.get(2), getResidAxis(), residualAmplitudeAxis);
    renderer = residualAmplitudeChart.getXYPlot().getRenderer();
    renderer.setSeriesPaint(1, COLORS[2]);
    residualPhaseChart = buildChart(xysc.get(3), getResidAxis(), residualPhaseAxis);
    renderer = residualPhaseChart.getXYPlot().getRenderer();
    renderer.setSeriesPaint(1, COLORS[2]);

    setSubtitles();
  }

}
