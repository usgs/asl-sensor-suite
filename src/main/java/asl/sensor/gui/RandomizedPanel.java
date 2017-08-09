package asl.sensor.gui;

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
import org.jfree.chart.annotations.XYTitleAnnotation;
import org.jfree.chart.axis.LogarithmicAxis;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.title.TextTitle;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.RectangleAnchor;

import asl.sensor.experiment.ExperimentEnum;
import asl.sensor.experiment.RandomizedExperiment;
import asl.sensor.experiment.ResponseExperiment;
import asl.sensor.input.DataStore;
import asl.sensor.utils.NumericUtils;

/**
 * Panel to display results from a randomized calibration experiment.
 * This includes plots of response magnitude and argument, selectable from a
 * drop-down combo box on the panel.
 * The inclusion of two selectable plots means that overrides are necessary
 * to produce output of both plots when creating a report of the results,
 * and that the typical means of assigning the visible chart cannot be used.
 * @author akearns
 *
 */
public class RandomizedPanel extends ExperimentPanel {

  public static final String MAGNITUDE = ResponseExperiment.MAGNITUDE;
  public static final String ARGUMENT = ResponseExperiment.ARGUMENT;
  private static final Color[] COLOR_LIST = 
      new Color[]{Color.RED, Color.BLUE, Color.GREEN};
  
  /**
   * 
   */
  private static final long serialVersionUID = -1791709117080520178L;
  /**
   * Utility function for formatting additional report pages from the
   * underlying experiment backend; can be called without constructing a
   * panel. Called by a non-static function in order to implement overrides, as
   * static functions do not get overridden by inheritance.
   * @param rnd RandomizedExperiment to pull data from (i.e., from a panel 
   * instance)
   * @return List of strings, each one representing a new page's worth of data
   */
  public static String[] getAdditionalReportPages(RandomizedExperiment rnd) {
    
    // TODO: refactor this now that period values are included in
    // the inset portion of the report text instead of merely in the extra data
    StringBuilder sb = new StringBuilder();
    
    StringBuilder csvPoles = new StringBuilder();
    StringBuilder csvZeros = new StringBuilder();
    StringBuilder csvTitle = new StringBuilder();
    DecimalFormat csvFormat = new DecimalFormat("+#.####;-#.####");
    NumericUtils.setInfinityPrintable(csvFormat);
    
    final int COL_WIDTH = 9;
    String[] columns = new String[]{"Init", "Fit", "Diff", "Mean", "PctDiff"};
    for (String column : columns) {
      StringBuilder paddedColumn = new StringBuilder(column);
      while ( paddedColumn.length() < COL_WIDTH ) {
        paddedColumn.append(" "); // add a space
      }
      csvTitle.append( paddedColumn );
    }
    
    List<Complex> fitP = rnd.getFitPoles();
    List<Complex> initP = rnd.getInitialPoles();
    List<Complex> fitZ = rnd.getFitZeros();
    List<Complex> initZ = rnd.getInitialZeros();
    
    boolean solverNotRun = rnd.getSolverState();
    
    if (solverNotRun) {
      return new String[]{};
    }
    
    // get statistics for differences between initial and solved parameters
    csvPoles = new StringBuilder("POLE VARIABLES, AS CSV:\n");
    csvPoles.append(csvTitle);
    csvPoles.append("\n");

    for (int i = 0; i < fitP.size(); ++i) {
      double realPartFit = fitP.get(i).getReal();
      double imagPartFit = fitP.get(i).getImaginary();

      double realPartInit = initP.get(i).getReal();
      double imagPartInit = initP.get(i).getImaginary();

      // make sure sign of the imaginary parts are the same
      if ( Math.signum(imagPartFit) != Math.signum(imagPartInit) ) {
        imagPartFit *= -1;
      }

      double realDiff = realPartInit - realPartFit;
      double imagDiff = imagPartInit - imagPartFit;

      double realAvg = (realPartInit + realPartFit) / 2.;
      double imagAvg = (imagPartInit + imagPartFit) / 2.;

      double realPct = realDiff * 100 / realPartFit;
      if ( realPartFit == 0. ) {
        realPct = 0.;
      }
      double imagPct = imagDiff * 100 / imagPartFit;
      if ( imagPartFit == 0. ) {
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
        while ( paddedColumn.length() < COL_WIDTH ) {
          paddedColumn.append(" "); // add a space
        }
        csvPoles.append( paddedColumn );
      }
      csvPoles.append("\n");

      for (double colNumber : imagRow) {
        String column = csvFormat.format(colNumber);
        StringBuilder paddedColumn = new StringBuilder(column);
        while ( paddedColumn.length() < COL_WIDTH ) {
          paddedColumn.append(" "); // add a space
        }
        csvPoles.append( paddedColumn );
      }
      csvPoles.append("\n");

      if (imagPartFit != 0.) {
        ++i; // skip complex conjugate
      }
    }

    // get statistics for differences between initial and solved parameters
    if ( fitZ.size() > 0 ) {
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
      if ( Math.signum(imagPartFit) != Math.signum(imagPartInit) ) {
        imagPartFit *= -1;
      }

      double realDiff = realPartInit - realPartFit;
      double imagDiff = imagPartInit - imagPartFit;

      double realAvg = (realPartInit + realPartFit) / 2.;
      double imagAvg = (imagPartInit + imagPartFit) / 2.;

      double realPct = realDiff * 100 / realPartFit;
      if ( realPartFit == 0. ) {
        realPct = 0.;
      }
      double imagPct = imagDiff * 100 / imagPartFit;
      if ( imagPartFit == 0. ) {
        imagPct = 0.;
      }

      double[] realRow = new double[]
          {realPartInit, realPartFit, realDiff, realAvg, realPct};

      double[] imagRow = new double[]
          {imagPartInit, imagPartFit, imagDiff, imagAvg, imagPct};

      for (double colNumber : realRow) {
        String column = csvFormat.format(colNumber);
        StringBuilder paddedColumn = new StringBuilder(column);
        while ( paddedColumn.length() < COL_WIDTH ) {
          paddedColumn.append(" "); // add a space
        }
        csvZeros.append( paddedColumn );
      }
      csvZeros.append("\n");

      for (double colNumber : imagRow) {
        String column = csvFormat.format(colNumber);
        StringBuilder paddedColumn = new StringBuilder(column);
        while ( paddedColumn.length() < COL_WIDTH ) {
          paddedColumn.append(" "); // add a space
        }
        csvZeros.append( paddedColumn );
      }
      csvZeros.append("\n");

      if (imagPartFit != 0.) {
        ++i; // skip complex conjugate
      }
    }

    sb.append(csvPoles);
    sb.append(csvZeros);

    String[] out = new String[]{sb.toString()}; // just a single new page
    return out;
  }
  /**
   * Static helper method for getting the formatted inset string directly
   * from a RandomizedExperiment
   * @param rnd RandomizedExperiment with data to be extracted
   * @return String format representation of data from the experiment
   */
  public static String getInsetString(RandomizedExperiment rnd) {
    
    final int MAX_LINE = 4; // maximum number of entries per line
    
    DecimalFormat df = new DecimalFormat("#.#####");
    NumericUtils.setInfinityPrintable(df);
    ComplexFormat cf = new ComplexFormat(df);
    
    List<Complex> fitP = rnd.getFitPoles();
    List<Complex> initP = rnd.getInitialPoles();
    List<Complex> fitZ = rnd.getFitZeros();
    List<Complex> initZ = rnd.getInitialZeros();
    
    if (fitP == null) {
      return "";
    }
    
    boolean solverNotRun = rnd.getSolverState();
    
    double initResid = rnd.getInitResidual();
    double fitResid = rnd.getFitResidual();
    
    StringBuilder sbInit = new StringBuilder();
    StringBuilder sbFit = new StringBuilder();
    // add poles, initial then fit (single loop, append the two builders)
    sbInit.append("Initial poles: \n");
    sbFit.append("Fit poles: \n");
    
    int numInLine = 0;
    
    for (int i = 0; i < fitP.size(); ++i) {
      
      Complex init = initP.get(i);
      Complex fit = fitP.get(i);
      double initPrd = NumericUtils.TAU / init.abs();
      double fitPrd = NumericUtils.TAU / fit.abs();
      
      sbInit.append( cf.format(init) );
      sbFit.append( cf.format(fit) );
      ++numInLine;
      // want to fit two to a line for paired values
      
      if ( init.getImaginary() != 0. ) {
        ++i; // INCREMENT I TO GET THE CONJUGATE AND NOT DO REDUNDANT OPERATION
        sbInit.append(";  ");
        sbFit.append(";  ");
        sbInit.append( cf.format( initP.get(i) ) );
        sbFit.append( cf.format( fitP.get(i) ) );
        sbInit.append(" (");
        sbInit.append( df.format(initPrd) );
        sbInit.append(" s)\n");
        sbFit.append(" (");
        sbFit.append( df.format(fitPrd) );
        sbFit.append(" s)\n");
        numInLine = 0;
      } else {
        // if data does not have conjugate pair, put period next to pole value
        // if there is still data, fit up to 4 poles in a line
        // but still separate conjugate pairs into their own line for space
        sbInit.append(" (");
        sbInit.append( df.format(initPrd) );
        sbInit.append(" s)");
        sbFit.append(" (");
        sbFit.append( df.format(fitPrd) );
        sbFit.append(" s)");
        
        if ( i + 1 < fitP.size() && 
             numInLine < MAX_LINE && 
             initP.get(i + 1).getImaginary() == 0. ) {
          // next value exists, is not part of a conjugate pair, and won't
          // put more than 4 pole values in a single line
          sbInit.append(";   ");
          sbFit.append(";   ");
        } else {
          sbInit.append("\n");
          sbFit.append("\n");
          numInLine = 0;
        }
      }
    }
    
    sbInit.append("\n");
    sbFit.append("\n");
    
    StringBuilder sbInitZ = new StringBuilder();
    StringBuilder sbFitZ = new StringBuilder();
    
    if ( fitZ.size() > 0 ) {
      sbInitZ.append("Initial zeros: \n");
      sbFitZ.append("Fit zeros: \n");
    }
    
    for (int i = 0; i < fitZ.size(); ++i) {
      
      Complex init = initZ.get(i);
      Complex fit = fitZ.get(i);
      double initPrd = NumericUtils.TAU / init.abs();
      double fitPrd = NumericUtils.TAU / fit.abs();
      
      sbInitZ.append( cf.format(init) );
      sbFitZ.append( cf.format(fit) );
      
      // want to fit two to a line for paired values
      if ( initZ.get(i).getImaginary() != 0. ) {
        ++i;
        sbInitZ.append("; ");
        sbFitZ.append("; ");
        sbInitZ.append( cf.format( initZ.get(i) ) );
        sbFitZ.append( cf.format( fitZ.get(i) ) );
        sbInitZ.append(" (");
        sbInitZ.append( df.format(initPrd) );
        sbInitZ.append(" s)\n");
        sbFitZ.append(" (");
        sbFitZ.append( df.format(fitPrd) );
        sbFitZ.append(" s)\n");
        
      } else {
        sbInitZ.append(" (");
        sbInitZ.append( df.format(initPrd) );
        sbInitZ.append(" s)");
        sbFitZ.append(" (");
        sbFitZ.append( df.format(fitPrd) );
        sbFitZ.append(" s)");
        if ( i + 1 < fitZ.size() ) {
          sbInitZ.append("\n");
          sbFitZ.append("\n");
        } else {
          sbInitZ.append(";   ");
          sbFitZ.append(";   ");
        }
      }
    }
    
    sbFit.append("\n");
    sbInit.append("\n");
    sbInitZ.append("\n");
    sbFitZ.append("\n");
    
    StringBuilder sb = new StringBuilder(sbInit);
    if (!solverNotRun) {
      sb.append(sbFit);
    }
    sb.append(sbInitZ);
    if (!solverNotRun) {
      sb.append(sbFitZ);
    }
    sb.append('\n');
    sb.append("Residuals:");
    sb.append('\n');
    sb.append("Initial (nom. resp curve): ");
    sb.append(initResid);
    sb.append('\n');
    if (!solverNotRun) {
      sb.append("Best fit: ");
      sb.append(fitResid);
    }
    
    return sb.toString();
  }
  private ValueAxis degreeAxis, residPhaseAxis, residAmpAxis;
  private JComboBox<String> plotSelection;
  private JCheckBox lowFreqBox, showParams, freqSpace;
  private JFreeChart magChart, argChart, residAmpChart, residPhaseChart;
  
  public RandomizedPanel(ExperimentEnum exp) {
    super(exp);
    
    channelType[0] = "Calibration input";
    channelType[1] = "Calibration output from sensor (RESP required)";
    
    String yAxisTitle = "10 * log10( RESP(f) )";
    String xAxisTitle = "Frequency (Hz)";
    String degreeAxisTitle = "phi(RESP(f))";
    
    xAxis = new LogarithmicAxis(xAxisTitle);
    
    yAxis = new NumberAxis(yAxisTitle);
    yAxis.setAutoRange(true);
    
    degreeAxis = new NumberAxis(degreeAxisTitle);
    degreeAxis.setAutoRange(true);    
    
    residPhaseAxis = new NumberAxis("Phase error (degrees)");
    residAmpAxis = new NumberAxis("Amplitude error (percentage)");
    
    ( (NumberAxis) yAxis).setAutoRangeIncludesZero(false);
    Font bold = xAxis.getLabelFont().deriveFont(Font.BOLD);
    xAxis.setLabelFont(bold);
    yAxis.setLabelFont(bold);
    degreeAxis.setLabelFont(bold);
    residPhaseAxis.setLabelFont(bold);
    residAmpAxis.setLabelFont(bold);
    
    lowFreqBox = new JCheckBox("Low frequency calibration");
    lowFreqBox.setSelected(true);
    
    showParams = new JCheckBox("Show params");
    showParams.setEnabled(false);
    showParams.addActionListener(this);
    
    freqSpace = new JCheckBox("Use Hz units (req. regen)");
    freqSpace.setSelected(true);
    
    applyAxesToChart(); // now that we've got axes defined
    
    magChart = buildChart(null, xAxis, yAxis);
    argChart = buildChart(null, xAxis, degreeAxis);
    residPhaseChart = buildChart(null, xAxis, residPhaseAxis);
    residAmpChart = buildChart(null, xAxis, residAmpAxis);
    
    // set the GUI components
    this.setLayout( new GridBagLayout() );
    GridBagConstraints gbc = new GridBagConstraints();
    
    gbc.fill = GridBagConstraints.BOTH;
    gbc.gridx = 0; gbc.gridy = 0;
    gbc.weightx = 1.0; gbc.weighty = 1.0;
    gbc.gridwidth = 3;
    gbc.anchor = GridBagConstraints.CENTER;
    this.add(chartPanel, gbc);
    
    // place the other UI elements in a single row below the chart
    gbc.gridwidth = 1;
    gbc.weighty = 0.0; gbc.weightx = 0.0;
    gbc.anchor = GridBagConstraints.WEST;
    gbc.fill = GridBagConstraints.NONE;
    gbc.gridy += 1; gbc.gridx = 0;
    JPanel checkBoxPanel = new JPanel();
    checkBoxPanel.setLayout( new BoxLayout(checkBoxPanel, BoxLayout.Y_AXIS) );
    checkBoxPanel.add(lowFreqBox);
    checkBoxPanel.add(showParams);
    checkBoxPanel.add(freqSpace);
    this.add(checkBoxPanel, gbc);
    
    gbc.gridx += 1;
    gbc.weightx = 1.0;
    gbc.fill = GridBagConstraints.NONE;
    gbc.anchor = GridBagConstraints.CENTER;
    // gbc.gridwidth = GridBagConstraints.REMAINDER;
    this.add(save, gbc);
    
    // plot selection combo box
    gbc.fill = GridBagConstraints.HORIZONTAL;
    gbc.gridx += 1;
    gbc.weightx = 0;
    gbc.anchor = GridBagConstraints.WEST;
    plotSelection = new JComboBox<String>();
    plotSelection.addItem(MAGNITUDE);
    plotSelection.addItem(ARGUMENT);
    plotSelection.addItem("Residual amplitude plot");
    plotSelection.addItem("Residual phase plot");
    plotSelection.addActionListener(this);
    this.add(plotSelection, gbc);
  }
  
  @Override
  public void actionPerformed(ActionEvent e) {
    
    super.actionPerformed(e);
    
    if ( e.getSource() == plotSelection ) {
      
      if (!set) {
        XYPlot xyp = chart.getXYPlot();
        String label = getXAxis().getLabel();
        xyp.getDomainAxis().setLabel(label);
        label = getYAxis().getLabel();
        xyp.getRangeAxis().setLabel(label);
        return;
      }
      
      int idx = plotSelection.getSelectedIndex();
      JFreeChart[] charts = 
          new JFreeChart[]{magChart, argChart, residAmpChart, residPhaseChart};
      chart = charts[idx];
      chartPanel.setChart(chart);
      
      return;
      
    }
    
    if ( e.getSource() == showParams ) {
      XYPlot xyp = magChart.getXYPlot();
      xyp.clearAnnotations();
      
      if ( showParams.isSelected() ) {    
        String inset = getInsetStrings();
        TextTitle result = new TextTitle();
        result.setText(inset);
        result.setBackgroundPaint(Color.white);
        double x;
        double y = 0.02;
        RectangleAnchor ra;
        // move text box left or right depending on which frequencies aren't
        // being fitted
        if ( lowFreqBox.isSelected() ) {
          x = 1;
          ra = RectangleAnchor.BOTTOM_RIGHT;
        } else {
          x = 0;
          ra = RectangleAnchor.BOTTOM_LEFT;
        }

        XYTitleAnnotation xyt = 
            new XYTitleAnnotation(x, y, result, ra);

        xyp.addAnnotation(xyt);
      }
      
      return;
    }
    
  }
  
  @Override
  protected void drawCharts() {
    
    // just force the active plot at the start to be the amplitude plot
    plotSelection.setSelectedIndex(0);
    showParams.setSelected(true);
    showParams.setEnabled(true);
    chart = magChart;
    chartPanel.setChart(chart);
    chartPanel.setMouseZoomable(true);
    
  }
  
  @Override
  public String[] getAdditionalReportPages() {
    // produce output of poles and zeros as period values in new report page
    
    RandomizedExperiment rnd = (RandomizedExperiment) expResult;
    
    return getAdditionalReportPages(rnd);
  }
  
  @Override
  public JFreeChart[] getCharts() {
    return new JFreeChart[]{magChart, argChart};
  }
  
  @Override
  /**
   * Get the index of the data holding the sensor output.
   * Note that the input data list is listed as CAL, OUT, RESP, so the
   * relevant index is the second one
   */
  protected int getIndexOfMainData() {
    return 1;
  }
  
  /**
   * Used to get the text that will populate the inset box for the plots
   * @return String to place in TextTitle
   */
  @Override
  public String getInsetStrings() {
    RandomizedExperiment rnd = (RandomizedExperiment) expResult;
    return getInsetString(rnd);
  }
  
  
  @Override
  public String getMetadataString() {
    RandomizedExperiment rnd = (RandomizedExperiment) expResult;
    StringBuilder sb = new StringBuilder();
    
    int iters = rnd.getIterations();
    sb.append("Iteration count from solver: ");
    sb.append(iters);
    sb.append("\n");
    
    sb.append( super.getMetadataString() );
    
    double[] weights = rnd.getWeights();
    sb.append("Residuals weighting:\n");
    sb.append("    Amplitude: ");
    sb.append(weights[0]);
    sb.append("\n");
    sb.append("    Phase: ");
    sb.append(weights[1]);
    return sb.toString();
  }
  
  @Override
  /**
   * Produce the filename of the report generated from this experiment.
   * Since response data is not directly associated with data at a given
   * time, rather than a sensor as a whole, we merely use the current date
   * and the first response used in the experiment.
   * @return String that will be default filename of PDF generated from data
   */
  public String getPDFFilename() {
    
    StringBuilder sb = new StringBuilder();
    if ( lowFreqBox.isSelected() ) {
      sb.append("Low_Frq_");
    } else {
      sb.append("High_Frq_");
    }
    
    sb.append( super.getPDFFilename() );
    
    return sb.toString();
  }
  
  @Override
  public JFreeChart[] getSecondPageCharts() {
    return new JFreeChart[]{residAmpChart, residPhaseChart};
  }
  
  @Override
  public ValueAxis getYAxis() {
    
    if ( null == plotSelection ) {
      return yAxis;
    }
    
    int idx = plotSelection.getSelectedIndex();
    ValueAxis[] out = 
        new ValueAxis[]{yAxis, degreeAxis, residAmpAxis, residPhaseAxis};
    return out[idx];
  }

  @Override
  public int panelsNeeded() {
    return 2;
  }

  @Override
  protected void updateData(DataStore ds) {
    
    set = true;
    showParams.setEnabled(false);
    
    final boolean isLowFreq = lowFreqBox.isSelected();
    seriesColorMap = new HashMap<String, Color>();
    
    RandomizedExperiment rndExp = (RandomizedExperiment) expResult;
    rndExp.setLowFreq(isLowFreq);
    rndExp.useFreqUnits( freqSpace.isSelected() );
    expResult.runExperimentOnData(ds);
    
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
      
      Color toColor = COLOR_LIST[i];
      String magName = (String) magSeries.getSeriesKey(i);
      seriesColorMap.put(magName, toColor);
      
      String argName = (String) argSeries.getSeriesKey(i);
      seriesColorMap.put(argName, toColor);
      
    }
    
    argChart = buildChart(argSeries, xAxis, degreeAxis);
    argChart.getXYPlot().getRangeAxis().setAutoRange(true);
    invertSeriesRenderingOrder( argChart );
    
    magChart = buildChart(magSeries, xAxis, yAxis);
    invertSeriesRenderingOrder( magChart );
    magChart.getXYPlot().getRangeAxis().setAutoRange(true);
    
    String inset = getInsetStrings();
    TextTitle result = new TextTitle();
    result.setText(inset);
    result.setBackgroundPaint(Color.white);
    double x;
    double y = 0.02;
    RectangleAnchor ra;
    // move text box left or right depending on which frequencies aren't
    // being fitted
    if ( lowFreqBox.isSelected() ) {
      x = 1;
      ra = RectangleAnchor.BOTTOM_RIGHT;
    } else {
      x = 0;
      ra = RectangleAnchor.BOTTOM_LEFT;
    }
    
    XYTitleAnnotation xyt = 
        new XYTitleAnnotation(x, y, result, ra);
    
    XYPlot xyp = magChart.getXYPlot();
    xyp.clearAnnotations();
    xyp.addAnnotation(xyt);
    
    appendChartTitle(argChart, appendFreqTitle);
    appendChartTitle(magChart, appendFreqTitle);
    
    // get residuals plot
    residAmpChart = buildChart(xysc.get(2), xAxis, residAmpAxis);
    /*
    double[] weights = rndExp.getWeights();
    StringBuilder sb = new StringBuilder();
    sb.append("Amplitude weighting: ");
    sb.append(weights[0]);
    sb.append("\nPhase weighting: ");
    sb.append(weights[1]);
    TextTitle weightInset = new TextTitle();
    weightInset.setText( sb.toString() );
    weightInset.setBackgroundPaint(Color.white);
    XYTitleAnnotation weightAnnot = 
        new XYTitleAnnotation(0, 1, weightInset, RectangleAnchor.TOP_LEFT);
    XYPlot residPlot = residChart.getXYPlot();
    residPlot.clearAnnotations();
    residPlot.addAnnotation(weightAnnot);
    */
    residPhaseChart = buildChart(xysc.get(3), xAxis, residPhaseAxis);
    /*
    double[] weights = rndExp.getWeights();
    StringBuilder sb = new StringBuilder();
    sb.append("Amplitude weighting: ");
    sb.append(weights[0]);
    sb.append("\nPhase weighting: ");
    sb.append(weights[1]);
    TextTitle weightInset = new TextTitle();
    weightInset.setText( sb.toString() );
    weightInset.setBackgroundPaint(Color.white);
    XYTitleAnnotation weightAnnot = 
        new XYTitleAnnotation(0, 1, weightInset, RectangleAnchor.TOP_LEFT);
    XYPlot residPlot = residChart.getXYPlot();
    residPlot.clearAnnotations();
    residPlot.addAnnotation(weightAnnot);
    */
  }

}
