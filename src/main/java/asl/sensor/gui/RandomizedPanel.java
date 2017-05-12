package asl.sensor.gui;

import java.awt.Color;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.text.DecimalFormat;
import java.util.HashMap;
import java.util.List;

import javax.swing.JCheckBox;
import javax.swing.JComboBox;

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
    StringBuilder sb = new StringBuilder();
    DecimalFormat df = new DecimalFormat("#.#####");
    
    List<Complex> fitP = rnd.getFitPoles();
    List<Complex> initP = rnd.getInitialPoles();
    List<Complex> fitZ = rnd.getFitZeros();
    List<Complex> initZ = rnd.getInitialZeros();
    
    StringBuilder initText = new StringBuilder("Initial:\n");
    StringBuilder fitText = new StringBuilder("Best fit:\n");
    
    sb.append("Pole and zero values, given as period (s):\n \n");
    sb.append("Poles:\n");
    for (int i = 0; i < fitP.size(); ++i) {
      double fitPrd = NumericUtils.TAU / fitP.get(i).abs();
      double initPrd = NumericUtils.TAU / initP.get(i).abs();
      
      fitText.append( df.format(fitPrd) );
      initText.append( df.format(initPrd) );
      fitText.append("\n");
      initText.append("\n");
      
      if ( fitP.get(i).getImaginary() != 0. ) {
        // complex conjugate pole is at same period value, don't report
        ++i;
      }
    }
    
    sb.append(initText);
    sb.append(fitText);
        
    initText = new StringBuilder();
    fitText = new StringBuilder();
    
    if ( fitZ.size() > 0 ) {
      sb.append(" \nZeros:\n");
      initText.append("Initial:\n");
      fitText.append("Best fit:\n");
    }

    for (int i = 0; i < fitZ.size(); ++i) {
      double fitPrd = NumericUtils.TAU / fitZ.get(i).abs();
      double initPrd = NumericUtils.TAU / initZ.get(i).abs();
      
      fitText.append( df.format(fitPrd) );
      initText.append( df.format(initPrd) );
      fitText.append("\n");
      initText.append("\n");
      
      if ( fitZ.get(i).getImaginary() != 0. ) {
        // complex conjugate pole is at same period value, don't report
        ++i;
      }
    }
    
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
    
    DecimalFormat df = new DecimalFormat("#.#####");
    ComplexFormat cf = new ComplexFormat(df);
    
    List<Complex> fitP = rnd.getFitPoles();
    List<Complex> initP = rnd.getInitialPoles();
    List<Complex> fitZ = rnd.getFitZeros();
    List<Complex> initZ = rnd.getInitialZeros();
    
    double initResid = rnd.getInitResidual();
    double fitResid = rnd.getFitResidual();
    
    StringBuilder sbInit = new StringBuilder();
    StringBuilder sbFit = new StringBuilder();
    // add poles, initial then fit (single loop, append the two builders)
    sbInit.append("Initial poles: \n");
    sbFit.append("Fit poles: \n");
    for (int i = 0; i < fitP.size(); ++i) {
      sbInit.append( cf.format( initP.get(i) ) );
      sbFit.append( cf.format( fitP.get(i) ) );
      
      // want to fit two to a line for paired values
      if ( initP.get(i).getImaginary() != 0. ) {
        ++i;
        sbInit.append("; ");
        sbFit.append("; ");
        sbInit.append( cf.format( initP.get(i) ) );
        sbInit.append("  ");
        sbFit.append( cf.format( fitP.get(i) ) );
        sbFit.append("  ");
      }
      sbInit.append("\n");
      sbFit.append("\n");
    }
    
    StringBuilder sbInitZ = new StringBuilder();
    StringBuilder sbFitZ = new StringBuilder();
    
    if ( fitZ.size() > 0 ) {
      sbInitZ.append("Initial zeros: \n");
      sbFitZ.append("Fit zeros: \n");
    }
    
    for (int i = 0; i < fitZ.size(); ++i) {

      sbInitZ.append( cf.format( initZ.get(i) ) );
      sbFitZ.append( cf.format( fitZ.get(i) ) );
      
      // want to fit two to a line for paired values
      if ( initZ.get(i).getImaginary() != 0. ) {
        ++i;
        sbInitZ.append("; ");
        sbFitZ.append("; ");
        sbInitZ.append( cf.format( initZ.get(i) ) );
        sbInitZ.append("\n");
        sbFitZ.append( cf.format( fitZ.get(i) ) );
        sbFitZ.append("\n");
        
      } else {
        if ( i + 1 < fitZ.size() ) {
          sbInitZ.append(";   ");
          sbFitZ.append(";   ");
        }
      }
    }
    
    sbFit.append("\n");
    sbInit.append("\n");
    sbInitZ.append("\n");
    sbFitZ.append("\n");
    
    StringBuilder sb = new StringBuilder( sbInit.append(sbFit) );
    sb.append( sbInitZ.append(sbFitZ) );
    sb.append('\n');
    sb.append("Residuals:");
    sb.append('\n');
    sb.append("Initial (nom. resp curve): ");
    sb.append(initResid);
    sb.append('\n');
    sb.append("Best fit: ");
    sb.append(fitResid);
    
    return sb.toString();
  }
  private ValueAxis degreeAxis, residAxis;
  
  private JComboBox<String> plotSelection;

  private JCheckBox lowFreqBox;

  private JFreeChart magChart, argChart, residChart;
  
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
    
    residAxis = new NumberAxis("Residual value per point");
    
    ( (NumberAxis) yAxis).setAutoRangeIncludesZero(false);
    Font bold = xAxis.getLabelFont().deriveFont(Font.BOLD);
    xAxis.setLabelFont(bold);
    yAxis.setLabelFont(bold);
    residAxis.setLabelFont(bold);
    
    lowFreqBox = new JCheckBox("Low frequency calibration");
    lowFreqBox.setSelected(true);
    
    applyAxesToChart(); // now that we've got axes defined
    
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
    this.add(lowFreqBox, gbc);
    
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
    plotSelection.addItem("Residuals plot");
    this.add(plotSelection, gbc);
    plotSelection.addActionListener(this);
  }
  
  @Override
  public void actionPerformed(ActionEvent e) {
    
    super.actionPerformed(e);
    
    if ( e.getSource() == plotSelection ) {
      
      if (!set) {
        applyAxesToChart();
        return;
      }
      
      int idx = plotSelection.getSelectedIndex();
      JFreeChart[] charts = 
          new JFreeChart[]{magChart, argChart, residChart};
      chart = charts[idx];
      chartPanel.setChart(chart);
      
      return;
      
    }
    
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
  
  /**
   * Used to get the text that will populate the inset box for the plots
   * @return String to place in TextTitle
   */
  @Override
  public String getInsetString() {
    RandomizedExperiment rnd = (RandomizedExperiment) expResult;
    return getInsetString(rnd);
  }
  
  @Override
  public String getMetadataString() {
    RandomizedExperiment rnd = (RandomizedExperiment) expResult;
    StringBuilder sb = new StringBuilder();
    sb.append("LOADED RESPONSE:");
    sb.append('\n');
    sb.append( rnd.getResponseName() );
    return sb.toString();
  }
  
  @Override
  public ValueAxis getYAxis() {
    
    if ( null == plotSelection ) {
      return yAxis;
    }
    
    int idx = plotSelection.getSelectedIndex();
    ValueAxis[] out = new ValueAxis[]{yAxis, degreeAxis, residAxis};
    return out[idx];
  }
  
  @Override
  public int panelsNeeded() {
    return 2;
  }
  
  @Override
  public JFreeChart[] getSecondPageCharts() {
    return new JFreeChart[]{residChart};
  }

  @Override
  public void updateData(DataStore ds) {

    final boolean isLowFreq = lowFreqBox.isSelected();
    seriesColorMap = new HashMap<String, Color>();
    
    RandomizedExperiment respExp = (RandomizedExperiment) expResult;
    respExp.setLowFreq(isLowFreq);
    expResult.setData(ds);
    
    String appendFreqTitle;
    
    if (isLowFreq) {
      appendFreqTitle = " (LOW FREQ.)";
    } else {
      appendFreqTitle = " (HIGH FREQ.)";
    }
    
    set = true;
    
    clearChartAndSetProgressData();
    
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
    
    // Range argRange = argSeries.getRangeBounds(true);
    
    String inset = getInsetString();
    TextTitle result = new TextTitle();
    result.setText( inset );
    result.setBackgroundPaint(Color.white);

    argChart = buildChart(argSeries, xAxis, degreeAxis);
    argChart.getXYPlot().getRangeAxis().setAutoRange(true);
    
    magChart = buildChart(magSeries, xAxis, yAxis);
    magChart.getXYPlot().getRangeAxis().setAutoRange(true);
    
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
    
    residChart = buildChart(xysc.get(2), xAxis, residAxis);
    xyp = residChart.getXYPlot();
    
    double[] weights = respExp.getWeights();
    StringBuilder sb = new StringBuilder();
    sb.append("Magnitude weighting: ");
    sb.append(weights[0]);
    sb.append("\nPhase weighting: ");
    sb.append(weights[1]);
    TextTitle weightInset = new TextTitle();
    result.setText( sb.toString() );
    result.setBackgroundPaint(Color.white);
    xyp.clearAnnotations();
    XYTitleAnnotation weightAnnot = 
        new XYTitleAnnotation(0, 0, weightInset, RectangleAnchor.BOTTOM_LEFT);
    xyp.addAnnotation(weightAnnot);
    
    plotSelection.setSelectedIndex(0);
    chart = magChart;
    
    /*
    if (idx == 0) {
      int idx = plotSelection.getSelectedIndex();
      chart = magChart;
    } else {
      chart = argChart;
      // chart.getXYPlot().getRangeAxis().setRange(argRange); 
    }
    */
    chartPanel.setChart(chart);
    chartPanel.setMouseZoomable(true);
    
  }

}
