package asl.sensor.gui;

import java.awt.Color;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.text.DecimalFormat;

import javax.swing.JComboBox;
import javax.swing.JPanel;

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
import asl.sensor.experiment.StepExperiment;
import asl.sensor.input.DataStore;

/**
 * Holds the plot results of a step experiment. Gets the timeseries data from it
 * as well as the corner and damping values gotten in the process of its
 * calculations. Aside from common elements to ExperimentPanels, also sets an
 * inset displaying parameters derived from backend fit calculation
 * @author akearns
 *
 */
public class StepPanel extends ExperimentPanel {

  /**
   * 
   */
  private static final long serialVersionUID = 3693391540945130688L;

  /**
   * Static helper method for getting the formatted inset string directly
   * from a StepExperiment
   * @param sp StepExperiment with data to be extracted
   * @return String format representation of data from the experiment
   */
  public static String getInsetString(StepExperiment sp) {  
    double[] rolloff = sp.getInitParams();
    double[] fit = sp.getFitParams();
    double corner = rolloff[0];
    double damping = rolloff[1];
    double fitCorner = fit[0];
    double fitDamping = fit[1];
    
    double cornerPrd = 1. / corner; 
    double fitCornerPrd = 1. / corner;
    
    DecimalFormat df = new DecimalFormat("#.######");
    
    
    StringBuilder sb = new StringBuilder();
    sb.append("RESP parameters\n");
    sb.append("corner frequency (Hz): ");
    sb.append( df.format(corner) );
    sb.append(" (");
    sb.append( df.format(cornerPrd) );
    sb.append( " secs)");
    sb.append("\n");
    sb.append("damping: ");
    sb.append( df.format(damping) );
    sb.append("\n");
    sb.append("Best-fit parameters\n");
    sb.append("corner frequency (Hz): ");
    sb.append( df.format(fitCorner) );
    sb.append(" (");
    sb.append( df.format(fitCornerPrd) );
    sb.append( " secs)");
    sb.append("\n");
    sb.append("damping: ");
    sb.append( df.format(fitDamping) );
    sb.append("\n");
    return sb.toString();
    
  }
  
  private JComboBox<String> plotSelection;
  private JFreeChart stepChart, magChart, phaseChart;
  private ValueAxis freqAxis, magAxis, phaseAxis;

  public StepPanel(ExperimentEnum exp) {
    super(exp);
    
    channelType[0] = "Calibration input";
    channelType[1] = "Calibration output from sensor (RESP required)";
    
    String xAxisTitle = "Time (s)";
    String yAxisTitle = "Normalized counts";
    xAxis = new NumberAxis(xAxisTitle);
    Font bold = xAxis.getLabelFont().deriveFont(Font.BOLD);
    xAxis.setLabelFont(bold);
    // xAxis.setAutoRange(true);
    
    yAxis = new NumberAxis(yAxisTitle);
    yAxis.setLabelFont(bold);
    // yAxis.setAutoRange(true);
    
    String freqAxisTitle = "Frequency (Hz)";
    freqAxis = new LogarithmicAxis(freqAxisTitle);
    freqAxis.setLabelFont(bold);
    freqAxis.setAutoRange(true);
    
    String magAxisTitle = "10 * log10( RESP(f) )";
    magAxis = new NumberAxis(magAxisTitle);
    magAxis.setLabelFont(bold);
    magAxis.setAutoRange(true);
    ((NumberAxis) magAxis).setAutoRangeIncludesZero(false);
    
    String phaseAxisTitle = "phi(RESP(f))";
    phaseAxis = new NumberAxis(phaseAxisTitle);
    phaseAxis.setLabelFont(bold);
    phaseAxis.setAutoRange(true);
    ((NumberAxis) phaseAxis).setAutoRangeIncludesZero(false);
    
    plotSelection = new JComboBox<String>();
    plotSelection.addItem("Step function");
    plotSelection.addItem("Response magnitude");
    plotSelection.addItem("Response argument");
    plotSelection.addActionListener(this);
    
    applyAxesToChart();
    
    this.setLayout( new GridBagLayout() );
    GridBagConstraints gbc = new GridBagConstraints();
    gbc.gridx = 0; gbc.gridy = 0;
    gbc.weightx = 1.0; gbc.weighty = 1.0;
    gbc.gridwidth = 3;
    gbc.fill = GridBagConstraints.BOTH;
    gbc.anchor = GridBagConstraints.CENTER;
    this.add(chartPanel, gbc);
    
    // add empty space on left side to space out other components
    JPanel space = new JPanel();
    space.setMaximumSize( plotSelection.getMaximumSize() );
    space.setPreferredSize( plotSelection.getPreferredSize() );
    gbc.weighty = 0.0; gbc.weightx = 1.0;
    gbc.fill = GridBagConstraints.BOTH;
    gbc.gridwidth = 1;
    gbc.gridy += 1;
    this.add(space, gbc);
    
    gbc.fill = GridBagConstraints.NONE;
    gbc.anchor = GridBagConstraints.CENTER;
    gbc.gridx += 1;
    gbc.weightx = 0.0;
    this.add(save, gbc);
    
    gbc.weightx = 1.0;
    gbc.gridx += 1;
    gbc.anchor = GridBagConstraints.LINE_END;
    this.add(plotSelection, gbc);
    
    plotTheseInBold = new String[]{};
    
    
  }
   
  @Override
  public void actionPerformed(ActionEvent e) {
    
    super.actionPerformed(e);
    
    if ( e.getSource() == plotSelection ) {
      if (!set) {    
        applyAxesToChart();
        return;
      }
      
      JFreeChart[] charts = getCharts();
      int idx = plotSelection.getSelectedIndex();
      chart = charts[idx];
      
      chartPanel.setChart(chart);
      chartPanel.restoreAutoBounds();
      chartPanel.validate();
      
      return;   
    }
    
  }
  
  @Override
  public JFreeChart[] getCharts() {
    return new JFreeChart[]{stepChart, magChart, phaseChart};
  }
  
  /**
   * Used to get the text that will populate the inset box for the plots
   * @return String to place in TextTitle
   */
  @Override
  public String getInsetString() {
    
    return getInsetString( (StepExperiment) expResult );
  
  }
  
  @Override
  public String getMetadataString() {
    StepExperiment stex = (StepExperiment) expResult;
    StringBuilder sb = new StringBuilder();
    sb.append("Residuals:\n");
    double[] resids = stex.getResiduals();
    sb.append("Initial:  ");
    sb.append(resids[0]);
    sb.append('\n');
    sb.append("Fit:  ");
    sb.append(resids[1]);
    sb.append('\n');
    sb.append( super.getMetadataString() );
    return sb.toString();
  }
  
  @Override
  public ValueAxis getXAxis() {
    
    if ( null == plotSelection ) {
      return xAxis;
    }
    
    ValueAxis[] array = new ValueAxis[]{xAxis, freqAxis, freqAxis};
    int idx = plotSelection.getSelectedIndex();
    return array[idx];
  }
  
  @Override
  public ValueAxis getYAxis() {
    
    if ( null == plotSelection ) {
      return yAxis;
    }
    
    ValueAxis[] array = new ValueAxis[]{yAxis, magAxis, phaseAxis};
    int idx = plotSelection.getSelectedIndex();
    return array[idx];
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

  @Override
  public int panelsNeeded() {
    return 2;
  }

  /**
   * Pass in and retrieve data from the step experiment backend, to plot;
   * this is both the timeseries data as well as a title inset displaying
   * the parameters used in the plot calculations
   */
  @Override
  public void updateData(final DataStore ds) {
    
    set = false;
    
    clearChartAndSetProgressData();
    
    expResult.setData(ds);
    
    set = true;
    
    TextTitle result = new TextTitle();
    result.setText( getInsetString() );
    result.setBackgroundPaint(Color.white);
    
    XYSeriesCollection stepData = expResult.getData().get(0);
    stepChart = buildChart(stepData, xAxis, yAxis);
    System.out.println("test1");
    XYPlot xyp = stepChart.getXYPlot();
    // xyp.setDomainAxis(xAxis);
    // xyp.setRangeAxis(yAxis);
    // xAxis.setRange( stepData.getDomainBounds(true) );
    // yAxis.setRange( stepData.getRangeBounds(true) );
    XYTitleAnnotation xyt = new XYTitleAnnotation(0.98, 0.5, result,
        RectangleAnchor.RIGHT);
    xyp.clearAnnotations();
    xyp.addAnnotation(xyt);
    
    XYSeriesCollection magData = expResult.getData().get(1);
    magChart = buildChart(magData, freqAxis, magAxis);
    System.out.println("test2");
    xyp = magChart.getXYPlot();
    //xyp.setDomainAxis(freqAxis);
    //xyp.setRangeAxis(magAxis);
    xyt = new XYTitleAnnotation(1., 0., result, RectangleAnchor.BOTTOM_RIGHT);
    xyp.clearAnnotations();
    xyp.addAnnotation(xyt);
    
    XYSeriesCollection phaseData = expResult.getData().get(2);
    phaseChart = buildChart(phaseData, freqAxis, phaseAxis);
    //xyp = phaseChart.getXYPlot();
    //xyp.setDomainAxis(freqAxis);
    //xyp.setRangeAxis(phaseAxis);
    
    JFreeChart[] charts = getCharts();
    int idx = plotSelection.getSelectedIndex();
    chart = charts[idx];
    
    chartPanel.setChart(charts[idx]);
    chartPanel.restoreAutoBounds();
    chartPanel.validate();
  }

}
