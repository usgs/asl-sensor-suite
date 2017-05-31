package asl.sensor.gui;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.util.HashSet;
import java.util.Set;

import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JSlider;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.jfree.chart.JFreeChart;
import org.jfree.chart.annotations.XYTitleAnnotation;
import org.jfree.chart.axis.LogarithmicAxis;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.Marker;
import org.jfree.chart.plot.ValueMarker;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.title.TextTitle;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.RectangleAnchor;

import asl.sensor.experiment.ExperimentEnum;
import asl.sensor.experiment.GainExperiment;
import asl.sensor.input.DataStore;

/**
 * Panel to display the results of the gain experiment calculations.
 * In addition to the expected components of the experiment panel,
 * this also includes selectors to set reference and calculated gain timeseries
 * targets, and sliders to set the window range over which to calculate
 * the gain statistics
 * @author akearns
 *
 */
public class GainPanel extends ExperimentPanel 
implements ChangeListener {

  
  private static final long serialVersionUID = 6697458429989867529L;
  
  /**
   * Max value of slider (ranges from 0 to 1000, converted to log10 scale)
   */
  public static final int SLIDER_MAX = 1000;
  
  /**
   * Static helper method for getting the formatted inset string directly
   * from a GainExperiment
   * @param gn GainExperiment with data to be extracted
   * @param refIdx Index of data to be loaded as reference (i.e., 0)
   * @param lowPrd low period boundary to take stats over
   * @param highPrd high period boundary to take stats over
   * @return String with data representation of experiment results (mean, sdev)
   */
  public static String 
  getInsetString(GainExperiment gn, int refIdx, double lowPrd, double highPrd) {

    double[] meanAndStdDev = 
        gn.getStatsFromFreqs(refIdx, 1/lowPrd, 1/highPrd);

    double mean = meanAndStdDev[0];
    double sDev = meanAndStdDev[1];
    double refGain = meanAndStdDev[2];
    double calcGain = meanAndStdDev[3];
    
    StringBuilder sb = new StringBuilder();
    sb.append("ratio: ");
    sb.append(mean);
    sb.append("\n");
    sb.append("sigma: ");
    sb.append(sDev);
    sb.append("\n");
    sb.append("ref. gain: ");
    sb.append(refGain);
    sb.append("\n");
    sb.append("** CALCULATED GAIN: ");
    sb.append(calcGain);
    return sb.toString();
  }
  protected JSlider leftSlider;
  protected JSlider rightSlider;
  protected JComboBox<String> refSeries;
  protected JButton recalcButton;
  
  protected double low, high;

  /**
   * Instantiate the panel, including sliders and stat calc button
   * @param exp
   */
  public GainPanel(ExperimentEnum exp) {
    // instantiate common components
    super(exp);
    
    for (int i = 0; i < 2; ++i) {
      channelType[i] = "Input data (RESP required)";
    }
    
    plotTheseInBold = new String[]{"NLNM"};
    
    String xAxisTitle = "Period (s)";
    String yAxisTitle = "Power (rel. 1 (m/s^2)^2/Hz)";
    xAxis = new LogarithmicAxis(xAxisTitle);
    yAxis = new NumberAxis(yAxisTitle);
    yAxis.setAutoRange(true);
    ( (NumberAxis) yAxis).setAutoRangeIncludesZero(false);
    Font bold = xAxis.getLabelFont().deriveFont(Font.BOLD);
    xAxis.setLabelFont(bold);
    yAxis.setLabelFont(bold);
    
    applyAxesToChart();
    
    // instantiate unique components
    leftSlider = new JSlider(0, SLIDER_MAX, 0);
    leftSlider.addChangeListener(this);
    leftSlider.setEnabled(false);
    rightSlider = new JSlider(0, SLIDER_MAX, SLIDER_MAX);
    rightSlider.addChangeListener(this);
    rightSlider.setEnabled(false);
    
    recalcButton = new JButton("Recalc over range");
    recalcButton.setEnabled(false);
    recalcButton.addActionListener(this);
    
    // add dummy entries to the combo box, but don't let them get filled
    refSeries = new JComboBox<String>();
    refSeries.addActionListener(this);
    refSeries.setEnabled(false);
    
    for (int i = 0; i < 2; ++i) {
      String out = "FILE NOT LOADED (" + i + ")";
      refSeries.addItem(out);
    }

    refSeries.setSelectedIndex(0);
    
    // create layout    
    this.setLayout( new GridBagLayout() );
    GridBagConstraints gbc = new GridBagConstraints();
    
    gbc.fill = GridBagConstraints.BOTH;
    gbc.gridx = 0; gbc.gridy = 0;
    gbc.weightx = 1.0; gbc.weighty = 1.0;
    gbc.gridwidth = 3;
    gbc.anchor = GridBagConstraints.CENTER;
    this.add(chartPanel, gbc);
    
    gbc.gridx = 0; gbc.gridy += 1;
    gbc.weighty = 0;
    gbc.gridwidth = 1;
    gbc.anchor = GridBagConstraints.EAST;
    this.add(leftSlider, gbc);
    gbc.fill = GridBagConstraints.NONE;
    gbc.gridx += 1;
    gbc.anchor = GridBagConstraints.CENTER;
    gbc.weightx = 0;
    this.add(recalcButton, gbc);
    gbc.fill = GridBagConstraints.BOTH;
    gbc.gridx += 1;
    gbc.anchor = GridBagConstraints.WEST;
    gbc.weightx = 1;
    this.add(rightSlider, gbc);
    
    gbc.gridx = 0; gbc.gridy += 1;
    gbc.fill = GridBagConstraints.BOTH;
    gbc.anchor = GridBagConstraints.CENTER;
    this.add(refSeries, gbc);
    gbc.weightx = 0;
    gbc.gridx += 1;
    gbc.fill = GridBagConstraints.NONE;
    this.add(save, gbc);
    
  }

  /**
   * Calls functions to do replotting and stat recalculations when different
   * timeseries are selected or the recalculate button is hit
   */
  @Override
  public void actionPerformed(ActionEvent e) {
    super.actionPerformed(e); // saving?
    
    if ( e.getSource() == recalcButton ) {
      
      setTitle();
      
      recalcButton.setEnabled(false);
      
      return;
    } 
    if ( e.getSource() == refSeries ) {
      
      // if we got here from removing the items from the list
      // (which happens when we load in new data)
      // don't do anything
      if ( !refSeries.isEnabled() ){
        return;
      }

      // if we selected a new series to plot, redraw the chart
      drawCharts();
      
    }
    
  }
  
  @Override
  /**
   * Given input data (including time series collection), get only the relevant
   * ones to display based on combo boxes and then do the statistics on those.
   * Because the range of the sliders is not necessarily set on switch,
   * the statistics are calculated over the octave centered at the plotted 
   * peak value's frequency. This function is called when new data is fed in 
   * or when the combo box active entries change
   */
  protected void drawCharts() {
    
    final int refIdx = refSeries.getSelectedIndex();
    final int idx1 = (refIdx + 1) % 2;


    double lowPrd, highPrd;
    double[] freqRange;

    int leftSliderValue, rightSliderValue;
    XYSeriesCollection xysc;

    // plot has 3 components: source, destination, NLNM line plot
    XYSeriesCollection xyscIn = expResult.getData().get(0);
    xysc = new XYSeriesCollection();
    xysc.addSeries( xyscIn.getSeries(refIdx) );
    xysc.addSeries( xyscIn.getSeries(idx1) );
    xysc.addSeries( xyscIn.getSeries("NLNM") );

    XYSeries xys = xysc.getSeries(0);
    if ( xysc.getSeriesKey(0).equals("NLNM") ) {
      xys = xysc.getSeries(1);
    }

    GainExperiment gn = (GainExperiment) expResult;

    // want to default to octave centered at highest value of fixed freq
    freqRange = gn.getOctaveCenteredAtPeak(refIdx);

    // get the locations (x-axis values) of frequency range as intervals
    lowPrd = Math.min(1/freqRange[0], 1/freqRange[1]);
    highPrd = Math.max(1/freqRange[0], 1/freqRange[1]);

    // since intervals of incoming data match, so too limits of plot
    // this is used in mapping scale of slider to x-axis values
    low = Math.log10( xys.getMinX() ); // value when slider is 0
    high = Math.log10( xys.getMaxX() ); // value when slider is 1000
    leftSliderValue = mapPeriodToSlider(lowPrd);
    rightSliderValue = mapPeriodToSlider(highPrd);

    setChart(xysc);

    // obviously, set the chart
    chartPanel.setChart(chart);
    chartPanel.setMouseZoomable(false);
    
    leftSlider.setValue(leftSliderValue);
    rightSlider.setValue(rightSliderValue);

    // set the domain to match the boundaries of the octave centered at peak
    chartPanel.setChart( setDomainMarkers(lowPrd, highPrd, chart) );

    // and now set the sliders to match where that window is
    leftSlider.setEnabled(true);
    rightSlider.setEnabled(true);

    // lastly, display the calculated statistics in a textbox in the corner
    setTitle();

  }
  
  @Override
  public String getInsetStrings() {
    int leftPos = leftSlider.getValue();
    double lowPrd = mapSliderToPeriod(leftPos);
    int rightPos = rightSlider.getValue();
    double highPrd = mapSliderToPeriod(rightPos);
    
    // remove old bars and draw the new ones
    // setDomainMarkers(lowPrd, highPrd, xyp);
    
    int refIdx = refSeries.getSelectedIndex();
    
    GainExperiment gn = (GainExperiment) expResult;
 
    return getInsetString(gn, refIdx, lowPrd, highPrd);
  }
  
  @Override
  public String getMetadataString() {

    // get range of data over which the gain statistics were calculated
    int leftPos = leftSlider.getValue();
    double lowPrd = mapSliderToPeriod(leftPos);
    int rightPos = rightSlider.getValue();
    double highPrd = mapSliderToPeriod(rightPos);
    
    StringBuilder sb = new StringBuilder();
    
    sb.append("Range used in stat calculation: ");
    sb.append(lowPrd);
    sb.append(" to ");
    sb.append(highPrd);
    
    sb.append( super.getMetadataString() );
    
    return sb.toString();
  }
  
    
  /**
   * Converts x-axis value from log scale to linear, to get slider position
   * @param prd period value marking data window boundary
   * @return value of slider (ranges from 0 to SLIDER_MAX)
   */
  public int mapPeriodToSlider(double prd) {
    double scale = (high - low)/SLIDER_MAX; // recall slider range is 0 to 1000
    return (int) ( ( Math.log10(prd) - low ) / scale );
  }

  /**
   * Converts the slider position to a logarithmic scale matching x-axis values
   * which is the period given in a rate of seconds
   * @param position value of slider
   * @return x-axis value corresponding to that position
   */
  public double mapSliderToPeriod(int position) {
    double scale = (high - low)/SLIDER_MAX; // slider range is 0 to 1000
    return Math.pow(10, low + (scale * position) );
  }
  

  @Override
  public int panelsNeeded() {
    return 2;
  }
  
  /**
   * Used to populate the comboboxes with the incoming data
   * @param ds DataStore object being processed 
   */
  private void setDataNames(DataStore ds) {
    
    refSeries.setEnabled(false);
    
    refSeries.removeAllItems();
    
    Set<String> preventDuplicates = new HashSet<String>();
    
    for (int i = 0; i < 2; ++i) {
      String name = ds.getBlock(i).getName();
      while ( preventDuplicates.contains(name) ) {
        name += "_";
      }
      preventDuplicates.add(name);
      refSeries.addItem(name);
    }
    
    refSeries.setSelectedIndex(0);
  }

  /**
   * Draws the lines marking the boundaries of the current window
   * @param lowPrd lower x-axis value (period, in seconds)
   * @param highPrd upper x-axis value (period, in seconds)
   * @param xyp plot displayed in this object's chart
   * @return XYPlot XYPlot with new domain markers set
   */
  protected static JFreeChart 
  setDomainMarkers(double lowPrd, double highPrd, JFreeChart chart) {
    XYPlot xyp = chart.getXYPlot();
    xyp.clearDomainMarkers();
    Marker startMarker = new ValueMarker(lowPrd);
    startMarker.setStroke( new BasicStroke( (float) 1.5 ) );
    Marker endMarker = new ValueMarker(highPrd);
    endMarker.setStroke( new BasicStroke( (float) 1.5 ) );
    xyp.addDomainMarker(startMarker);
    xyp.addDomainMarker(endMarker);
    return chart;
  }

  /**
   * Displays the statistic results when the calculate button is hit
   * in an inset box on the chart, also used as text in report generation
   */
  private void setTitle() {
    XYPlot xyp = (XYPlot) chartPanel.getChart().getPlot();
    TextTitle result = new TextTitle();
    String temp = getInsetStrings();
    result.setText(temp);
    result.setBackgroundPaint(Color.white);
    XYTitleAnnotation xyt = new XYTitleAnnotation(0.98, 0.98, result,
        RectangleAnchor.TOP_RIGHT);
    xyp.clearAnnotations();
    xyp.addAnnotation(xyt);
  }
  
  @Override
  public void stateChanged(ChangeEvent e) {
    
    // enforce slider boundaries
    if ( e.getSource() == leftSlider ) {
      if ( leftSlider.getValue() > rightSlider.getValue() - 10 ) {
        leftSlider.setValue( rightSlider.getValue() - 10 );
        if ( leftSlider.getValue() < 0 ) {
          leftSlider.setValue(0);
          rightSlider.setValue(10);
        }
      }
    } else if ( e.getSource() == rightSlider ) {
      if ( leftSlider.getValue() + 10 > rightSlider.getValue() ) {
        rightSlider.setValue( leftSlider.getValue() + 10 );
        if ( rightSlider.getValue() > SLIDER_MAX ) {
          rightSlider.setValue(SLIDER_MAX);
          leftSlider.setValue(SLIDER_MAX - 10);
        }
      }
    }
    
    if (e.getSource() == leftSlider || e.getSource() == rightSlider) {
      
      // now we need to redraw the vertical bars
      
      // new slider window means new results on calculation
      recalcButton.setEnabled(true);
      
      // get plot (where we put the vertical bars)
      JFreeChart chart = chartPanel.getChart();
      
      // clear out annotations to prevent issues with misleading data
      chart.getXYPlot().clearAnnotations();
      
      // convert slider locations to (log-scale) frequency
      int leftPos = leftSlider.getValue();
      double lowPrd = mapSliderToPeriod(leftPos);
      int rightPos = rightSlider.getValue();
      double highPrd = mapSliderToPeriod(rightPos);
      
      // remove old bars and draw the new ones
      chartPanel.setChart( setDomainMarkers(lowPrd, highPrd, chart) );
    }
  }
  
  @Override
  protected void updateData(final DataStore ds) {
    
    set = true;
    
    setDataNames(ds);
    
    expResult.runExperimentOnData(ds);

    // need to have 2 series for relative gain
    refSeries.setEnabled(true);

  }
  
}
