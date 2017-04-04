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
import javax.swing.SwingWorker;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

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

  
  /**
   * 
   */
  private static final long serialVersionUID = 6697458429989867529L;
  private JSlider leftSlider;
  private JSlider rightSlider;
  private JComboBox<String> firstSeries;
  private JComboBox<String> secondSeries;
  private JButton recalcButton;
  
  private double low, high;
  /**
   * Max value of slider (ranges from 0 to 1000, converted to log10 scale)
   */
  private static final int SLIDER_MAX = 1000;
  
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
    
    xAxisTitle = "Period (s)";
    yAxisTitle = "Power (rel. 1 (m/s^2)^2/Hz)";
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
    firstSeries = new JComboBox<String>();
    firstSeries.addActionListener(this);
    firstSeries.setEnabled(false);
    secondSeries = new JComboBox<String>();
    secondSeries.addActionListener(this);
    secondSeries.setEnabled(false);
    
    for (int i = 0; i < 2; ++i) {
      String out = "FILE NOT LOADED (" + i + ")";
      firstSeries.addItem(out);
      secondSeries.addItem(out);
    }

    firstSeries.setSelectedIndex(0);
    secondSeries.setSelectedIndex(1);
    
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
    this.add(firstSeries, gbc);
    gbc.weightx = 0;
    gbc.gridx += 1;
    gbc.fill = GridBagConstraints.NONE;
    this.add(save, gbc);
    gbc.gridx += 1;
    gbc.fill = GridBagConstraints.BOTH;
    this.add(secondSeries, gbc);
    
  }

  /**
   * Calls functions to do replotting and stat recalculations when different
   * timeseries are selected or the recalculate button is hit
   */
  @Override
  public void actionPerformed(ActionEvent e) {
    super.actionPerformed(e); // saving?
    
    int idx0 = firstSeries.getSelectedIndex();
    int idx1 = secondSeries.getSelectedIndex();
    
    if ( e.getSource() == recalcButton ) {
      
      // we need to get the values of the sliders again, convert to frequency
      int leftPos = leftSlider.getValue();
      double lowPrd = mapSliderToPeriod(leftPos);
      int rightPos = rightSlider.getValue();
      double highPrd = mapSliderToPeriod(rightPos);
      
      double[] meanAndStdDev = 
          ((GainExperiment) expResult).getStatsFromFreqs(
              idx0, idx1, 1/lowPrd, 1/highPrd);
      
      double mean = meanAndStdDev[0];
      double sDev = meanAndStdDev[1];
      double refGain = meanAndStdDev[2];
      double calcGain = meanAndStdDev[3];
      
      setTitle(mean, sDev, refGain, calcGain);
      
      recalcButton.setEnabled(false);
      
      return;
    } 
    if ( e.getSource() == firstSeries ) {
      
      // if we got here from removing the items from the list
      // (which happens when we load in new data)
      // don't do anything
      if ( !firstSeries.isEnabled() ){
        return;
      }
      
      // don't allow the same value for the two indices (plot behaves badly)
      // assume the user is setting firstSeries selection correctly,
      // and thus make sure that the secondSeries doesn't have a collision
      // we disable the combo box to prevent a second event from triggering
      if (idx0 == idx1) {
        secondSeries.setEnabled(false);
        idx1 = (idx1 + 1) % secondSeries.getItemCount();
        secondSeries.setSelectedIndex(idx1);
        secondSeries.setEnabled(true);
      }
      
    } else if ( e.getSource() == secondSeries ) {      
      
      // same as above, do nothing
      if ( !secondSeries.isEnabled() ) {
        return;
      }
      
      // same as with the above, but assume secondSeries selection intentional
      if (idx0 == idx1) {
        firstSeries.setEnabled(false);
        idx0 = (idx0 + 1) % firstSeries.getItemCount();
        firstSeries.setSelectedIndex(idx0);
        firstSeries.setEnabled(true);
      }
      
    }
    // now that we have a guarantee of no collision, update data accordingly    
    if ( e.getSource() == firstSeries || e.getSource() == secondSeries) {
      // if we selected a new series to plot, redraw the chart
      updateDataDriver(idx0, idx1);
    }
    
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
  
  /**
   * Used to populate the comboboxes with the incoming data
   * @param ds DataStore object being processed 
   */
  private void setDataNames(DataStore ds) {
    
    firstSeries.setEnabled(false);
    secondSeries.setEnabled(false);
    
    firstSeries.removeAllItems();
    secondSeries.removeAllItems();
    
    Set<String> preventDuplicates = new HashSet<String>();
    
    for (int i = 0; i < 2; ++i) {
      String name = ds.getBlock(i).getName();
      while ( preventDuplicates.contains(name) ) {
        name += "_";
      }
      preventDuplicates.add(name);
      firstSeries.addItem(name);
      secondSeries.addItem(name);
    }
    
    firstSeries.setSelectedIndex(0);
    secondSeries.setSelectedIndex(1);
  }
  
    
  /**
   * Draws the lines marking the boundaries of the current window
   * @param lowPrd lower x-axis value (period, in seconds)
   * @param highPrd upper x-axis value (period, in seconds)
   * @param xyp plot displayed in this object's chart
   */
  private void setDomainMarkers(double lowPrd, double highPrd, XYPlot xyp) {
    xyp.clearDomainMarkers();
    Marker startMarker = new ValueMarker( lowPrd );
    startMarker.setStroke( new BasicStroke( (float) 1.5 ) );
    Marker endMarker = new ValueMarker( highPrd );
    endMarker.setStroke( new BasicStroke( (float) 1.5 ) );
    xyp.addDomainMarker(startMarker);
    xyp.addDomainMarker(endMarker);
  }

  /**
   * Displays the statistic results when the calculate button is hit
   * in an inset box on the chart
   * @param mean Calculated mean value
   * @param sDev Calculated standard deviation value
   */
  private void 
  setTitle(double mean, double sDev, double refGain, double calcGain) {
    XYPlot xyp = (XYPlot) chartPanel.getChart().getPlot();
    TextTitle result = new TextTitle();
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
    String temp = sb.toString();
    result.setText(temp);
    result.setBackgroundPaint(Color.white);
    XYTitleAnnotation xyt = new XYTitleAnnotation(0.98, 0.98, result,
        RectangleAnchor.TOP_RIGHT);
    xyp.clearAnnotations();
    xyp.addAnnotation(xyt);
  }
  
  /**
   * Called when one of the sliders has been adjusted, used to
   * 
   */
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
      XYPlot xyp = chartPanel.getChart().getXYPlot();
      
      // clear out annotations to prevent issues with misleading data
      xyp.clearAnnotations();
      
      // convert slider locations to (log-scale) frequency
      int leftPos = leftSlider.getValue();
      double lowPrd = mapSliderToPeriod(leftPos);
      int rightPos = rightSlider.getValue();
      double highPrd = mapSliderToPeriod(rightPos);
      
      // remove old bars and draw the new ones
      setDomainMarkers(lowPrd, highPrd, xyp);
    }
  }
  
  /**
   * Sets the data and then calls functions to produce backend calculations
   * when data is loaded in
   */
  @Override
  public void updateData(final DataStore ds) {
    
    setDataNames(ds);
    
    // TODO: this is going to require a refactoring, I know already
    
    SwingWorker<Integer, Void> worker = new SwingWorker<Integer, Void>() {
      
      boolean errorTriggered;
      
      public Integer doInBackground() {
        
        try{
          expResult.setData(ds);
        } catch (IndexOutOfBoundsException e) {
          errorTriggered = true;
          return 1; // TODO: use return value to get error occurent
        }
        
        // need to have 2 series for relative gain
        firstSeries.setEnabled(true);
        secondSeries.setEnabled(true);
        int idx0 = firstSeries.getSelectedIndex();
        int idx1 = secondSeries.getSelectedIndex();
        updateDataDriver(idx0, idx1);
        
        return 0;
      }
      
      public void done() {
        
        if (errorTriggered) {
          displayErrorMessage("INSUFFICIENT DATA LOADED");
          return;
        }
        
      }
      
    };
    
    new Thread(worker).run();

  }

  /**
   * Given input data (including time series collection), get only the relevant
   * ones to display based on combo boxes and then do the statistics on those.
   * Because the range of the sliders is not necessarily set on switch,
   * the statistics are calculated over the octave centered at the plotted 
   * peak value's frequency. This function is cal0led when new data is fed in 
   * or when the combo box active entries change
   */
  private void updateDataDriver(int index0, int index1) {
    
    final int idx0 = index0;
    final int idx1 = index1;
    
    displayInfoMessage("Calculating data...");
    
    SwingWorker<Integer, Void> worker = new SwingWorker<Integer, Void>() {
      
      double lowPrd, highPrd;
      double[] freqRange;
      double[] gainStatistics;
      
      int leftSliderValue, rightSliderValue;
      XYSeriesCollection xysc;
      
      public Integer doInBackground() {

        // plot has 3 components: source, destination, NLNM line plot
        XYSeriesCollection xyscIn = expResult.getData().get(0);
        xysc = new XYSeriesCollection();
        xysc.addSeries( xyscIn.getSeries(idx0) );
        xysc.addSeries( xyscIn.getSeries(idx1) );
        xysc.addSeries( xyscIn.getSeries("NLNM") );
        
        XYSeries xys = xysc.getSeries(0);
        if ( xysc.getSeriesKey(0).equals("NLNM") ) {
          xys = xysc.getSeries(1);
        }
        
        GainExperiment gn = (GainExperiment) expResult;

        // want to default to octave centered at highest value of fixed freq
        freqRange = gn.getOctaveCenteredAtPeak(idx0);
        
        // get the locations (x-axis values) of frequency range as intervals
        lowPrd = Math.min(1/freqRange[0], 1/freqRange[1]);
        highPrd = Math.max(1/freqRange[0], 1/freqRange[1]);
        
        // since intervals of incoming data match, so too limits of plot
        // this is used in mapping scale of slider to x-axis values
        low = Math.log10( xys.getMinX() ); // value when slider is 0
        high = Math.log10( xys.getMaxX() ); // value when slider is 1000
        leftSliderValue = mapPeriodToSlider(lowPrd);
        rightSliderValue = mapPeriodToSlider(highPrd);
        
        
        // get statistics from frequency (convert from period)
        gainStatistics = gn.getStatsFromFreqs(idx0, idx1, 
            freqRange[0], freqRange[1]);
        
        return 0;
      }
      
      public void done() {
        
        displayInfoMessage("Data loaded...drawing chart");
        
        setChart(xysc);
        
        // obviously, set the chart
        chartPanel.setChart(chart);
        chartPanel.setMouseZoomable(false);
        
        // set vertical bars and enable sliders
        XYPlot xyp = chartPanel.getChart().getXYPlot();

        leftSlider.setValue(leftSliderValue);
        rightSlider.setValue(rightSliderValue);
        
        // set the domain to match the boundaries of the octave centered at peak
        setDomainMarkers(lowPrd, highPrd, xyp);
        
        // and now set the sliders to match where that window is
        leftSlider.setEnabled(true);
        rightSlider.setEnabled(true);
        
        // lastly, display the calculated statistics in a textbox in the corner
        double mean = gainStatistics[0];
        double sDev = gainStatistics[1];
        double refGain = gainStatistics[2];
        double calcGain = gainStatistics[3];
        
        setTitle(mean, sDev, refGain, calcGain);
      }
      
    };
    
    worker.execute();

  }

  @Override
  public int panelsNeeded() {
    return 2;
  }
  
}
