package asl.sensor.gui;

import java.awt.Color;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.text.DecimalFormat;

import javax.swing.JComboBox;
import javax.swing.event.ChangeEvent;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.annotations.XYTitleAnnotation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.title.TextTitle;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.RectangleAnchor;

import asl.sensor.experiment.ExperimentEnum;
import asl.sensor.experiment.ExperimentFactory;
import asl.sensor.experiment.GainSixExperiment;
import asl.sensor.input.DataStore;

public class GainSixPanel extends GainPanel {
  
  
  /**
   * 
   */
  private static final long serialVersionUID = 6140615094847886109L;

  /**
   * Static helper method for getting the formatted inset string directly
   * from a GainExperiment
   * @param gn GainSixExperiment with data to be extracted
   * @param plotIdx Plot to have this inset applied to 
   * @param refIdx Index of data to be loaded as reference (i.e., 0)
   * @param lowPrd low period boundary to take stats over
   * @param highPrd high period boundary to take stats over
   * @return String with data representation of experiment results (mean, sdev)
   */
  public static String 
  getInsetString(GainSixExperiment gn, int plotIdx, 
      int refIdx, double lowPrd, double highPrd) {

    double[][] meanAndStdDevAll = 
        gn.getStatsFromFreqs(refIdx, 1/lowPrd, 1/highPrd);
    
    double[] meanAndStdDev = meanAndStdDevAll[plotIdx];

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
    
    DecimalFormat df = new DecimalFormat("#.###");
    
    if (plotIdx == 0) {
      sb.append("\nNorth azimuth (deg): ");
      sb.append( df.format( Math.toDegrees( gn.getNorthAzimuth() ) ) );
    } else if (plotIdx == 1) {
      sb.append("\nEast azimuth (rad): ");
      sb.append( df.format( Math.toDegrees( gn.getEastAzimuth() ) ) );
    }
    
    return sb.toString();
  }

  protected JFreeChart northChart, eastChart, vertChart;
  protected JComboBox<String> plotSelection;

  /**
   * Instantiate the panel, including sliders and stat calc button
   * @param exp
   */
  public GainSixPanel(ExperimentEnum exp) {
    // instantiate common components
    super(exp);
    // make sure the experiment is gain-six
    expResult = ExperimentFactory.createExperiment(exp);
    
    for (int i = 0; i < 2; ++i) {
      int num = i + 1;
      channelType[3 * i] = "North sensor " + num + " (RESP required)";
      channelType[(3 * i) + 1] = "East sensor " + num  + " (RESP required)";
      channelType[(3 * i) + 2] = "Vertical sensor " + num + " (RESP required)";
    }
    
    plotTheseInBold = new String[]{"NLNM"};
    
    String xTitle = getXAxis().getLabel();
    String yTitle = getYAxis().getLabel();
    
    northChart = 
        ChartFactory.createXYLineChart( expType.getName() + " (North)",
        xTitle, yTitle, null);
    eastChart = 
        ChartFactory.createXYLineChart( expType.getName() + " (East)",
        xTitle, yTitle, null);
    vertChart = 
        ChartFactory.createXYLineChart( expType.getName() + " (Vertical)",
        xTitle, yTitle, null);
    
    // create layout    
    removeAll();
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
    
    gbc.fill = GridBagConstraints.HORIZONTAL;
    gbc.gridx += 1;
    gbc.weightx = 0;
    gbc.anchor = GridBagConstraints.WEST;
    plotSelection = new JComboBox<String>();
    plotSelection.addItem("North component plot");
    plotSelection.addItem("East component plot");
    plotSelection.addItem("Vertical component plot");
    plotSelection.addActionListener(this);
    add(plotSelection, gbc);
    
  }

  /**
   * Calls functions to do replotting and stat recalculations when different
   * timeseries are selected or the recalculate button is hit
   */
  @Override
  public void actionPerformed(ActionEvent e) {
    
    if ( e.getSource() == recalcButton ) {
      
      // TODO: set title for each chart;
      setTitle();
      
      recalcButton.setEnabled(false);
      
      return;
    } 
    
    super.actionPerformed(e); // saving?
    
    if ( e.getSource() == plotSelection ) {
      int idx = plotSelection.getSelectedIndex();
      
      JFreeChart[] charts = getCharts();
      chart = charts[idx];
      
      chartPanel.setChart(chart);
      chartPanel.setMouseZoomable(true);
      
      return;
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

    GainSixExperiment gn = (GainSixExperiment) expResult;

    // want to default to octave centered at highest value of fixed freq
    freqRange = gn.getOctaveCenteredAtPeak(refIdx);

    // get the locations (x-axis values) of frequency range as intervals
    lowPrd = Math.min(1/freqRange[0], 1/freqRange[1]);
    highPrd = Math.max(1/freqRange[0], 1/freqRange[1]);
    
    double[] minMax = gn.getMinMaxFrequencies();

    // since intervals of incoming data match, so too limits of plot
    // this is used in mapping scale of slider to x-axis values
    low = Math.log10(minMax[0]); // value when slider is 0
    high = Math.log10(minMax[1]); // value when slider is 1000
    leftSliderValue = mapPeriodToSlider(lowPrd);
    rightSliderValue = mapPeriodToSlider(highPrd);
    
    // plot has 3 components: source, destination, NLNM line plot
    JFreeChart[] charts = getCharts();
    for (int i = 0; i < charts.length; ++i) {
      XYSeriesCollection xyscIn = expResult.getData().get(i);
      XYSeriesCollection xysc = new XYSeriesCollection();
      xysc.addSeries( xyscIn.getSeries(refIdx) );
      xysc.addSeries( xyscIn.getSeries(idx1) );
      xysc.addSeries( xyscIn.getSeries("NLNM") );

      charts[i] = buildChart(xysc);
      
      // set vertical bars and enable sliders
     
      leftSlider.setValue(leftSliderValue);
      rightSlider.setValue(rightSliderValue);

      // set the domain to match the boundaries of the octave centered at peak
      charts[i] = setDomainMarkers(lowPrd, highPrd, charts[i]);
      
      // and now set the sliders to match where that window is
      leftSlider.setEnabled(true);
      rightSlider.setEnabled(true);

    }
    
    // make sure pointers to each chart are set properly
    northChart = charts[0];
    eastChart = charts[1];
    vertChart = charts[2];

    // lastly, display the calculated statistics in a textbox in the corner
    setTitle();
    
    int idx = plotSelection.getSelectedIndex();
    chart = charts[idx];

    // obviously, set the chart
    chartPanel.setChart(chart);
    chartPanel.setMouseZoomable(false);
    
  }
  
  @Override
  public JFreeChart[] getCharts() {
    return new JFreeChart[]{northChart, eastChart, vertChart};
  }
  

  /**
   * Since each chart has a unique inset, this method ensures each one has
   * its own unique inset; all inset strings are compiled in the getInsetString
   * method.
   * @param idx Index of experiment sub-index to get data from
   * @return String of results returned by that experiment
   */
  private String getInsetStringPerChart(int idx) {
    int leftPos = leftSlider.getValue();
    double lowPrd = mapSliderToPeriod(leftPos);
    int rightPos = rightSlider.getValue();
    double highPrd = mapSliderToPeriod(rightPos);
    
    // remove old bars and draw the new ones
    // setDomainMarkers(lowPrd, highPrd, xyp);
    
    int refIdx = refSeries.getSelectedIndex();
    
    GainSixExperiment gn = (GainSixExperiment) expResult;
 
    return getInsetString(gn, idx, refIdx, lowPrd, highPrd);
  }
  
  @Override
  public String getInsetStrings() {
    // displays all inset strings as according to report
    StringBuilder sb = new StringBuilder();
    String[] labels = new String[]{"North data:", "East data:", "Vert. data"};
    for (int i = 0; i < 3; ++i) {
      sb.append(labels[i]);
      sb.append('\n');
      sb.append( getInsetStringPerChart(i) );
      sb.append('\n');
    }
    
    return sb.toString();
  }
  
  @Override
  public int panelsNeeded() {
    return 6;
  }
  
  /**
   * Used to populate the comboboxes with the incoming data
   * @param ds DataStore object being processed 
   */
  private void setDataNames(DataStore ds) {
    
    refSeries.setEnabled(false);
    
    refSeries.removeAllItems();
    
    refSeries.addItem("Data from sensor set 1");
    refSeries.addItem("Data from sensor set 2");
    
    refSeries.setSelectedIndex(0);
  }

  /**
   * Displays the statistic results when the calculate button is hit
   * in an inset box on the chart, also used as text in report generation
   */
  private void setTitle() {
    JFreeChart[] charts = getCharts();
    for (int i = 0; i < charts.length; ++i) {
      XYPlot xyp = charts[i].getXYPlot();
      TextTitle result = new TextTitle();
      String temp = getInsetStringPerChart(i);
      result.setText(temp);
      result.setBackgroundPaint(Color.white);
      XYTitleAnnotation xyt = new XYTitleAnnotation(0.98, 0.98, result,
          RectangleAnchor.TOP_RIGHT);
      xyp.clearAnnotations();
      xyp.addAnnotation(xyt);
    }
    
  }
  
  @Override
  public void stateChanged(ChangeEvent e) {
    
    // enforce slider boundaries, same as deriving class
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
      
      // get plots (where we put the vertical bars) 
      // ** THIS IS WHAT IS DIFFERENT **
      for ( JFreeChart chart : getCharts() ) {
        XYPlot xyp = chart.getXYPlot();
        
        // clear out annotations to prevent issues with misleading data
        xyp.clearAnnotations();
        
        // convert slider locations to (log-scale) frequency
        int leftPos = leftSlider.getValue();
        double lowPrd = mapSliderToPeriod(leftPos);
        int rightPos = rightSlider.getValue();
        double highPrd = mapSliderToPeriod(rightPos);
        
        // remove old bars and draw the new ones
        setDomainMarkers(lowPrd, highPrd, chart);
      }
      
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
