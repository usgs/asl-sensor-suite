package asl.sensor.gui;

import asl.sensor.ExperimentFactory;
import asl.sensor.experiment.GainExperiment;
import asl.sensor.input.DataStore;
import java.awt.BasicStroke;
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

/**
 * Panel to display the results of the gain experiment calculations.
 * In addition to the expected components of the experiment panel,
 * this also includes selectors to set reference and calculated gain timeseries
 * targets, and sliders to set the window range over which to calculate
 * the gain statistics
 *
 * @author akearns - KBRWyle
 */
public class GainPanel extends ExperimentPanel
    implements ChangeListener {


  private static final long serialVersionUID = 6697458429989867529L;
  /**
   * Max value of slider (ranges from 0 to 1000, converted to log10 scale)
   */
  private static final int SLIDER_MAX = 1000;
  final JSlider leftSlider;
  final JSlider rightSlider;
  final JComboBox<String> referenceSeries;
  final JButton recalcButton;
  double lowPeriod, highPeriod;

  /**
   * Instantiate the panel, including sliders and stat calc button
   */
  public GainPanel(ExperimentFactory experiment) {
    // instantiate common components
    super(experiment);

    for (int i = 0; i < 2; ++i) {
      channelType[i] = "Input data (RESP required)";
    }

    plotTheseInBold = new String[]{"NLNM"};

    String xAxisTitle = "Period (s)";
    String yAxisTitle = "Power (rel. 1 (m/s^2)^2/Hz)";
    xAxis = new LogarithmicAxis(xAxisTitle);
    yAxis = new NumberAxis(yAxisTitle);
    yAxis.setAutoRange(true);
    ((NumberAxis) yAxis).setAutoRangeIncludesZero(false);
    Font temp = xAxis.getLabelFont();
    Font bold = temp.deriveFont(Font.BOLD, temp.getSize() + 2);
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
    referenceSeries = new JComboBox<>();
    referenceSeries.addActionListener(this);
    referenceSeries.setEnabled(false);

    for (int i = 0; i < 2; ++i) {
      String out = "FILE NOT LOADED (" + i + ")";
      referenceSeries.addItem(out);
    }

    referenceSeries.setSelectedIndex(0);

    // create layout
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

    constraints.gridx = 0;
    constraints.gridy += 1;
    constraints.weighty = 0;
    constraints.gridwidth = 1;
    constraints.anchor = GridBagConstraints.EAST;
    this.add(leftSlider, constraints);
    constraints.fill = GridBagConstraints.NONE;
    constraints.gridx += 1;
    constraints.anchor = GridBagConstraints.CENTER;
    constraints.weightx = 0;
    this.add(recalcButton, constraints);
    constraints.fill = GridBagConstraints.BOTH;
    constraints.gridx += 1;
    constraints.anchor = GridBagConstraints.WEST;
    constraints.weightx = 1;
    this.add(rightSlider, constraints);

    constraints.gridx = 0;
    constraints.gridy += 1;
    constraints.fill = GridBagConstraints.BOTH;
    constraints.anchor = GridBagConstraints.CENTER;
    this.add(referenceSeries, constraints);
    constraints.weightx = 0;
    constraints.gridx += 1;
    constraints.fill = GridBagConstraints.NONE;
    this.add(save, constraints);

  }

  /**
   * Draws the lines marking the boundaries of the current window
   *
   * @param lowPrd lower x-axis value (period, in seconds)
   * @param highPrd upper x-axis value (period, in seconds)
   * @param chart Chart to add domain markers to
   * @return XYPlot XYPlot with new domain markers set
   */
  static JFreeChart setDomainMarkers(double lowPrd, double highPrd, JFreeChart chart) {

    // make sure lower value is assigned to left, higher to right
    double tempLow = Math.min(lowPrd, highPrd);
    highPrd = Math.max(lowPrd, highPrd);
    lowPrd = tempLow;

    XYPlot plot = chart.getXYPlot();
    plot.clearDomainMarkers();
    Marker startMarker = new ValueMarker(lowPrd);
    startMarker.setStroke(new BasicStroke((float) 1.5));
    Marker endMarker = new ValueMarker(highPrd);
    endMarker.setStroke(new BasicStroke((float) 1.5));
    plot.addDomainMarker(startMarker);
    plot.addDomainMarker(endMarker);
    return chart;
  }

  /**
   * Calls functions to do replotting and stat recalculations when different
   * timeseries are selected or the recalculate button is hit
   */
  @Override
  public void actionPerformed(ActionEvent event) {
    super.actionPerformed(event); // saving?

    if (event.getSource() == recalcButton) {
      setTitle();
      recalcButton.setEnabled(false);
      return;
    }
    if (event.getSource() == referenceSeries) {
      // if we got here from removing the items from the list
      // (which happens when we load in new data)
      // don't do anything
      if (!referenceSeries.isEnabled()) {
        return;
      }
      // if we selected a new series to plot, redraw the chart
      drawCharts();
    }
  }


  /**
   * Given input data (including time series collection), get only the relevant
   * ones to display based on combo boxes and then do the statistics on those.
   * Because the range of the sliders is not necessarily set on switch,
   * the statistics are calculated over the octave centered at the plotted
   * peak value's frequency. This function is called when new data is fed in
   * or when the combo box active entries change
   */
  @Override
  protected void drawCharts() {

    final int referenceIndex = referenceSeries.getSelectedIndex();
    final int index1 = (referenceIndex + 1) % 2;

    int leftSliderValue, rightSliderValue;
    XYSeriesCollection timeSeries;

    // plot has 3 components: source, destination, NLNM line plot
    XYSeriesCollection timeSeriesIn = expResult.getData().get(0);
    timeSeries = new XYSeriesCollection();
    timeSeries.addSeries(timeSeriesIn.getSeries(referenceIndex));
    timeSeries.addSeries(timeSeriesIn.getSeries(index1));
    timeSeries.addSeries(timeSeriesIn.getSeries("NLNM"));

    XYSeries xys = timeSeries.getSeries(0);
    if (timeSeries.getSeriesKey(0).equals("NLNM")) {
      xys = timeSeries.getSeries(1);
    }

    GainExperiment gain = (GainExperiment) expResult;
    // since intervals of incoming data match, so too limits of plot
    // this is used in mapping scale of slider to x-axis values
    lowPeriod = Math.log10(xys.getMinX()); // value when slider is 0
    highPeriod = Math.log10(xys.getMaxX()); // value when slider is 1000
    leftSliderValue = mapPeriodToSlider(GainExperiment.DEFAULT_LOW_BOUND);
    rightSliderValue = mapPeriodToSlider(GainExperiment.DEFAULT_UP_BOUND);
    updateReference(referenceIndex);
    gain.setRangeForStatistics(GainExperiment.DEFAULT_LOW_BOUND, GainExperiment.DEFAULT_UP_BOUND);

    setChart(timeSeries);

    chartPanel.setChart(chart);
    chartPanel.setMouseZoomable(false);

    setSliderValues(leftSliderValue, rightSliderValue);

    // set the domain to match the boundaries of the octave centered at peak
    chartPanel.setChart(setDomainMarkers(GainExperiment.DEFAULT_LOW_BOUND,
        GainExperiment.DEFAULT_UP_BOUND, chart));

    // and now set the sliders to match where that window is
    leftSlider.setEnabled(true);
    rightSlider.setEnabled(true);

    // lastly, display the calculated statistics in a textbox in the corner
    setTitle();

  }

  @Override
  public String getMetadataString() {
    // get range of data over which the gain statistics were calculated
    double lowPeriod = mapSliderToPeriod(leftSlider.getValue());
    double highPeriod = mapSliderToPeriod(rightSlider.getValue());

    return "Range used in stat calculation: " + lowPeriod + " to " + highPeriod + '\n'
        + super.getMetadataString();
  }


  /**
   * Converts x-axis value from log scale to linear, to get slider position
   *
   * @param period period value marking data window boundary
   * @return value of slider (ranges from 0 to SLIDER_MAX)
   */
  int mapPeriodToSlider(double period) {
    double scale = (highPeriod - lowPeriod) / SLIDER_MAX; // recall slider range is 0 to 1000
    return (int) ((Math.log10(period) - lowPeriod) / scale);
  }

  /**
   * Converts the slider position to a logarithmic scale matching x-axis values
   * which is the period given in a rate of seconds
   *
   * @param position value of slider
   * @return x-axis value corresponding to that position
   */
  double mapSliderToPeriod(int position) {
    double scale = (highPeriod - lowPeriod) / SLIDER_MAX; // slider range is 0 to 1000
    return Math.pow(10, lowPeriod + (scale * position));
  }


  @Override
  public int panelsNeeded() {
    return 2;
  }

  /**
   * Used to populate the comboboxes with the incoming data
   *
   * @param dataStore DataStore object being processed
   */
  private void setDataNames(DataStore dataStore) {
    referenceSeries.setEnabled(false);
    referenceSeries.removeAllItems();
    Set<String> preventDuplicates = new HashSet<>();
    for (int i = 0; i < 2; ++i) {
      String name = dataStore.getBlock(i).getName();
      while (preventDuplicates.contains(name)) {
        name += "_";
      }
      preventDuplicates.add(name);
      referenceSeries.addItem(name);
    }
    referenceSeries.setSelectedIndex(0);
  }

  void setSliderValues(int leftSliderValue, int rightSliderValue) {
    // enforce constraint that left slider value is the smaller one
    int leftTemp = Math.min(leftSliderValue, rightSliderValue);
    rightSliderValue = Math.max(leftSliderValue, rightSliderValue);
    leftSliderValue = leftTemp;

    leftSlider.setValue(leftSliderValue);
    rightSlider.setValue(rightSliderValue);
  }

  /**
   * Displays the statistic results when the calculate button is hit
   * in an inset box on the chart, also used as text in report generation
   */
  private void setTitle() {
    XYPlot plot = (XYPlot) chartPanel.getChart().getPlot();
    TextTitle result = getDefaultTextTitle();
    result.setText(expResult.getInsetStrings()[0]);
    XYTitleAnnotation xyt = new XYTitleAnnotation(0.98, 0.98, result,
        RectangleAnchor.TOP_RIGHT);
    plot.clearAnnotations();
    plot.addAnnotation(xyt);
  }

  @Override
  public void stateChanged(ChangeEvent event) {

    super.stateChanged(event);

    // enforce slider boundaries
    if (event.getSource() == leftSlider) {
      if (leftSlider.getValue() > rightSlider.getValue() - 10) {
        leftSlider.setValue(rightSlider.getValue() - 10);
        if (leftSlider.getValue() < 0) {
          leftSlider.setValue(0);
          rightSlider.setValue(10);
        }
      }
    } else if (event.getSource() == rightSlider) {
      if (leftSlider.getValue() + 10 > rightSlider.getValue()) {
        rightSlider.setValue(leftSlider.getValue() + 10);
        if (rightSlider.getValue() > SLIDER_MAX) {
          rightSlider.setValue(SLIDER_MAX);
          leftSlider.setValue(SLIDER_MAX - 10);
        }
      }
    }

    if (event.getSource() == leftSlider || event.getSource() == rightSlider) {
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
      chartPanel.setChart(setDomainMarkers(lowPrd, highPrd, chart));
      updateStatistics(lowPrd, highPrd);
    }
  }

  @Override
  protected void updateData(final DataStore dataStore) {
    set = true;
    setDataNames(dataStore);
    expResult.runExperimentOnData(dataStore);
    // need to have 2 series for relative gain
    referenceSeries.setEnabled(true);
  }

  protected void updateStatistics(double lowPrd, double highPrd) {
    GainExperiment gain = (GainExperiment) expResult;
    gain.setRangeForStatistics(lowPrd, highPrd);
  }

  protected void updateReference(int referenceIndex) {
    GainExperiment gain = (GainExperiment) expResult;
    gain.setReferenceIndex(referenceIndex);
  }

}
