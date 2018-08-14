package asl.sensor.gui;

import asl.sensor.ExperimentFactory;
import asl.sensor.experiment.GainSixExperiment;
import asl.sensor.experiment.GainExperiment;
import asl.sensor.input.DataStore;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import javax.swing.JComboBox;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.annotations.XYTitleAnnotation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.title.TextTitle;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.RectangleAnchor;

/**
 * Six-input extension for RelativeGain gui
 *
 * @author akearns - KBRWyle
 */
public class GainSixPanel extends GainPanel {

  private static final long serialVersionUID = 6140615094847886109L;
  private final JComboBox<String> plotSelection; // which chart to display in window?
  private JFreeChart northChart;
  private JFreeChart eastChart;
  private JFreeChart verticalChart; // gain result per dimensional component
  /**
   * Instantiate the panel, including sliders and stat calc button
   */
  public GainSixPanel(ExperimentFactory experiment) {
    // instantiate common components
    super(experiment);
    // make sure the experiment is gain-six
    expResult = experiment.createExperiment();

    for (int i = 0; i < 2; ++i) {
      int num = i + 1;
      channelType[3 * i] = "North sensor " + num + " (RESP required)";
      channelType[(3 * i) + 1] = "East sensor " + num + " (RESP required)";
      channelType[(3 * i) + 2] = "Vertical sensor " + num + " (RESP required)";
    }

    plotTheseInBold = new String[]{"NLNM"};


    northChart =
        ChartFactory.createXYLineChart(expType.getName() + " (North)",
            "", "", null);
    eastChart =
        ChartFactory.createXYLineChart(expType.getName() + " (East)",
            "", "", null);
    verticalChart =
        ChartFactory.createXYLineChart(expType.getName() + " (Vertical)",
            "", "", null);
    for (JFreeChart chart : getCharts()) {
      chart.getXYPlot().setDomainAxis(getXAxis());
      chart.getXYPlot().setRangeAxis(getYAxis());
    }

    // create layout
    removeAll();
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

    constraints.fill = GridBagConstraints.HORIZONTAL;
    constraints.gridx += 1;
    constraints.weightx = 0;
    constraints.anchor = GridBagConstraints.WEST;
    plotSelection = new JComboBox<>();
    plotSelection.addItem("North component plot");
    plotSelection.addItem("East component plot");
    plotSelection.addItem("Vertical component plot");
    plotSelection.addActionListener(this);
    add(plotSelection, constraints);

  }

  /**
   * Calls functions to do replotting and stat recalculations when different
   * timeseries are selected or the recalculate button is hit
   */
  @Override
  public void actionPerformed(ActionEvent event) {
    if (event.getSource() == recalcButton) {
      setTitle();
      recalcButton.setEnabled(false);
      return;
    }

    super.actionPerformed(event); // saving?

    if (event.getSource() == plotSelection) {
      int index = plotSelection.getSelectedIndex();

      JFreeChart[] charts = getCharts();
      chart = charts[index];
      chartPanel.setChart(chart);
      chartPanel.setMouseZoomable(true);
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

    GainSixExperiment gainBackend = (GainSixExperiment) expResult;
    final int referenceIndex = referenceSeries.getSelectedIndex();
    final int index1 = (referenceIndex + 1) % 2;
    int leftSliderValue, rightSliderValue;

    double[] minMax = gainBackend.getMinMaxFrequencies();
    updateReference(referenceIndex);

    // since intervals of incoming data match, so too limits of plot
    // this is used in mapping scale of slider to x-axis values
    lowPeriod = Math.log10(minMax[0]); // value when slider is 0
    highPeriod = Math.log10(minMax[1]); // value when slider is 1000
    leftSliderValue = mapPeriodToSlider(GainExperiment.DEFAULT_LOW_BOUND);
    rightSliderValue = mapPeriodToSlider(GainExperiment.DEFAULT_UP_BOUND);

    // plot has 3 components: source, destination, NLNM line plot
    JFreeChart[] charts = getCharts();
    for (int i = 0; i < charts.length; ++i) {
      XYSeriesCollection timeseriesIn = expResult.getData().get(i);
      XYSeriesCollection timeseries = new XYSeriesCollection();
      timeseries.addSeries(timeseriesIn.getSeries(referenceIndex));
      timeseries.addSeries(timeseriesIn.getSeries(index1));
      timeseries.addSeries(timeseriesIn.getSeries("NLNM"));

      charts[i] = buildChart(timeseries);

      // set vertical bars and enable sliders
      setSliderValues(leftSliderValue, rightSliderValue);

      // set the domain to match the boundaries of the octave centered at peak
      charts[i] = setDomainMarkers(GainExperiment.DEFAULT_LOW_BOUND,
          GainExperiment.DEFAULT_UP_BOUND, charts[i]);

      // and now set the sliders to match where that window is
      leftSlider.setEnabled(true);
      rightSlider.setEnabled(true);

    }

    // make sure pointers to each chart are set properly
    northChart = charts[0];
    eastChart = charts[1];
    verticalChart = charts[2];

    // lastly, display the calculated statistics in a textbox in the corner
    setTitle();

    chart = charts[plotSelection.getSelectedIndex()];

    // obviously, set the chart
    chartPanel.setChart(chart);
    chartPanel.setMouseZoomable(false);

  }

  @Override
  public JFreeChart[] getCharts() {
    return new JFreeChart[]{northChart, eastChart, verticalChart};
  }

  @Override
  public int panelsNeeded() {
    return 6;
  }

  /**
   * Used to populate the comboboxes with the incoming data
   */
  private void setDataNames() {

    referenceSeries.setEnabled(false);

    referenceSeries.removeAllItems();

    referenceSeries.addItem("Data from sensor set 1");
    referenceSeries.addItem("Data from sensor set 2");

    referenceSeries.setSelectedIndex(0);
  }

  /**
   * Displays the statistic results when the calculate button is hit
   * in an inset box on the chart, also used as text in report generation
   */
  private void setTitle() {
    JFreeChart[] charts = getCharts();
    String[] results = expResult.getInsetStrings();
    for (int i = 0; i < charts.length; ++i) {
      XYPlot plot = charts[i].getXYPlot();
      TextTitle result = getDefaultTextTitle();
      result.setText(results[i]);
      XYTitleAnnotation title = new XYTitleAnnotation(0.98, 0.98, result,
          RectangleAnchor.TOP_RIGHT);
      plot.clearAnnotations();
      plot.addAnnotation(title);
    }

  }

  @Override
  protected void updateData(final DataStore dataStore) {

    set = true;

    setDataNames();

    expResult.runExperimentOnData(dataStore);

    // need to have 2 series for relative gain
    referenceSeries.setEnabled(true);

  }

  @Override
  protected void updateStatistics(double lowPrd, double highPrd) {
    GainSixExperiment gain = (GainSixExperiment) expResult;
    gain.setRangeForStatistics(lowPrd, highPrd);
  }

  @Override
  void updateReference(int referenceIndex) {
    GainSixExperiment gain = (GainSixExperiment) expResult;
    gain.setReferenceIndex(referenceIndex);
  }

}
