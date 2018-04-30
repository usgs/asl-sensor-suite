package asl.sensor.gui;

import asl.sensor.experiment.ExperimentEnum;
import asl.sensor.experiment.StepExperiment;
import asl.sensor.input.DataStore;
import java.awt.Color;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import javax.swing.JComboBox;
import javax.swing.JPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.DateAxis;
import org.jfree.chart.axis.LogarithmicAxis;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.block.BlockContainer;
import org.jfree.chart.block.FlowArrangement;
import org.jfree.chart.title.CompositeTitle;
import org.jfree.chart.title.TextTitle;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.RectangleEdge;
import org.jfree.ui.VerticalAlignment;

/**
 * Holds the plot results of a step experiment. Gets the timeseries data from it
 * as well as the corner and damping values gotten in the process of its
 * calculations. Aside from common elements to ExperimentPanels, also sets an
 * inset displaying parameters derived from backend fit calculation
 *
 * @author akearns
 */
public class StepPanel extends ExperimentPanel {

  private static final long serialVersionUID = 3693391540945130688L;
  private static final int TITLE_IDX = 0;

  /**
   * Static helper method for getting the formatted inset string directly
   * from a StepExperiment
   *
   * @param experiment StepExperiment with data to be extracted
   * @return String format representation of data from the experiment
   */
  private static String getInsetString(StepExperiment experiment) {
    String[] strings = getInsetStringList(experiment);
    StringBuilder sb = new StringBuilder();
    for (String str : strings) {
      sb.append(str);
      sb.append("\n");
    }
    return sb.toString();
  }

  private static String[] getInsetStringList(StepExperiment experiment) {
    double[] rolloff = experiment.getInitParams();
    double[] fit = experiment.getFitParams();
    double corner = rolloff[0];
    double damping = rolloff[1];
    double fitCorner = fit[0];
    double fitDamping = fit[1];

    double cornerPrd = 1. / corner;
    double fitCornerPrd = 1. / fitCorner;

    String sb = "RESP parameters"
        + "\nCorner frequency (Hz): "
        + DECIMAL_FORMAT.get().format(corner)
        + " ("
        + DECIMAL_FORMAT.get().format(cornerPrd)
        + " secs)"
        + "\nDamping: "
        + DECIMAL_FORMAT.get().format(damping)
        + "\n";

    String sb2 = "Best-fit parameters"
        + "\nCorner frequency (Hz): "
        + DECIMAL_FORMAT.get().format(fitCorner)
        + " ("
        + DECIMAL_FORMAT.get().format(fitCornerPrd)
        + " secs)"
        + "\nDamping: "
        + DECIMAL_FORMAT.get().format(fitDamping)
        + "\n";
    return new String[]{sb, sb2};
  }

  private final JComboBox<String> plotSelection;
  private JFreeChart stepChart, magChart, phaseChart;
  private final ValueAxis freqAxis;
  private final ValueAxis magAxis;
  private final ValueAxis phaseAxis;

  StepPanel(ExperimentEnum experiment) {
    super(experiment);

    channelType[0] = "Calibration input";
    channelType[1] = "Calibration output from sensor (RESP required)";

    String xAxisTitle = "Time";
    String yAxisTitle = "Normalized counts";
    xAxis = new DateAxis(xAxisTitle);
    ((DateAxis) xAxis).setDateFormatOverride(ExperimentPanel.DATE_TIME_FORMAT.get());
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

    plotSelection = new JComboBox<>();
    plotSelection.addItem("Step function");
    plotSelection.addItem("Response amplitude");
    plotSelection.addItem("Response phase");
    plotSelection.addActionListener(this);

    applyAxesToChart();

    this.setLayout(new GridBagLayout());
    GridBagConstraints constraints = new GridBagConstraints();
    constraints.gridx = 0;
    constraints.gridy = 0;
    constraints.weightx = 1.0;
    constraints.weighty = 1.0;
    constraints.gridwidth = 3;
    constraints.fill = GridBagConstraints.BOTH;
    constraints.anchor = GridBagConstraints.CENTER;
    this.add(chartPanel, constraints);

    // add empty space on left side to space out other components
    JPanel space = new JPanel();
    space.setMaximumSize(plotSelection.getMaximumSize());
    space.setPreferredSize(plotSelection.getPreferredSize());
    constraints.weighty = 0.0;
    constraints.weightx = 1.0;
    constraints.fill = GridBagConstraints.BOTH;
    constraints.gridwidth = 1;
    constraints.gridy += 1;
    this.add(space, constraints);

    constraints.fill = GridBagConstraints.NONE;
    constraints.anchor = GridBagConstraints.CENTER;
    constraints.gridx += 1;
    constraints.weightx = 0.0;
    this.add(save, constraints);

    constraints.weightx = 1.0;
    constraints.gridx += 1;
    constraints.anchor = GridBagConstraints.LINE_END;
    this.add(plotSelection, constraints);

    plotTheseInBold = new String[]{};


  }

  @Override
  public void actionPerformed(ActionEvent event) {

    super.actionPerformed(event);

    if (event.getSource() == plotSelection) {
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
    }
  }

  @Override
  protected void drawCharts() {
    JFreeChart[] charts = getCharts();
    int index = plotSelection.getSelectedIndex();
    chart = charts[index];

    chartPanel.setChart(charts[index]);
    chartPanel.restoreAutoBounds();
    chartPanel.validate();
  }

  @Override
  public JFreeChart[] getCharts() {
    return new JFreeChart[]{stepChart, magChart, phaseChart};
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
   * Used to get the text that will populate the inset box for the plots
   *
   * @return String to place in TextTitle
   */
  @Override
  public String getInsetStrings() {

    return getInsetString((StepExperiment) expResult);

  }

  @Override
  public String getMetadataString() {
    StepExperiment experiment = (StepExperiment) expResult;
    StringBuilder sb = new StringBuilder();
    sb.append("Residuals:\n");
    double[] residuals = experiment.getResiduals();
    sb.append("Initial:  ");
    sb.append(residuals[0]);
    sb.append("\nFit:  ");
    sb.append(residuals[1]);
    sb.append('\n');
    sb.append(super.getMetadataString());
    return sb.toString();
  }

  @Override
  public ValueAxis getXAxis() {
    if (null == plotSelection) {
      return xAxis;
    }

    ValueAxis[] array = new ValueAxis[]{xAxis, freqAxis, freqAxis};
    return array[plotSelection.getSelectedIndex()];
  }

  @Override
  public ValueAxis getYAxis() {
    if (null == plotSelection) {
      return yAxis;
    }

    ValueAxis[] array = new ValueAxis[]{yAxis, magAxis, phaseAxis};
    return array[plotSelection.getSelectedIndex()];
  }

  @Override
  public int panelsNeeded() {
    return 2;
  }

  private void setSubtitles() {
    BlockContainer container = new BlockContainer(new FlowArrangement());
    CompositeTitle title = new CompositeTitle(container);
    String[] insets = getInsetStringList((StepExperiment) expResult);
    for (String inset : insets) {
      TextTitle result = new TextTitle();
      result.setText(inset);
      result.setBackgroundPaint(Color.white);
      container.add(result);
    }

    title.setVerticalAlignment(VerticalAlignment.BOTTOM);
    title.setPosition(RectangleEdge.BOTTOM);
    for (JFreeChart chart : getCharts()) {
      chart.addSubtitle(TITLE_IDX, title);
    }
  }

  /**
   * Pass in and retrieve data from the step experiment backend, to plot;
   * this is both the timeseries data as well as a title inset displaying
   * the parameters used in the plot calculations
   */
  @Override
  protected void updateData(final DataStore dataStore) {

    set = true;

    expResult.runExperimentOnData(dataStore);

    XYSeriesCollection stepData = expResult.getData().get(0);
    stepChart = buildChart(stepData, xAxis, yAxis);

    XYSeriesCollection magData = expResult.getData().get(1);
    magChart = buildChart(magData, freqAxis, magAxis);

    XYSeriesCollection phaseData = expResult.getData().get(2);
    phaseChart = buildChart(phaseData, freqAxis, phaseAxis);

    setSubtitles();
  }

}
