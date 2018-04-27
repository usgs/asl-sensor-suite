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
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.title.TextTitle;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.RectangleAnchor;
import asl.sensor.experiment.ExperimentEnum;
import asl.sensor.experiment.SineExperiment;
import asl.sensor.input.DataStore;

public class SinePanel extends ExperimentPanel {

  /**
   *
   */
  private static final long serialVersionUID = 2453757553804095685L;

  private JFreeChart sinesChart, linearChart;
  private NumberAxis calAxis, outAxis;
  private JComboBox<String> plotSelection;

  public SinePanel(ExperimentEnum exp) {
    super(exp);
    channelType[0] = "Calibration input";
    channelType[1] = "Calibration output (no RESP needed)";
    String xAxisTitle = "Time (s)";
    String yAxisTitle = "Normalized sine wave signal";
    xAxis = new NumberAxis(xAxisTitle);
    xAxis.setAutoRange(true);
    yAxis = new NumberAxis(yAxisTitle);
    yAxis.setAutoRange(true);
    Font bold = xAxis.getLabelFont().deriveFont(Font.BOLD);
    xAxis.setLabelFont(bold);
    yAxis.setLabelFont(bold);
    applyAxesToChart();

    calAxis = new NumberAxis("Calibration input sample");
    outAxis = new NumberAxis("Output sample");
    calAxis.setAutoRange(true);
    outAxis.setAutoRange(true);
    calAxis.setLabelFont(bold);
    outAxis.setLabelFont(bold);

    plotSelection = new JComboBox<String>();
    plotSelection.addItem("Sine plot");
    plotSelection.addItem("Linearity plot");
    plotSelection.addActionListener(this);
    plotSelection.setEnabled(false);

    this.setLayout(new GridBagLayout());
    GridBagConstraints gbc = new GridBagConstraints();

    gbc.fill = GridBagConstraints.BOTH;
    gbc.gridx = 0;
    gbc.gridy = 0;
    gbc.weightx = 1.0;
    gbc.weighty = 1.0;
    gbc.gridwidth = 3;
    gbc.anchor = GridBagConstraints.CENTER;
    this.add(chartPanel, gbc);

    // place the other UI elements in a single row below the chart
    gbc.gridwidth = 1;
    gbc.weighty = 0.0;
    gbc.weightx = 0.0;
    gbc.anchor = GridBagConstraints.WEST;
    gbc.fill = GridBagConstraints.NONE;
    gbc.gridy += 1;
    gbc.gridx = 0;
    JPanel space = new JPanel();
    space.setPreferredSize(plotSelection.getPreferredSize());
    this.add(space, gbc);

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
    this.add(plotSelection, gbc);
  }

  @Override
  public void actionPerformed(ActionEvent e) {

    super.actionPerformed(e);

    if (e.getSource() == plotSelection) {

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
          new JFreeChart[]{sinesChart, linearChart};
      chart = charts[idx];
      chartPanel.setChart(chart);

      return;

    }

  }

  @Override
  protected void drawCharts() {
    plotSelection.setSelectedIndex(0);
    plotSelection.setEnabled(true);
    XYSeriesCollection xysc = expResult.getData().get(0);
    sinesChart = buildChart(xysc);
    XYPlot xyp = (XYPlot) sinesChart.getPlot();

    TextTitle result = new TextTitle();
    result.setText(getInsetStrings());
    result.setBackgroundPaint(Color.white);
    XYTitleAnnotation xyt = new XYTitleAnnotation(0.98, 0.98, result,
        RectangleAnchor.TOP_RIGHT);
    xyp.clearAnnotations();
    xyp.addAnnotation(xyt);

    xysc = expResult.getData().get(1);
    // may need to make this something other than an xylinechart
    linearChart = buildChart(xysc, calAxis, outAxis);
    appendChartTitle(linearChart, " (Linearity plot)");
    chart = sinesChart;
    chartPanel.setChart(chart);
    chartPanel.setMouseZoomable(true);
  }

  @Override
  public int panelsNeeded() {
    return expResult.blocksNeeded();
  }

  @Override
  protected void updateData(DataStore ds) {
    expResult.runExperimentOnData(ds);

    XYSeriesCollection xysc = expResult.getData().get(0);
    for (int i = 0; i < xysc.getSeriesCount(); ++i) {
      Color toColor = COLORS[i % COLORS.length];
      String curve = (String) xysc.getSeriesKey(i);
      seriesColorMap.put(curve, toColor);
    }
    set = true;
  }

  @Override
  public JFreeChart[] getCharts() {
    return new JFreeChart[] {sinesChart, linearChart};
  }

  /**
   * Used to get the text that will represent the title text in the PDF result
   */
  @Override
  public String getInsetStrings() {
    return getInsetString((SineExperiment) expResult);
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

  public static String getInsetString(SineExperiment rexp) {
    double cAmp = rexp.getCalAmplitude();
    double oAmp = rexp.getOutAmplitude();
    String calAmp = DECIMAL_FORMAT.get().format(cAmp);
    String outAmp = DECIMAL_FORMAT.get().format(oAmp);
    String ratio = DECIMAL_FORMAT.get().format(cAmp / oAmp);
    String pFreq = DECIMAL_FORMAT.get().format(rexp.getEstSineFreq());
    StringBuilder sb = new StringBuilder();
    sb.append("Calculated calibration amplitude: ");
    sb.append(calAmp);
    sb.append('\n');
    sb.append("Calculated output amplitude: ");
    sb.append(outAmp);
    sb.append('\n');
    sb.append("Amplitude ratio: ");
    sb.append(ratio);
    sb.append('\n');
    sb.append("Estimated sine frequency: ");
    sb.append(pFreq);
    return sb.toString();
  }

}
