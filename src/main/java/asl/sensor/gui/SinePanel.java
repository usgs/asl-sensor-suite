package asl.sensor.gui;

import asl.sensor.experiment.ExperimentEnum;
import asl.sensor.experiment.SineExperiment;
import asl.sensor.input.DataStore;
import java.awt.Color;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.text.DecimalFormat;
import org.jfree.chart.annotations.XYTitleAnnotation;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.title.TextTitle;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.RectangleAnchor;

public class SinePanel extends ExperimentPanel {

  /**
   *
   */
  private static final long serialVersionUID = 2453757553804095685L;

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

    this.setLayout(new GridBagLayout());
    GridBagConstraints gbc = new GridBagConstraints();
    gbc.gridx = 0;
    gbc.gridy = 0;
    gbc.weightx = 1.0;
    gbc.weighty = 1.0;
    gbc.fill = GridBagConstraints.BOTH;
    gbc.anchor = GridBagConstraints.CENTER;
    this.add(chartPanel, gbc);
    gbc.weighty = 0.0;
    gbc.fill = GridBagConstraints.NONE;
    gbc.gridy += 1;
    this.add(save, gbc);
  }

  @Override
  protected void drawCharts() {
    XYSeriesCollection xysc = expResult.getData().get(0);

    setChart(xysc);
    XYPlot xyp = (XYPlot) chart.getPlot();

    TextTitle result = new TextTitle();
    result.setText(getInsetStrings());
    result.setBackgroundPaint(Color.white);
    XYTitleAnnotation xyt = new XYTitleAnnotation(0.98, 0.98, result,
        RectangleAnchor.TOP_RIGHT);
    xyp.clearAnnotations();
    xyp.addAnnotation(xyt);

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

  /**
   * Used to get the text that will represent the title text in the PDF result
   */
  @Override
  public String getInsetStrings() {
    return getInsetString((SineExperiment) expResult);
  }

  public static String getInsetString(SineExperiment rexp) {
    DecimalFormat df = new DecimalFormat("#.#######");
    double cAmp = rexp.getCalAmplitude();
    double oAmp = rexp.getOutAmplitude();
    String calAmp = df.format(cAmp);
    String outAmp = df.format(oAmp);
    String ratio = df.format(cAmp / oAmp);
    String pFreq = df.format(rexp.getEstSineFreq());
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
