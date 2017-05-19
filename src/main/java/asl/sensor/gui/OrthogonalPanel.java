package asl.sensor.gui;

import java.awt.Color;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.util.Arrays;

import org.jfree.chart.annotations.XYTitleAnnotation;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.title.TextTitle;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.RectangleAnchor;

import asl.sensor.experiment.ExperimentEnum;
import asl.sensor.experiment.OrthogonalExperiment;
import asl.sensor.input.DataStore;

/**
 * Panel to display results of Orthogonal Experiment.
 * This plots the difference between reference and rotated sensor signals
 * and also displays the azimuth angles in an inset box.
 * No additional interface components are created for this panel.
 * @author akearns
 *
 */
public class OrthogonalPanel extends ExperimentPanel {

  /**
   * 
   */
  private static final long serialVersionUID = -2749224338484110043L;

  public static String getInsetString(OrthogonalExperiment ort) {  
    double[] fit = ort.getSolutionParams();
    double angle = ort.getFitAngle();
    
    StringBuilder sb = new StringBuilder();
    sb.append("Calculated angle between non-reference sensors:\n");
    sb.append(angle);
    sb.append('\n');
    sb.append("Offset angles for LH1 and LH2 sensor outputs:\n");
    sb.append( Arrays.toString(fit) );
    
    return sb.toString();
  }

  public OrthogonalPanel(ExperimentEnum exp) {
    super(exp);
    
    channelType[0] = "North reference sensor";
    channelType[1] = "East reference sensor";
    channelType[2] = "Assumed-north test sensor";
    channelType[3] = "Assumed-east test sensor";
    
    String xAxisTitle = "Time (s)";
    String yAxisTitle = "Amplitude difference (counts)";
    xAxis = new NumberAxis(xAxisTitle);
    xAxis.setAutoRange(true);
    //SimpleDateFormat sdf = new SimpleDateFormat("HH:mm");
    //sdf.setTimeZone( TimeZone.getTimeZone("UTC") );
    //xAxis.setLabel("UTC Time");
    yAxis = new NumberAxis(yAxisTitle);
    // yAxis.setAutoRange(true);
    Font bold = xAxis.getLabelFont().deriveFont(Font.BOLD);
    xAxis.setLabelFont(bold);
    yAxis.setLabelFont(bold);
    
    applyAxesToChart();
    
    this.setLayout( new GridBagLayout() );
    GridBagConstraints gbc = new GridBagConstraints();
    gbc.gridx = 0; gbc.gridy = 0;
    gbc.weightx = 1.0; gbc.weighty = 1.0;
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
    result.setText( getInsetString() );
    result.setBackgroundPaint(Color.white);
    XYTitleAnnotation xyt = new XYTitleAnnotation(0.98, 0.98, result,
        RectangleAnchor.TOP_RIGHT);
    xyp.clearAnnotations();
    xyp.addAnnotation(xyt);
    
    chartPanel.setChart(chart);
    
  }
  
  @Override
  public String getInsetString() {
    return getInsetString( (OrthogonalExperiment) expResult );
  }
    
    
  @Override
  public int panelsNeeded() {
    return 4;
  }
  
  @Override
  protected void updateData(final DataStore ds) {
    
    set = true;
    
    expResult.setData(ds);
  }

}
