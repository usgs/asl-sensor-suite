package asl.sensor.gui;

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
    double[] rolloff = sp.getCornerAndDamping();
    double[] fit = sp.getFitCornerAndDamping();
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

  public StepPanel(ExperimentEnum exp) {
    super(exp);
    
    channelType[0] = "Calibration input";
    channelType[1] = "Calibration output from sensor (RESP required)";
    
    xAxisTitle = "Time (s)";
    yAxisTitle = "Normalized counts";
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
    
    plotTheseInBold = new String[]{};
    
    
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
    sb.append("LOADED RESPONSE:");
    sb.append('\n');
    sb.append( stex.getResponseName() );
    return sb.toString();
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
    
    set = true;
    
    displayInfoMessage("Running stepcal testing...");
    
    expResult.setData(ds);
    XYSeriesCollection xysc = expResult.getData().get(0);
    
    // here's the stuff that needs to stay here, not moved to experiment class
    setChart(xysc);

    TextTitle result = new TextTitle();
    result.setText( getInsetString() );
    result.setBackgroundPaint(Color.white);
    XYTitleAnnotation xyt = new XYTitleAnnotation(0.98, 0.5, result,
        RectangleAnchor.RIGHT);
    XYPlot xyp = (XYPlot) chart.getPlot();
    xyp.clearAnnotations();
    xyp.addAnnotation(xyt);
    
    chartPanel.setChart(chart);
  }

}
