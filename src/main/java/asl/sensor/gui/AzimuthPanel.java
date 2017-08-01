package asl.sensor.gui;

import java.awt.Color;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.text.DecimalFormat;
import java.util.List;

import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JSpinner;
import javax.swing.SpinnerModel;
import javax.swing.SpinnerNumberModel;
import javax.swing.SwingConstants;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.annotations.XYTitleAnnotation;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.PolarPlot;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.DefaultXYItemRenderer;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.chart.renderer.xy.XYSplineRenderer;
import org.jfree.chart.title.TextTitle;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.RectangleAnchor;

import asl.sensor.experiment.AzimuthExperiment;
import asl.sensor.experiment.ExperimentEnum;
import asl.sensor.input.DataStore;

/**
 * Wrapper class to display result from Azimuth. Overrides some parent
 * functions because the main plot uses polar orientation rather than typical
 * x-y plotting. 
 * @author akearns
 *
 */
public class AzimuthPanel extends ExperimentPanel {

  /**
   * 
   */
  private static final long serialVersionUID = 4088024342809622854L;
  private static final DecimalFormat df = new DecimalFormat("#.###");
  JSpinner offsetSpinner;
  JFreeChart angleChart, coherenceChart, estimChart;
  
  JComboBox<String> chartSelector;

  public AzimuthPanel(ExperimentEnum exp) {
    super(exp);
    
    SpinnerModel spinModel = new SpinnerNumberModel(0, -360, 360, 0.1);
    offsetSpinner = new JSpinner(spinModel);
    
    JLabel jbl = new JLabel("Offset angle (deg.):");
    jbl.setLabelFor(offsetSpinner);
    jbl.setHorizontalTextPosition(SwingConstants.RIGHT);
    jbl.setHorizontalAlignment(SwingConstants.RIGHT);
    JPanel labelPanel = new JPanel();
    labelPanel.add(jbl);
    
    chartSelector = new JComboBox<String>();
    chartSelector.addItem("Azimuth angle");
    chartSelector.addItem("Coherence");
    chartSelector.addItem("Estimation");
    chartSelector.setSelectedItem(0);
    chartSelector.addActionListener(this);
    
    plotTheseInBold = new String[]{}; // shouldn't be used anyway
    
    channelType[0] = "North test sensor";
    channelType[1] = "East test sensor";
    channelType[2] = "Reference sensor " + 
                     "(use offset to specify degrees from north)";
    
    // don't bother instantiating axes, we need to build a custom polar plot
    // and so will just use the chartfactory methods to do our building anyway    
    
    angleChart = ChartFactory.createPolarChart( expType.getName(), 
        null, false, false, false);
    chart = angleChart;
    chartPanel.setChart(chart);
    
    coherenceChart = 
        ChartFactory.createXYLineChart( expType.getName() + " Coherence",
        "Frequency (Hz)", "Coherence of best-fit angle", null);
    
    estimChart = 
        ChartFactory.createXYLineChart( expType.getName() + " Windowing",
        "Window start", "Coherence, Angle", null);
    
    this.setLayout( new GridBagLayout() );
    
    GridBagConstraints gbc = new GridBagConstraints();
    gbc.gridx = 0; gbc.gridy = 0;
    gbc.weightx = 1; gbc.weighty = 0;
    gbc.fill = GridBagConstraints.NONE;
    gbc.anchor = GridBagConstraints.EAST;
    
    gbc.anchor = GridBagConstraints.CENTER;
    gbc.gridx = 0; gbc.gridy = 0;
    gbc.gridwidth = 3;
    gbc.weightx = 1.0; gbc.weighty = 1.0;
    gbc.fill = GridBagConstraints.BOTH;
    this.add(chartPanel, gbc);
    
    gbc.weighty = 0.0;
    gbc.gridy += 1;
    gbc.gridwidth = 1;
    gbc.fill = GridBagConstraints.NONE;
    gbc.anchor = GridBagConstraints.EAST;
    this.add(jbl, gbc);
    
    gbc.gridx += 1;
    gbc.anchor = GridBagConstraints.WEST;
    this.add(offsetSpinner, gbc);
    
    gbc.gridx += 1;
    gbc.fill = GridBagConstraints.HORIZONTAL;
    gbc.anchor = GridBagConstraints.CENTER;
    this.add(chartSelector, gbc);
    
    gbc.gridx = 0; gbc.gridy += 1;
    gbc.gridwidth = 3;
    gbc.fill = GridBagConstraints.NONE;
    gbc.anchor = GridBagConstraints.CENTER;
    this.add(save, gbc);
  }
  
  @Override
  public void actionPerformed(ActionEvent e) {
    
    if (e.getSource() == chartSelector) {
      JFreeChart[] charts = getCharts();
      chart = charts[chartSelector.getSelectedIndex()];
      chartPanel.setChart(chart);
      return;
    }
    
    super.actionPerformed(e);
  }

  @Override
  protected void clearChartAndSetProgressData() {
    chartSelector.setSelectedIndex(0);
    angleChart = ChartFactory.createPolarChart( expType.getName(), 
        null, false, false, false);
    chart = angleChart;
    chartPanel.setChart(chart);
    displayInfoMessage("Running calculation...");
  }
  
  @Override
  public void displayErrorMessage(String errMsg) {
    
    if (chartSelector.getSelectedIndex() == 0) {
      PolarPlot plot = (PolarPlot) angleChart.getPlot();
      plot.clearCornerTextItems();
      plot.addCornerTextItem(errMsg);
    } else {
      super.displayInfoMessage(errMsg);
    }
    
  }

  public void displayInfoMessage(String infoMsg) {
    
    if (chartSelector.getSelectedIndex() == 0) {
      PolarPlot plot = (PolarPlot) angleChart.getPlot();
      plot.clearCornerTextItems();
      plot.addCornerTextItem(infoMsg);
    } else {
      super.displayInfoMessage(infoMsg);
    }
    
  }

  @Override
  protected void drawCharts() {
    chartSelector.setEnabled(true);
    chartSelector.setSelectedIndex(0);
    chart = angleChart;
    chartPanel.setChart(chart);
  }
  
  @Override
  public JFreeChart[] getCharts() {
    return new JFreeChart[]{angleChart, coherenceChart, estimChart};
  }
  
  @Override
  public String getInsetStrings() {
    AzimuthExperiment az = (AzimuthExperiment) expResult;
    double value = az.getOffset();
    double angle = az.getFitAngle();
    String angleStr = "FIT ANGLE: " + df.format(angle);
    double result = ( (value + angle) % 360 + 360) % 360;

    angleStr += " + " + df.format(value) + " = " + df.format(result);
    angleStr += " (+/- " + df.format( az.getUncertainty() ) + ")";
    return angleStr;
  }
  
  @Override
  public int panelsNeeded() {
    return 3;
  }
  
  @Override
  protected void updateData(DataStore ds) {
    
    set = true;
    
    double value = (double) offsetSpinner.getValue();
    
    if (value < 0) {
      value += 360;
    }
    
    AzimuthExperiment az = (AzimuthExperiment) expResult;
    az.setOffset(value);
    
    expResult.runExperimentOnData(ds);
    
    List<XYSeriesCollection> allData = expResult.getData();
    
    XYSeriesCollection polars = allData.get(0);
    XYSeriesCollection xysc = allData.get(1);
    
    coherenceChart = ChartFactory.createXYLineChart( 
        expType.getName() + " Coherence", "Frequency (Hz)", "Coherence", xysc);
    
    angleChart = ChartFactory.createPolarChart( expType.getName(),
        polars, true, true, false);
    
    String angleStr = getInsetStrings();
    
    XYPlot xyp = (XYPlot) coherenceChart.getPlot();
    TextTitle title = new TextTitle(angleStr);
    title.setBackgroundPaint(Color.white);
    XYTitleAnnotation xyt = new XYTitleAnnotation(0.98, 0.02, title,
        RectangleAnchor.BOTTOM_RIGHT);
    xyp.clearAnnotations();
    xyp.addAnnotation(xyt);
    // plot.addCornerTextItem(angleStr);
    
    PolarPlot plot = (PolarPlot) angleChart.getPlot();
    plot.clearCornerTextItems();
    plot.addCornerTextItem(angleStr);
    
    String titleEst = expType.getName() + " Accuracy Estimation";
    estimChart = ChartFactory.createXYLineChart( titleEst,
        "xaxis", "yaxis", xysc);
    XYSeriesCollection angleEstim = allData.get(2);
    XYSeriesCollection coherEstim = allData.get(3);
    xyp = estimChart.getXYPlot();
    xyp.setDataset(0, angleEstim);
    xyp.setDataset(1, coherEstim);
    xyp.setRenderer( 0, new DefaultXYItemRenderer() );
    // set color of second dataset to be blue
    XYItemRenderer renderer = new DefaultXYItemRenderer();
    renderer.setSeriesPaint(0, Color.BLUE);
    xyp.setRenderer(1, renderer);
    xyp.setRangeAxis( 0, new NumberAxis("Angle est. (deg)") );
    xyp.setRangeAxis( 1, new NumberAxis("Coherence est.") );
    NumberAxis xAx = new NumberAxis("Time window start");
    xAx.setAutoRangeIncludesZero(false);
    xyp.setDomainAxis(xAx);
    xyp.mapDatasetToRangeAxis(0, 0);
    xyp.mapDatasetToRangeAxis(1, 1);
    
    /*
    estimChart = 
        new JFreeChart(expType.getName() + " Windowing", getFont(), xyp, true);
    */
    
    chartSelector.setSelectedIndex(0);
  }
  
}
