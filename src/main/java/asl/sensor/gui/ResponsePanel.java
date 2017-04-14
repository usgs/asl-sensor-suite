package asl.sensor.gui;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.image.BufferedImage;
import java.util.HashMap;
import java.util.List;

import javax.swing.JCheckBox;
import javax.swing.JComboBox;

import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.LogarithmicAxis;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.data.xy.XYSeriesCollection;

import asl.sensor.experiment.ExperimentEnum;
import asl.sensor.experiment.ResponseExperiment;
import asl.sensor.input.DataStore;

/**
 * Panel used to display response curves. Includes plots of response magnitude
 * and argument (rotation in complex space). Because this includes two charts
 * and also does not require any input timeseries, there are more overridden
 * methods in this class than in most other panels, including the image output.
 * The chart to be displayed at any moment is selected with a combobox.
 * 
 * @author akearns
 *
 */
public class ResponsePanel extends ExperimentPanel {

  private ValueAxis freqAxis, degreeAxis;
  private String freqAxisTitle, degreeAxisTitle;
  private JCheckBox freqSpaceBox;
  private JComboBox<String> plotSelection;
  
  private final static Color[] COLOR_LIST = 
      new Color[]{Color.RED, Color.BLUE, Color.GREEN};
  
  private JFreeChart magChart, argChart;
  
  // values
  public static final String MAGNITUDE = ResponseExperiment.MAGNITUDE;
  public static final String ARGUMENT = ResponseExperiment.ARGUMENT;
  
  public ResponsePanel(ExperimentEnum exp) {
    super(exp);
    
    for (int i = 0; i < 3; ++i) {
      channelType[i] = "Response data (SEED data not used)";
    }
    
    xAxisTitle = "Period (s)";
    freqAxisTitle = "Frequency (Hz)";
    yAxisTitle = "10 * log10( RESP(f) )";
    degreeAxisTitle = "phi(RESP(f))";
    
    xAxis = new LogarithmicAxis(xAxisTitle);
    freqAxis = new LogarithmicAxis(freqAxisTitle);
    
    yAxis = new NumberAxis(yAxisTitle);
    yAxis.setAutoRange(true);
    
    degreeAxis = new NumberAxis(degreeAxisTitle);
    degreeAxis.setAutoRange(true);
    
    ( (NumberAxis) yAxis).setAutoRangeIncludesZero(false);
    Font bold = xAxis.getLabelFont().deriveFont(Font.BOLD);
    xAxis.setLabelFont(bold);
    yAxis.setLabelFont(bold);
    freqAxis.setLabelFont(bold);
    
    freqSpaceBox = new JCheckBox("Use Hz units (requires regen)");
    freqSpaceBox.setSelected(true);
    
    applyAxesToChart(); // now that we've got axes defined
    
    // set the GUI components
    this.setLayout( new GridBagLayout() );
    GridBagConstraints gbc = new GridBagConstraints();
    
    gbc.fill = GridBagConstraints.BOTH;
    gbc.gridx = 0; gbc.gridy = 0;
    gbc.weightx = 1.0; gbc.weighty = 1.0;
    gbc.gridwidth = 3;
    gbc.anchor = GridBagConstraints.CENTER;
    this.add(chartPanel, gbc);
    
    // place the other UI elements in a single row below the chart
    gbc.gridwidth = 1;
    gbc.weighty = 0.0; gbc.weightx = 0.0;
    gbc.anchor = GridBagConstraints.WEST;
    gbc.fill = GridBagConstraints.NONE;
    gbc.gridy += 1; gbc.gridx = 0;
    this.add(freqSpaceBox, gbc);
    
    gbc.gridx += 1;
    gbc.weightx = 1.0;
    gbc.fill = GridBagConstraints.NONE;
    gbc.anchor = GridBagConstraints.CENTER;
    // gbc.gridwidth = GridBagConstraints.REMAINDER;
    this.add(save, gbc);
    
    // add an empty panel as a spacer to keep the save button in the center
    gbc.fill = GridBagConstraints.HORIZONTAL;
    gbc.gridx += 1;
    gbc.weightx = 0;
    gbc.anchor = GridBagConstraints.WEST;
    plotSelection = new JComboBox<String>();
    plotSelection.addItem(MAGNITUDE);
    plotSelection.addItem(ARGUMENT);
    this.add(plotSelection, gbc);
    plotSelection.addActionListener(this);
    
  }

  @Override
  public int plotsToShow() {
    return 0;
  }
  
  @Override
  public ValueAxis getXAxis() {
    
    // true if using Hz units
    if ( freqSpaceBox.isSelected() ) {
        return freqAxis;
    }
    
    return xAxis;
    
  }
  
  @Override
  public ValueAxis getYAxis() {
    
    if ( null == plotSelection ) {
      return yAxis;
    }
    
    if ( plotSelection.getSelectedItem().equals(MAGNITUDE) ) {
      return yAxis;
    } else {
      return degreeAxis;
    }
  }
  
  /**
   * 
   */
  private static final long serialVersionUID = 1L;

  @Override
  public void updateData(DataStore ds) {

    seriesColorMap = new HashMap<String, Color>();
    
    boolean freqSpace = freqSpaceBox.isSelected();
    ResponseExperiment respExp = (ResponseExperiment) expResult;
    respExp.setFreqSpace(freqSpace);
    expResult.setData(ds);
    
    set = true;
    
    List<XYSeriesCollection> xysc = expResult.getData();
    XYSeriesCollection magSeries = xysc.get(0);
    XYSeriesCollection argSeries = xysc.get(1);
    
    for (int i = 0; i < magSeries.getSeriesCount(); ++i) {
        Color toColor = COLOR_LIST[i % COLOR_LIST.length];
        String magName = (String) magSeries.getSeriesKey(i);
        String argName = (String) argSeries.getSeriesKey(i);
        seriesColorMap.put(magName, toColor);
        seriesColorMap.put(argName, toColor);
    }
    
    int idx = plotSelection.getSelectedIndex();
    
    argChart = buildChart(argSeries);
    argChart.getXYPlot().setRangeAxis(degreeAxis);
    
    magChart = buildChart(magSeries);
    magChart.getXYPlot().setRangeAxis(yAxis);

    if (idx == 0) {
      chart = magChart;
    } else {
      chart = argChart;
    }
    chartPanel.setChart(chart);
    chartPanel.setMouseZoomable(true);
    
  }

  @Override
  public JFreeChart[] getCharts() {
    return new JFreeChart[]{magChart, argChart};
  }
  
  @Override
  public void actionPerformed(ActionEvent e) {
    
    super.actionPerformed(e);
    
    if ( e.getSource() == plotSelection ) {
      if (!set) {
        return;
      }
      
      int idx = plotSelection.getSelectedIndex();
      if (idx == 0) {
        chartPanel.setChart(magChart);
      } else {
        chartPanel.setChart(argChart);
      }
      
      return;
      
    }
    
  }
  
  @Override
  public int panelsNeeded() {
    return 3;
  }

}
