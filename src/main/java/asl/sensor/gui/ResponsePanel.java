package asl.sensor.gui;

import java.awt.Color;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;

import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JPanel;

import org.jfree.chart.axis.LogarithmicAxis;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.data.xy.XYSeriesCollection;

import asl.sensor.experiment.ExperimentEnum;
import asl.sensor.experiment.ResponseExperiment;
import asl.sensor.input.DataStore;

public class ResponsePanel extends ExperimentPanel {

  public ValueAxis freqAxis, degreeAxis;
  public String freqAxisTitle, degreeAxisTitle;
  public JCheckBox freqSpaceBox;
  public JComboBox<String> plotSelection;
  
  public static final String MAGNITUDE = ResponseExperiment.MAGNITUDE;
  public static final String ARGUMENT = ResponseExperiment.ARGUMENT;
  
  public ResponsePanel(ExperimentEnum exp) {
    super(exp);
    // TODO Auto-generated constructor stub
    channelType[0] = "Response data (SEED data not used)";
    
    xAxisTitle = "Period (s)";
    freqAxisTitle = "Frequency (Hz)";
    yAxisTitle = "10 * log10( RESP(f) )";
    degreeAxisTitle = "phi(RESP(f))";
    xAxis = new LogarithmicAxis(xAxisTitle);
    freqAxis = new LogarithmicAxis(freqAxisTitle);
    yAxis = new NumberAxis(yAxisTitle);
    yAxis.setAutoRange(true);
    
    degreeAxis = new NumberAxis(degreeAxisTitle);
    
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
    gbc.fill = GridBagConstraints.NONE;
    gbc.gridx += 1;
    gbc.weightx = 0;
    gbc.anchor = GridBagConstraints.WEST;
    plotSelection = new JComboBox<String>();
    plotSelection.addItem(MAGNITUDE);
    plotSelection.addItem(ARGUMENT);
    this.add(plotSelection, gbc);
    
    seriesColorMap.put(MAGNITUDE, Color.RED);
    seriesColorMap.put(ARGUMENT, Color.BLUE);
    
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
    // TODO Auto-generated method stub
    boolean freqSpace = freqSpaceBox.isEnabled();
    expResult.setData(ds, freqSpace);
    
    XYSeriesCollection xysc = new XYSeriesCollection();
    XYSeriesCollection fromExp = (XYSeriesCollection) expResult.getData();
    
    if ( plotSelection.getSelectedItem().equals(MAGNITUDE) ) {
      xysc.addSeries( fromExp.getSeries(MAGNITUDE) );
    } else {
      xysc.addSeries( fromExp.getSeries(ARGUMENT) );
    }
    
    populateChart( xysc );

    chartPanel.setChart(chart);
    chartPanel.setMouseZoomable(true);
    
  }

  @Override
  public int panelsNeeded() {
    return 1;
  }

}
