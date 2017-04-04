package asl.sensor.gui;

import java.awt.Color;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
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
import org.jfree.chart.plot.PolarPlot;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.title.TextTitle;
import org.jfree.chart.title.Title;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.RectangleAnchor;

import asl.sensor.experiment.AzimuthExperiment;
import asl.sensor.experiment.ExperimentEnum;
import asl.sensor.input.DataStore;

public class AzimuthPanel extends ExperimentPanel {

  JSpinner offsetSpinner;
  JFreeChart angleChart, coherenceChart;
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
    chartSelector.setSelectedItem(0);
    chartSelector.addActionListener(this);
    
    plotTheseInBold = new String[]{}; // shouldn't be used anyway
    
    channelType[0] = "North test sensor";
    channelType[1] = "East test sensor";
    channelType[2] = "Reference sensor " + 
                     "(use offset to specify degrees from north)";
    
    
    angleChart = ChartFactory.createPolarChart( expType.getName(), 
        null, false, false, false);
    chart = angleChart;
    chartPanel.setChart(chart);
    
    coherenceChart = 
        ChartFactory.createXYLineChart( expType.getName() + " Coherence",
        "Frequency (Hz)", "Coherence", null);
    
    this.setLayout( new GridBagLayout() );
    
    GridBagConstraints gbc = new GridBagConstraints();
    gbc.gridx = 0; gbc.gridy = 0;
    gbc.weightx = 1; gbc.weighty = 0;
    gbc.fill = GridBagConstraints.NONE;
    gbc.anchor = GridBagConstraints.EAST;
    /*
    this.add(jbl, gbc);
    
    gbc.gridx += 1;
    gbc.fill = GridBagConstraints.NONE;
    gbc.anchor = GridBagConstraints.WEST;
    this.add(offsetSpinner, gbc);
    */

    
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

  /**
   * 
   */
  private static final long serialVersionUID = 4088024342809622854L;

  @Override
  public void actionPerformed(ActionEvent e) {
    
    if (e.getSource() == chartSelector) {
      if ( chartSelector.getSelectedIndex() == 0 ) {
        chart = angleChart;
      } else {
        chart = coherenceChart;
      }
      chartPanel.setChart(chart);
      return;
    }
    
    super.actionPerformed(e);
  }
  
  @Override
  public void updateData(DataStore ds) {
    // TODO Auto-generated method stub
    
    double value = (double) offsetSpinner.getValue();
    
    if (value < 0) {
      value += 360;
    }
    
    AzimuthExperiment az = (AzimuthExperiment) expResult;
    az.setOffset(value);
    
    expResult.setData(ds);
    
    List<XYSeriesCollection> allData = expResult.getData();
    
    XYSeriesCollection polars = allData.get(0);
    XYSeriesCollection xysc = allData.get(1);
    
    coherenceChart = ChartFactory.createXYLineChart( 
        expType.getName() + " Coherence", "Frequency (Hz)", "Coherence", xysc);
    
    angleChart = ChartFactory.createPolarChart( expType.getName(),
        polars, true, true, false);
    
    double angle = az.getFitAngle();
    String angleStr = "FIT ANGLE: " + angle;
    double result = ( (value + angle) % 360 + 360) % 360;
    
    angleStr += " + " + value + " = " + result;
    
    XYPlot xyp = (XYPlot) coherenceChart.getPlot();
    TextTitle title = new TextTitle(angleStr);
    title.setBackgroundPaint(Color.white);
    XYTitleAnnotation xyt = new XYTitleAnnotation(0.98, 0.02, title,
        RectangleAnchor.BOTTOM_RIGHT);
    xyp.clearAnnotations();
    xyp.addAnnotation(xyt);
    // plot.addCornerTextItem(angleStr);
    
    PolarPlot plot = (PolarPlot) angleChart.getPlot();
    plot.addCornerTextItem(angleStr);
    
    if ( chartSelector.getSelectedIndex() == 0 ) {
      chart = angleChart;
    } else {
      chart = coherenceChart;
    }
    
    chartPanel.setChart(chart);
  }

  @Override
  public int panelsNeeded() {
    return 3;
  }

}
