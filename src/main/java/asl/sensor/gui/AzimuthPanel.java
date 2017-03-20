package asl.sensor.gui;

import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;

import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JSpinner;
import javax.swing.SpinnerModel;
import javax.swing.SpinnerNumberModel;
import javax.swing.SwingConstants;

import org.apache.commons.math3.fitting.leastsquares.MultivariateJacobianFunction;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PolarPlot;

import asl.sensor.experiment.AzimuthExperiment;
import asl.sensor.experiment.ExperimentEnum;
import asl.sensor.input.DataStore;

public class AzimuthPanel extends ExperimentPanel {

  JSpinner offsetSpinner;
  
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
    
    plotTheseInBold = new String[]{}; // shouldn't be used anyway
    
    channelType[0] = "North test sensor";
    channelType[1] = "East test sensor";
    channelType[2] = "Reference sensor " + 
                     "(use offset to specify degrees from north)";
    
    
    JFreeChart chart = ChartFactory.createPolarChart( expType.getName(), 
        null, false, false, false);
    chartPanel.setChart(chart);
    
    this.setLayout( new GridBagLayout() );
    GridBagConstraints gbc = new GridBagConstraints();
    gbc.gridx = 0; gbc.gridy = 0;
    gbc.weightx = 1; gbc.weighty = 0;
    gbc.fill = GridBagConstraints.NONE;
    gbc.anchor = GridBagConstraints.EAST;
    this.add(jbl, gbc);
    
    gbc.gridx += 1;
    gbc.fill = GridBagConstraints.NONE;
    gbc.anchor = GridBagConstraints.WEST;
    this.add(offsetSpinner, gbc);
    
    gbc.anchor = GridBagConstraints.CENTER;
    gbc.gridx = 0; gbc.gridy = 0;
    gbc.gridwidth = 2;
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
    
    gbc.gridx = 0; gbc.gridy += 1;
    gbc.gridwidth = 2;
    gbc.anchor = GridBagConstraints.CENTER;
    this.add(save, gbc);
  }

  /**
   * 
   */
  private static final long serialVersionUID = 4088024342809622854L;

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
    
    chart = ChartFactory.createPolarChart( expType.getName(),
        expResult.getData(), true, true, false);
    
    double angle = az.getFitAngle();
    String angleStr = "FIT ANGLE: -" + angle;
    double result = (value - angle) % 360;
    
    angleStr += " + " + value + " = " + result;
    
    PolarPlot plot = (PolarPlot) chart.getPlot();
    plot.addCornerTextItem(angleStr);
    
    chartPanel.setChart(chart);
    
  }

  @Override
  public int panelsNeeded() {
    return 3;
  }

}
