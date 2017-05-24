package asl.sensor.gui;

import java.awt.Color;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;

import javax.swing.JComboBox;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.JFreeChart;
import org.jfree.data.xy.XYSeriesCollection;

import asl.sensor.experiment.ExperimentEnum;
import asl.sensor.experiment.ExperimentFactory;
import asl.sensor.experiment.NoiseNineExperiment;
import asl.sensor.input.DataStore;

/**
 * Panel for 9-input self noise. Similar to 3-input self noise (NoisePanel)
 * but includes multiple plots, one for each linear axis in 3D space
 * (north-south, east-west, up-down) and a combo box to select them
 * @author akearns
 *
 */
public class NoiseNinePanel extends NoisePanel {
  
  /**
   * 
   */
  private static final long serialVersionUID = -8049021432657749975L;
  JComboBox<String> plotSelection;
  boolean set;
  
  JFreeChart northChart, eastChart, vertChart;

  public NoiseNinePanel(ExperimentEnum exp) {
    super(exp);
    
    expResult = ExperimentFactory.createExperiment(exp);
    
    set = false;
    
    for (int i = 0; i < 3; ++i) {
      int num = i + 1;
      channelType[3 * i] = "North sensor " + num + " (RESP required)";
      channelType[(3 * i) + 1] = "East sensor " + num  + " (RESP required)";
      channelType[(3 * i) + 2] = "Vertical sensor " + num + " (RESP required)";
    }
    
    this.setLayout( new GridBagLayout() );
    GridBagConstraints gbc = new GridBagConstraints();
    
    String xTitle = getXAxis().getLabel();
    String yTitle = getYAxis().getLabel();
    
    northChart = 
        ChartFactory.createXYLineChart( expType.getName() + " (North)",
        xTitle, yTitle, null);
    eastChart = 
        ChartFactory.createXYLineChart( expType.getName() + " (East)",
        xTitle, yTitle, null);
    vertChart = 
        ChartFactory.createXYLineChart( expType.getName() + " (Vertical)",
        xTitle, yTitle, null);
    
    chart = northChart;
    chartPanel.setChart(chart);
    
    removeAll(); // get rid of blank spacer jpanel from super
    // (everything else will get redrawn)
    
    gbc.fill = GridBagConstraints.BOTH;
    gbc.gridx = 0; gbc.gridy = 0;
    gbc.weightx = 1.0; gbc.weighty = 1.0;
    gbc.gridwidth = 3;
    gbc.anchor = GridBagConstraints.CENTER;
    add(chartPanel, gbc);
    
    // place the other UI elements in a single row below the chart
    gbc.gridwidth = 1;
    gbc.weighty = 0.0; gbc.weightx = 0.0;
    gbc.anchor = GridBagConstraints.WEST;
    gbc.fill = GridBagConstraints.NONE;
    gbc.gridy += 1; gbc.gridx = 0;
    add(freqSpaceBox, gbc);
    
    gbc.gridx += 1;
    gbc.weightx = 1.0;
    gbc.fill = GridBagConstraints.NONE;
    gbc.anchor = GridBagConstraints.CENTER;
    // gbc.gridwidth = GridBagConstraints.REMAINDER;
    add(save, gbc);
    
    // combo box to select items
    gbc.fill = GridBagConstraints.HORIZONTAL;
    gbc.gridx += 1;
    gbc.weightx = 0;
    gbc.anchor = GridBagConstraints.WEST;
    plotSelection = new JComboBox<String>();
    plotSelection.addItem("North component plot");
    plotSelection.addItem("East component plot");
    plotSelection.addItem("Vertical component plot");
    plotSelection.addActionListener(this);
    add(plotSelection, gbc);
    
    revalidate();
    
  }

  @Override
  public void actionPerformed(ActionEvent e) {
    
    super.actionPerformed(e);
    
    if ( e.getSource() == plotSelection ) {
      
      int idx = plotSelection.getSelectedIndex();
      if (idx == 0) {
        chart = northChart;
      } else if (idx == 1){
        chart = eastChart;
      } else {
        chart = vertChart;
      }
      
      chartPanel.setChart(chart);
      chartPanel.setMouseZoomable(true);
      
      return;
      
    }
    
  }
  
  protected void drawCharts() {
    
    int idx = plotSelection.getSelectedIndex();
    
    if (idx == 0) {
      chart = northChart;
    } else if (idx == 1){
      chart = eastChart;
    } else {
      chart = vertChart;
    }

    chartPanel.setChart(chart);
    chartPanel.setMouseZoomable(true);

  }
  
  @Override
  public JFreeChart[] getCharts() {
    return new JFreeChart[]{northChart, eastChart, vertChart};
  }
  
  @Override
  public int panelsNeeded() {
    return 9;
  }
  
  @Override
  public void updateData(final DataStore ds) {
    
    set = true;
    
    boolean freqSpace = freqSpaceBox.isSelected();
    
    final boolean freqSpaceImmutable = freqSpace;

    NoiseNineExperiment noisExp = (NoiseNineExperiment) expResult;
    noisExp.setFreqSpace(freqSpaceImmutable);
    
    expResult.runExperimentOnData(ds);
    
    for (int j = 0; j < 3; ++j) {
      XYSeriesCollection xysc = expResult.getData().get(j);

      for (int i = 0; i < NOISE_PLOT_COUNT; ++i) {
        String name = (String) xysc.getSeriesKey(i);
        System.out.println(name);
        Color plotColor = COLORS[i % 3];
        seriesColorMap.put(name, plotColor);
        System.out.println(name+","+plotColor);
        if (i >= 3) {
          System.out.println(name+","+i);
          seriesDashedSet.add(name);
        }

      }
    }

    set = true;

    System.out.println("Charts being set!");
    
    northChart = buildChart( expResult.getData().get(0) );
    northChart.setTitle("Self-noise (NORTH)");
    System.out.println("North chart set!");
    
    eastChart = buildChart( expResult.getData().get(1) );
    eastChart.setTitle("Self-noise (EAST)");
    System.out.println("East chart set!");
    
    vertChart = buildChart( expResult.getData().get(2) );
    vertChart.setTitle("Self-noise (VERTICAL)");
    System.out.println("Vert chart set!");
  
  }
  
}
