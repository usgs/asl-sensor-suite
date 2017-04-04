package asl.sensor.gui;

import java.awt.Color;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;

import javax.swing.JComboBox;
import javax.swing.JPanel;
import javax.swing.SwingWorker;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.JFreeChart;
import org.jfree.data.xy.XYSeriesCollection;

import asl.sensor.experiment.ExperimentEnum;
import asl.sensor.experiment.NoiseExperiment;
import asl.sensor.input.DataStore;

public class NoiseNinePanel extends NoisePanel {
  
  JComboBox<String> plotSelection;
  boolean set;
  JFreeChart northChart, eastChart, vertChart;
  
  public NoiseNinePanel(ExperimentEnum exp) {
    super(exp);
    
    set = false;
    
    for (int i = 0; i < 3; ++i) {
      int num = i + 1;
      channelType[3 * i] = "North sensor " + num + " (RESP required)";
      channelType[(3 * i) + 1] = "East sensor " + num  + " (RESP required)";
      channelType[(3 * i) + 2] = "Vertical sensor " + num + " (RESP required)";
    }
    
    this.setLayout( new GridBagLayout() );
    GridBagConstraints gbc = new GridBagConstraints();
    
    northChart = chart;
    eastChart = chart;
    vertChart = chart;
    
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
    
    // add an empty panel as a spacer to keep the save button in the center
    gbc.fill = GridBagConstraints.NONE;
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

  /**
   * 
   */
  private static final long serialVersionUID = -8049021432657749975L;

  @Override
  public void actionPerformed(ActionEvent e) {
    
    super.actionPerformed(e);
    
    if ( e.getSource() == plotSelection ) {
      if (!set) {
        return;
      }
      
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
  
  @Override
  protected void updateDriver(final DataStore ds, boolean freqSpace) {
    
    final boolean freqSpaceImmutable = freqSpace;
    
    displayInfoMessage("Calculating data...");
    
    SwingWorker<Integer, Void> worker = new SwingWorker<Integer, Void>() {
      @Override
      public Integer doInBackground() {
        NoiseExperiment noisExp = (NoiseExperiment) expResult;
        noisExp.setFreqSpace(freqSpaceImmutable);
        expResult.setData(ds);
        
        for (int j = 0; j < 3; ++j) {
          XYSeriesCollection xysc = expResult.getData().get(j);
          
          for (int i = 0; i < NOISE_PLOT_COUNT; ++i) {
            String name = (String) xysc.getSeriesKey(i);
            Color plotColor = COLORS[i % 3];
            seriesColorMap.put(name, plotColor);
            if (i >= 3) {
              seriesDashedSet.add(name);
            }

          }
        }
        
        set = true;
        return 0;
      }

      @Override
      public void done() {
        
        displayInfoMessage("Data loaded...drawing chart");
        
        
        northChart = buildChart( expResult.getData().get(0) );
        northChart.setTitle("Self-noise (NORTH)");
        eastChart = buildChart( expResult.getData().get(1) );
        eastChart.setTitle("Self-noise (EAST)");
        vertChart = buildChart( expResult.getData().get(2) );
        vertChart.setTitle("Self-noise (VERTICAL)");

        
        
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

    };
    
    worker.execute();
    
  }
  
  @Override
  public int panelsNeeded() {
    return 9;
  }
  
}
