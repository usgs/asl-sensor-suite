package asl.sensor.gui;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics2D;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.image.BufferedImage;

import javax.swing.JComboBox;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
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
    
    
    northChart = 
        ChartFactory.createXYLineChart( expType.getName() + " (North)",
        getXTitle(), getYTitle(), null);
    eastChart = 
        ChartFactory.createXYLineChart( expType.getName() + " (East)",
        getXTitle(), getYTitle(), null);
    vertChart = 
        ChartFactory.createXYLineChart( expType.getName() + " (Vertical)",
        getXTitle(), getYTitle(), null);
    
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

  /**
   * 
   */
  private static final long serialVersionUID = -8049021432657749975L;

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
  
  @Override
  public void updateData(final DataStore ds) {
    

    if (ds.numberFullySet() < 9) {
      displayErrorMessage("INSUFFICIENT DATA LOADED");
      return;
    }
    
    boolean freqSpace = freqSpaceBox.isSelected();
    
    updateDriver(ds, freqSpace);
    // setting the new chart is enough to update the plots
    
  }
  
  @Override
  protected void updateDriver(final DataStore ds, boolean freqSpace) {
    
    final boolean freqSpaceImmutable = freqSpace;

    displayInfoMessage("Calculating data...");

    NoiseNineExperiment noisExp = (NoiseNineExperiment) expResult;
    noisExp.setFreqSpace(freqSpaceImmutable);
    expResult.setData(ds);

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
  public BufferedImage getAsImage(int width, int height) {
    
    if (!set) {
      ChartPanel cp = new ChartPanel(chart);
      BufferedImage bi =  
          new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
      cp.setSize( new Dimension(width, height) );
      Graphics2D g = bi.createGraphics();
      cp.printAll(g);
      g.dispose();
      return bi;
      
    }
    
    height = (height * 3) / 3;
    
    // Dimension outSize = new Dimension(width, height);
    Dimension chartSize = new Dimension(width, height / 3);
    
    ChartPanel outCPanel = new ChartPanel(northChart);
    outCPanel.setSize(chartSize);
    outCPanel.setPreferredSize(chartSize);
    outCPanel.setMinimumSize(chartSize);
    outCPanel.setMaximumSize(chartSize);
    
    ChartPanel outCPanel2 = new ChartPanel(eastChart);
    outCPanel2.setSize(chartSize);
    outCPanel2.setPreferredSize(chartSize);
    outCPanel2.setMinimumSize(chartSize);
    outCPanel2.setMaximumSize(chartSize);
    
    ChartPanel outCPanel3 = new ChartPanel(vertChart);
    outCPanel3.setSize(chartSize);
    outCPanel3.setPreferredSize(chartSize);
    outCPanel3.setMinimumSize(chartSize);
    outCPanel3.setMaximumSize(chartSize);
    
    int totalHeight = 
        outCPanel.getHeight() + outCPanel2.getHeight() + outCPanel3.getHeight();
    
    BufferedImage bi = new BufferedImage(
        (int) outCPanel.getWidth(), 
        totalHeight, 
        BufferedImage.TYPE_INT_ARGB);
    
    BufferedImage northBuff = new BufferedImage(
        (int) outCPanel.getWidth(),
        (int) outCPanel.getHeight(),
        BufferedImage.TYPE_INT_ARGB);
    
    Graphics2D g = northBuff.createGraphics();
    outCPanel.printAll(g);
    g.dispose();
    
    BufferedImage eastBuff = new BufferedImage(
        (int) outCPanel2.getWidth(),
        (int) outCPanel2.getHeight(),
        BufferedImage.TYPE_INT_ARGB);

    g = eastBuff.createGraphics();
    outCPanel2.printAll(g);
    g.dispose();
    
    BufferedImage vertBuff = new BufferedImage(
        (int) outCPanel3.getWidth(),
        (int) outCPanel3.getHeight(),
        BufferedImage.TYPE_INT_ARGB);

    g = vertBuff.createGraphics();
    outCPanel3.printAll(g);
    g.dispose();
    
    int vertHeight = northBuff.getHeight() + eastBuff.getHeight();
    
    g = bi.createGraphics();
    g.drawImage(northBuff, null, 0, 0);
    g.drawImage( eastBuff, null, 0, northBuff.getHeight() );
    g.drawImage(vertBuff, null, 0, vertHeight);
    g.dispose();
    
    return bi;
  }
  
  @Override
  public int panelsNeeded() {
    return 9;
  }
  
}
