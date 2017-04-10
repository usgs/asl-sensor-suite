package asl.sensor.gui;

import java.awt.Color;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;

import javax.swing.JCheckBox;
import javax.swing.JPanel;

import org.jfree.chart.axis.LogarithmicAxis;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.data.xy.XYSeriesCollection;

import asl.sensor.experiment.ExperimentEnum;
import asl.sensor.experiment.NoiseExperiment;
import asl.sensor.input.DataStore;

/**
 * Panel for displaying the results of the self-noise experiment (3-input).
 * In addition to general requirements of output panels, also includes
 * a checkbox to choose between frequency and interval x-axis and
 * the variant axes for when that box is checked.
 * @author akearns
 *
 */
public class NoisePanel extends ExperimentPanel {

  /**
   * Auto-generated serialize ID
   */
  private static final long serialVersionUID = 9018553361096758354L;
  
  protected JCheckBox freqSpaceBox;
  
  protected final int NOISE_PLOT_COUNT = 6;
  // three PSDs, three self-noise calcs
  
  protected final Color[] COLORS = {Color.RED, Color.BLUE, Color.GREEN};
  
  private NumberAxis freqAxis;
  private String freqAxisTitle; // name to use when plotting in units of Hz
  
  /**
   * Constructs a new panel and lays out all the components in it
   * @param exp
   */
  public NoisePanel(ExperimentEnum exp) {
    
    // create chart, chartPanel, save button & file chooser, 
    super(exp);
    
    for (int i = 0; i < 3; ++i) {
      channelType[i] = "Input data (RESP required)";
    }
    
    plotTheseInBold = new String[]{"NLNM","NHNM"};
    
    // instantiate local fields
    xAxisTitle = "Period (s)";
    freqAxisTitle = "Frequency (Hz)";
    yAxisTitle = "Power (rel. 1 (m/s^2)^2/Hz)";
    xAxis = new LogarithmicAxis(xAxisTitle);
    freqAxis = new LogarithmicAxis(freqAxisTitle);
    yAxis = new NumberAxis(yAxisTitle);
    yAxis.setAutoRange(true);
    ( (NumberAxis) yAxis).setAutoRangeIncludesZero(false);
    Font bold = xAxis.getLabelFont().deriveFont(Font.BOLD);
    xAxis.setLabelFont(bold);
    yAxis.setLabelFont(bold);
    freqAxis.setLabelFont(bold);
    
    freqSpaceBox = new JCheckBox("Use Hz units (requires regen)");
    freqSpaceBox.setSelected(false);
    
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
    JPanel spacer = new JPanel();
    spacer.setPreferredSize( freqSpaceBox.getPreferredSize() );
    this.add(spacer, gbc);
  }
  
  @Override
  public void actionPerformed(ActionEvent e) {
    // overridden in the event we add more stuff to this panel
    super.actionPerformed(e); // only actionlistener here
  }
  
  /**
   * Gets the x-axis for this panel based on whether or not the
   * selection box to plot in units of Hz is selected. If it is, this
   * plot will have frequency units of Hz in the x-axis, otherwise it will have
   * interval units of seconds in it
   */
  @Override
  public ValueAxis getXAxis() {
    
    // true if using Hz units
    if ( freqSpaceBox.isSelected() ) {
        return freqAxis;
    }
    
    return xAxis;
    
  }
  
  /**
   * As with the getXAxis function, if the selection box to use Hz is selected,
   * this will not return the default x-axis title but instead will return one
   * reflecting the use of frequency/Hz units (the default is using the 
   * interval (time)).
   */
  @Override
  public String getXTitle() {
    if (freqSpaceBox.isSelected()) {
      return freqAxisTitle;
    }
    return xAxisTitle;
  }

  /**
   * Initially called function to calculate self-noise when data is passed in
   */
  @Override
  public void updateData(final DataStore ds) {
    
    set = true;
    
    boolean freqSpace = freqSpaceBox.isSelected();
    
    updateDriver(ds, freqSpace);
    // setting the new chart is enough to update the plots
    
    
  }
  

  /**
   * Uses a threaded call to run the self-noise calculations in the background,
   * to be plotted. In the process of getting these, the names of the data
   * are used to populate the color of the input data, as well as whether or not
   * data should be dashed or solid (the color choices are designed to
   * match the first three input plots, where the data is expected from, and
   * the dashes are used to distinguish the PSD from the self-noise plot)
   * @param ds 
   *  DataStore object that contains the seed and resp files to calculate
   * @param freqSpace Boolean matching whether or not to plot in units of
   * frequency if true (Hz) or in units of interval if false (s) 
   */
  protected void updateDriver(final DataStore ds, boolean freqSpace) {
    
    final boolean freqSpaceImmutable = freqSpace;
    
    displayInfoMessage("Calculating data...");
    
    NoiseExperiment noisExp = (NoiseExperiment) expResult;
    noisExp.setFreqSpace(freqSpaceImmutable);
    expResult.setData(ds);

    XYSeriesCollection xysc = expResult.getData().get(0);

    for (int i = 0; i < NOISE_PLOT_COUNT; ++i) {
      String name = (String) xysc.getSeriesKey(i);
      Color plotColor = COLORS[i % 3];
      seriesColorMap.put(name, plotColor);
      if (i >= 3) {
        seriesDashedSet.add(name);
      }

    }

    displayInfoMessage("Data loaded...drawing chart");

    setChart( expResult.getData().get(0) );

    chartPanel.setChart(chart);
    chartPanel.setMouseZoomable(true);

    
  }

  @Override
  public int panelsNeeded() {
    
    return 3;
  }
  

}
