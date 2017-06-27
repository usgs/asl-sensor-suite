package asl.sensor.gui;

import java.awt.Color;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Set;
import java.util.TimeZone;

import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JDialog;
import javax.swing.JOptionPane;

import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.LogarithmicAxis;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.data.xy.XYSeriesCollection;

import asl.sensor.experiment.ExperimentEnum;
import asl.sensor.experiment.ResponseExperiment;
import asl.sensor.input.DataStore;
import asl.sensor.input.InstrumentResponse;

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

  private final static Color[] COLOR_LIST = 
      new Color[]{Color.RED, Color.BLUE, Color.GREEN};
  // values
  public static final String MAGNITUDE = ResponseExperiment.MAGNITUDE;
  public static final String ARGUMENT = ResponseExperiment.ARGUMENT;
  /**
   * 
   */
  private static final long serialVersionUID = 1L;
  
  private ValueAxis freqAxis, degreeAxis;
  
  private JCheckBox freqSpaceBox;
  private JComboBox<String> plotSelection;
  
  private JButton copyEmbedResp;
  
  private JFreeChart magChart, argChart;

  public ResponsePanel(ExperimentEnum exp) {
    super(exp);
    
    for (int i = 0; i < 3; ++i) {
      channelType[i] = "Response data (SEED data not used)";
    }
    
    String xAxisTitle = "Period (s)";
    String freqAxisTitle = "Frequency (Hz)";
    String yAxisTitle = "10 * log10( RESP(f) )";
    String degreeAxisTitle = "phi(RESP(f))";
    
    xAxis = new LogarithmicAxis(xAxisTitle);
    freqAxis = new LogarithmicAxis(freqAxisTitle);
    xAxis.setAutoRange(true);
    freqAxis.setAutoRange(true);
    
    yAxis = new NumberAxis(yAxisTitle);
    yAxis.setAutoRange(true);
    
    degreeAxis = new NumberAxis(degreeAxisTitle);
    degreeAxis.setAutoRange(true);
    ((NumberAxis) degreeAxis).setAutoRangeIncludesZero(false);
    
    ( (NumberAxis) yAxis).setAutoRangeIncludesZero(false);
    Font bold = xAxis.getLabelFont().deriveFont(Font.BOLD);
    xAxis.setLabelFont(bold);
    yAxis.setLabelFont(bold);
    freqAxis.setLabelFont(bold);
    degreeAxis.setLabelFont(bold);
    
    freqSpaceBox = new JCheckBox("Use Hz units (requires regen)");
    freqSpaceBox.setSelected(true);
    
    copyEmbedResp = new JButton("Extract an embedded response for editing");
    copyEmbedResp.addActionListener(this);
    
    plotSelection = new JComboBox<String>();
    plotSelection.addItem(MAGNITUDE);
    plotSelection.addItem(ARGUMENT);
    plotSelection.addActionListener(this);
    
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
    
    gbc.gridy += 1;
    gbc.weighty = 0;
    this.add(copyEmbedResp, gbc);
    
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
    this.add(plotSelection, gbc);
    
  }
  
  @Override
  public void actionPerformed(ActionEvent e) {
    
    super.actionPerformed(e);
    
    if ( e.getSource() == plotSelection ) {
      if (!set) {
        return;
      }
      
      JFreeChart[] charts = new JFreeChart[]{magChart, argChart};
      int idx = plotSelection.getSelectedIndex();
      chartPanel.setChart(charts[idx]);
      
      return;
      
    }
    
    if ( e.getSource() == copyEmbedResp ) {
      Set<String> respFilenames = InstrumentResponse.parseInstrumentList();
      
      List<String> names = new ArrayList<String>(respFilenames);
      Collections.sort(names);
      
      JDialog dialog = new JDialog();
      Object result = JOptionPane.showInputDialog(
          dialog,
          "Select a response to copy:",
          "RESP File Selection",
          JOptionPane.PLAIN_MESSAGE,
          null, names.toArray(),
          names.get( names.size() - 1 ) );
      
      String resultStr = (String) result;
      
      try {
        // copy response file out of embedded set and into responses folder
        
        File respDir = new File("responses/");
        
        if ( !respDir.exists() ) {
          respDir.mkdir();
        }
        
        ClassLoader cl = InputPanel.class.getClassLoader();
        InputStream is = cl.getResourceAsStream(resultStr);
        BufferedReader fr = new BufferedReader( new InputStreamReader(is) );
        StringBuilder sb = new StringBuilder();
        String line = fr.readLine();
        while ( line != null ) {
          sb.append(line);
          sb.append("\n");
          line = fr.readLine();
        }
        fr.close();
        
        StringBuilder fileNameOut = new StringBuilder();
        fileNameOut.append( respDir.getCanonicalPath() );
        fileNameOut.append("/");
        fileNameOut.append(resultStr);
        fileNameOut.append("_");
        
        SimpleDateFormat sdf = new SimpleDateFormat("YYYY.DD");
        Calendar cCal = Calendar.getInstance( sdf.getTimeZone() );
        String date = sdf.format( cCal.getTime() );
        fileNameOut.append(date);
        
        File respOut = new File( respDir.getCanonicalPath() + "/" + resultStr );
        FileWriter respWriter = new FileWriter(respOut, false);
        respWriter.write( sb.toString() ); 
        respWriter.close();
        
      } catch (IOException e1) {
        e1.printStackTrace();
      }
      
      return; // end of 'if response button clicked'
    }
    
  }
  
  protected void drawCharts() {

    plotSelection.setSelectedIndex(0);
    chart = magChart;
    chartPanel.setChart(chart);
    chartPanel.setMouseZoomable(true);
    
  }
  
  @Override
  public String[] getAdditionalReportPages() {
    ResponseExperiment respExp = (ResponseExperiment) expResult;
    InstrumentResponse[] irs = respExp.getResponses();
    String[] pages = new String[irs.length];
    for (int i = 0; i < pages.length; ++i) {
      pages[i] = irs[i].toString();
    }
    return pages;
  }
  
  @Override
  public JFreeChart[] getCharts() {
    return new JFreeChart[]{magChart, argChart};
  }

  @Override
  /**
   * Produce the filename of the report generated from this experiment.
   * Since response data is not directly associated with data at a given
   * time, rather than a sensor as a whole, we merely use the current date
   * and the first response used in the experiment.
   * @return String that will be default filename of PDF generated from data
   */
  public String getPDFFilename() {
    
    SimpleDateFormat sdf = new SimpleDateFormat("YYYY.DDD");
    sdf.setTimeZone( TimeZone.getTimeZone("UTC") );
    Calendar cCal = Calendar.getInstance( sdf.getTimeZone() );
    // experiment has no time metadata to be associated with it, get time now
    String date = sdf.format( cCal.getTime() );
    
    String test = expType.getName().replace(' ', '_');
    
    int idx = getIndexOfMainData(); // first resp in list
    String name = expResult.getInputNames().get(idx);
    
    StringBuilder sb = new StringBuilder();
    sb.append(test);
    sb.append('_');
    sb.append(name);
    sb.append('_');
    sb.append(date);
    sb.append(".pdf");
    return sb.toString();
    
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
    
    ValueAxis[] axes = new ValueAxis[]{yAxis, degreeAxis};
    
    if ( null == plotSelection ) {
      return yAxis;
    }
    
    return axes[plotSelection.getSelectedIndex()];
  }
  
  @Override
  public int panelsNeeded() {
    return 3;
  }
  
  @Override
  public int plotsToShow() {
    return 0;
  }
  
  @Override
  protected void updateData(DataStore ds) {

    set = true;
    
    seriesColorMap = new HashMap<String, Color>();
    
    boolean freqSpace = freqSpaceBox.isSelected();
    ResponseExperiment respExp = (ResponseExperiment) expResult;
    respExp.setFreqSpace(freqSpace);
    expResult.runExperimentOnData(ds);
    
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
    
    argChart = buildChart(argSeries, xAxis, degreeAxis);
    argChart.getXYPlot().getRangeAxis().setAutoRange(true);
    magChart = buildChart(magSeries, xAxis, yAxis);
    magChart.getXYPlot().getRangeAxis().setAutoRange(true);
  }

}
