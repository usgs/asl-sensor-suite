package asl.sensor.gui;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.TimeZone;

import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JPanel;
import javax.swing.filechooser.FileNameExtensionFilter;

import org.apache.pdfbox.pdmodel.PDDocument;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.annotations.XYTitleAnnotation;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.chart.title.TextTitle;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.RectangleAnchor;

import asl.sensor.experiment.Experiment;
import asl.sensor.experiment.ExperimentEnum;
import asl.sensor.experiment.ExperimentFactory;
import asl.sensor.input.DataStore;
import asl.sensor.utils.ReportingUtils;

/**
 * Panel used to display the data produced from a specified sensor test.
 * This primarily exists as a chartpanel, plus a filechooser and button used
 * in saving the chart held in the panel to a file.
 * 
 * A default construction of the GUI components exists in this class, but
 * implementing methods are suggested to override this in their constructor
 * (see existing classes such as GainPanel and NoisePanel for examples).
 * 
 * @author akearns
 *
 */
public abstract class ExperimentPanel extends JPanel implements ActionListener {
 
  private static final long serialVersionUID = -5591522915365766604L;

  /**
   * Append text to a chart's title (used for distinguishing random cal types).
   * @param chart Chart whose title will be modified
   * @param appendText Text to append to chart's current title
   */
  public static void appendChartTitle(JFreeChart chart, String appendText) {
    String titleText = chart.getTitle().getText();
    chart.getTitle().setText(titleText + appendText);
  }
  
  /**
   * Get start and end times of data for experiments that use time series data
   * @param expResult experiment with data already added
   * @return string representing the start and end of the 
   * experiment's data range
   */
  public static String getTimeStampString(Experiment expResult) {
    StringBuilder sb = new StringBuilder();
    SimpleDateFormat sdf = new SimpleDateFormat("Y.DDD.HH:mm:ss");
    sdf.setTimeZone( TimeZone.getTimeZone("UTC") );
    
    Calendar cCal = Calendar.getInstance( sdf.getTimeZone() );
    
    sb.append("Time of report generation:\n");
    sb.append( sdf.format( cCal.getTime() ) );
    sb.append('\n');
    
    long startTime = expResult.getStart();
    long endTime = expResult.getEnd();
    if ( !(startTime == 0L && endTime == 0L) ) {
      cCal.setTimeInMillis( startTime / 1000 );
      
      sb.append("Data start time:\n");
      sb.append( sdf.format( cCal.getTime() ) );
      sb.append('\n');
      
      cCal.setTimeInMillis( endTime / 1000 );
      
      sb.append("Data end time:\n");
      sb.append( sdf.format( cCal.getTime() ) );
      sb.append('\n');
    }
    return sb.toString();
  }
  protected JButton save;
  
  protected JFreeChart chart; // replace with plot object
  
  protected ChartPanel chartPanel;
  
  protected JFileChooser fc; // save image when image save button clicked
  
  public final ExperimentEnum expType; 
          // used to define experiment of each plot object
  
  protected Experiment expResult;
  
  protected ValueAxis xAxis, yAxis;
  
  public String[] channelType;
  
  protected boolean set;
  protected String[] plotTheseInBold; // given in the implementing function
  
  // this is a String because bolded names are intended to be fixed
  protected Map<String, Color> seriesColorMap;
  
  protected Set<String> seriesDashedSet;
  // these are map/set because they are based on the data read in, not fixed
  
  /**
   * Construct a new panel, using a backend defined by the passed-in enum
   * @param exp Experiment enum with corresponding backend for factory
   * instantiation
   */
  public ExperimentPanel(ExperimentEnum exp) {
    
    set = false;
    
    channelType = new String[DataStore.FILE_COUNT];
    
    // default initialization for channel type string
    for (int i = 0; i < channelType.length; ++i) {
      channelType[i] = "NOT USED";
    }
    
    seriesColorMap = new HashMap<String, Color>();
    seriesDashedSet = new HashSet<String>();
    plotTheseInBold = new String[]{};
    
    expType = exp;
    expResult = ExperimentFactory.createExperiment(exp);
    
    chart = ChartFactory.createXYLineChart( expType.getName(), 
        "", "", null);
    chartPanel = new ChartPanel(chart);
    // chartPanel.setMouseZoomable(false);
    
    fc = new JFileChooser();
    
    save = new JButton("Save plot (PNG)");
    save.addActionListener(this);
    
    // basic layout for components (recommended to override in concrete class)
    this.setLayout( new BoxLayout(this, BoxLayout.Y_AXIS) );
    this.add(chartPanel);
    this.add(save);
  }
  
  /**
   * Handle's saving this plot's chart to file (PNG image) 
   * when the save button is clicked.
   */
  @Override
  public void actionPerformed(ActionEvent e) {
    
    if( e.getSource() == save ) {
      String ext = ".png";
      fc.addChoosableFileFilter(
          new FileNameExtensionFilter("PNG image (.png)",ext) );
      fc.setFileFilter(fc.getChoosableFileFilters()[1]);
      int returnVal = fc.showSaveDialog(save);
      if (returnVal == JFileChooser.APPROVE_OPTION) {
        File selFile = fc.getSelectedFile();
        if( !selFile.getName().endsWith( ext.toLowerCase() ) ) {
          selFile = new File( selFile.toString() + ext);
        }
        try {
          ChartUtilities.saveChartAsPNG(selFile,chart,640,480);
        } catch (IOException e1) {
          e1.printStackTrace();
        }
      }
    }
  }
  
  /**
   * Gets the axes to be used to plot the data 
   */
  protected void applyAxesToChart() {
    XYPlot xyp = chart.getXYPlot();
    xyp.setDomainAxis( getXAxis() );
    xyp.setRangeAxis( getYAxis() );
  }
  
  /**
   * Function to construct a chart from the XYSeriesCollection produced
   * from this panel's backend. Any data that requires a specific plot color,
   * dashed line, or bold line have their corresponding properties applied
   * @param xyDataset Data to be plotted
   * @return XY Line Chart with the corresponding data in it
   */
  public JFreeChart buildChart(XYSeriesCollection xyDataset) {
    return buildChart( xyDataset, getXAxis(), getYAxis() );
  }
  
  /**
   * Function to construct a chart from the XYSeriesCollection produced
   * from this panel's backend. Any data that requires a specific plot color,
   * dashed line, or bold line have their corresponding properties applied.
   * This function should be used for charts in cases where a chartPanel uses
   * a combo box or similar menu item to select one of multiple charts where
   * each chart may use a different axis for either domain or range.
   * For example, step calibration has a panel with time series data of the step
   * using x-axis of seconds and y of counts, and has two other plots with
   * axes matching the charts of the response panel (x is frequency and 
   * y is magnitude and phase).
   * @param xyDataset Data to be plotted
   * @param x X-axis to be applied to the chart
   * @param y Y-axis to be applied to the chart
   * @return XY Line Chart with the corresponding data in it
   */
  public JFreeChart 
  buildChart(XYSeriesCollection xyDataset, ValueAxis x, ValueAxis y) {
    
    JFreeChart chart = ChartFactory.createXYLineChart(
        expType.getName(),
        x.getLabel(),
        y.getLabel(),
        xyDataset,
        PlotOrientation.VERTICAL,
        true, // include legend
        false, 
        false);

    // apply effects to the components that require it (i.e., NLNM time series)
    XYPlot xyPlot = chart.getXYPlot();
    XYItemRenderer xyir = xyPlot.getRenderer();


    // force certain colors and whether or not a line should be dashed

    for ( String series : seriesColorMap.keySet() ) {
      int seriesIdx = xyDataset.getSeriesIndex(series);
      if (seriesIdx >= 0) {
        xyir.setSeriesPaint( seriesIdx, seriesColorMap.get(series) );
      } else {
        continue;
      }

      if ( seriesDashedSet.contains(series) ) {
        xyir.setSeriesPaint( seriesIdx, seriesColorMap.get(series).darker() );

        BasicStroke stroke = (BasicStroke) xyir.getSeriesStroke(seriesIdx);
        if (stroke == null) {
          stroke = (BasicStroke) xyir.getBaseStroke();
        }
        float width = stroke.getLineWidth();
        int join = stroke.getLineJoin();
        int cap = stroke.getEndCap();

        float[] dashing = new float[]{1,4};

        stroke = new BasicStroke(width, cap, join, 10f, dashing, 0f);
        xyir.setSeriesStroke(seriesIdx, stroke);
      }

    }

    if ( !(plotTheseInBold.length == 0) ) {
      for (String series : plotTheseInBold) {
        int seriesIdx = xyDataset.getSeriesIndex(series);
        if (seriesIdx < 0) {
          continue;
        }

        BasicStroke stroke = (BasicStroke) xyir.getSeriesStroke(seriesIdx);
        if (stroke == null) {
          stroke = (BasicStroke) xyir.getBaseStroke();
        }
        stroke = new BasicStroke( stroke.getLineWidth()*2 );
        xyir.setSeriesStroke(seriesIdx, stroke);
        xyir.setSeriesPaint(seriesIdx, new Color(0,0,0) );
      }
    }
    
    xyPlot.setDomainAxis(x);
    xyPlot.setRangeAxis(y);
    
    return chart;
  }
  
  /**
   * Clear chart data and display text that it is loading new data
   */
  protected void clearChartAndSetProgressData() {
    chart = 
        ChartFactory.createXYLineChart( expType.getName(), "",  "",  null );
    applyAxesToChart();
    chartPanel.setChart(chart);
    displayInfoMessage("Running calculation...");
  }
  
  /**
   * Overlay an error message in the event of an exception or other issue
   * @param errMsg Text of the message to be displayed
   */
  public void displayErrorMessage(String errMsg) {
    XYPlot xyp = (XYPlot) chartPanel.getChart().getPlot();
    TextTitle result = new TextTitle();
    result.setText(errMsg);
    result.setBackgroundPaint(Color.red);
    result.setPaint(Color.white);
    XYTitleAnnotation xyt = new XYTitleAnnotation(0.5, 0.5, result,
        RectangleAnchor.CENTER);
    xyp.clearAnnotations();
    xyp.addAnnotation(xyt);
  }

  /**
   * Overlay informational text, such as extra results and statistics for plots
   * @param infoMsg
   */
  public void displayInfoMessage(String infoMsg) {
    XYPlot xyp = (XYPlot) chartPanel.getChart().getPlot();
    TextTitle result = new TextTitle();
    result.setText(infoMsg);
    XYTitleAnnotation xyt = new XYTitleAnnotation(0.5, 0.5, result,
        RectangleAnchor.CENTER);
    xyp.clearAnnotations();
    xyp.addAnnotation(xyt);
  }
  
  /**
   * Used to return more detailed information from the experiment, such
   * as a full-page report of a best-fit response file. Most experiments
   * will not need to override this method, but it may be useful to add
   * more detailed or verbose information that cannot be fit into a single
   * report page.
   * @return Array of strings, each one to be written to a new report page
   */
  public String[] getAdditionalReportPages() {
    return new String[]{};
  }

  /**
   * Return image of panel's plots with specified dimensions
   * Used to compile PNG image of all charts contained in this panel
   * @param width Width of output image in pixels
   * @param height Height of output image in pixels
   * @return buffered image of this panel's chart
   */
  public BufferedImage getAsImage(int width, int height) {
    
    JFreeChart[] jfcs = getCharts();
    return ReportingUtils.chartsToImage(width, height, jfcs);
    
  }


  /**
   * Returns the identifiers of each input plot being used, such as 
   * "calibration input" for the calibration tests.
   * @return Strings used to populate channel type identifiers in input panel
   */
  public String[] getChannelTypes() {
    return channelType;
  }
  
  /**
   * Return all chart panels used in this object;
   * to be overridden by implementing experiment panels that contain multiple
   * charts.
   * Primary use of this function is to enumerate charts to save as images/PDF
   * @return All chartpanels used in this object
   */
  public JFreeChart[] getCharts() {
    return new JFreeChart[]{chart};
  }
  
  /**
   * Used to return any title insets as text format for saving in PDF,
   * to be overridden by any panel that uses an inset
   * @return String with any relevant parameters in it
   */
  public String getInsetString() {
    return "";
  }

  /**
   * Used to return any metadata from the experiment to be saved in PDF
   * to be overridden by panels with data that should be included in the report
   * that it would not make sense to display in the inset, such as response
   * filenames
   * @return A string with any additional data to be included in the PDF report
   */
  public String getMetadataString() {
    return "";
  }
  
  /**
   * For report generation, give the list of indices of response files used
   * in the plot
   * @return list where each index is a relevant response file
   */
  public int[] getResponseIndices() {
    return expResult.listActiveResponseIndices();
  }
  
  /**
   * Default x-axis return function.
   * Though the x-axis is a local variable, some panels may have multiple unit
   * types for the x-axis (i.e., for units of seconds vs. Hz); accessing
   * the x-axis object through this function allows for overrides allowing for
   * more flexibility.
   * @return ValueAxis to be applied to chart
   */
  public ValueAxis getXAxis() {
    return xAxis;
  }
  
  /**
   * Default y-axis return function. As with getXAxis, designed to be overridden
   * for charts that may use multiple scales.
   * @return ValueAxis to be applied to chart
   */
  public ValueAxis getYAxis() {
    return yAxis;
  }
  
  /**
   * Function used to query backend on whether or not a datastore has all the
   * data that a backend needs to calculate. This is used mainly to inform
   * the main window (see SensorSuite class) that the generate result button
   * can be set active
   * @param ds Datastore to run data check on
   * @return True if the backend can run with the data provided
   */
  public boolean hasEnoughData(final DataStore ds) {
    return expResult.hasEnoughData(ds);
  }
  
  /**
   * True if data has been loaded into the experiment backend yet
   * @return True if there is data to process
   */
  public boolean hasRun() {
    return set;
  }
  
  /**
   * Function template for informing main window of number of panels to display
   * to fit all data needed by the program
   * @return Number of plots to show in the input panel
   */
  public abstract int panelsNeeded();
  
  /**
   * Number of panels to return in an output report
   * @return number of panels to include 
   */
  public int plotsToShow() {
    return panelsNeeded();
  }
  
  /**
   * Takes a PDF and adds a page dedicated to string data related to the
   * data that has been passed into the experiment backend,
   * including any text that might be included in chart title insets and
   * the input data start and end timestamps
   * @param pdf PDF document to append data to
   */
  public void saveInsetDataText(PDDocument pdf) {
    
    StringBuilder sb = new StringBuilder( getInsetString() );
    if ( sb.length() > 0 ) {
      sb.append("\n\n");
    }
    String metadata = getMetadataString();
    if ( metadata.length() > 0 ) {
      sb.append(metadata);
      sb.append("\n \n");
    }
    sb.append( getTimeStampString(expResult) );
    ReportingUtils.textToPDFPage( sb.toString(), pdf );
    ReportingUtils.textListToPDFPages( pdf, getAdditionalReportPages() );
    return;
  }
  
  /**
   * Loads in charts used in this panel and prints them out in a PDF document
   * and includes any relevant metadata / analysis results as plain text
   * @param pdf Document to save data to
   */
  public void savePDFResults(PDDocument pdf) {

    int width = 1280;
    int height = 960;
    JFreeChart[] cps = getCharts();
    
    ReportingUtils.chartsToPDFPage(width, height, pdf, cps);

    saveInsetDataText(pdf);
    
  }

  /**
   * Used to plot the results of a backend function from an experiment
   * using a collection of XYSeries mapped by strings. This will be set to
   * the default chart object held by the panel.
   * @param xyDataset collection of XYSeries to plot
   */
  protected void setChart(XYSeriesCollection xyDataset) {

     chart = buildChart(xyDataset);

  }
  
  /**
   * Function template for sending input to a backend fucntion and collecting
   * the corresponding data
   * @param ds DataStore object containing seed and resp files
   */
  public abstract void updateData(final DataStore ds);
  // details of how to run updateData are left up to the implementing panel
  // however, the boolean "set" should be set to true to enable PDF saving

}
