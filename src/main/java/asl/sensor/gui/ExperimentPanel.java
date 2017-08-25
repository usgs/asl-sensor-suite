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
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TimeZone;

import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JPanel;
import javax.swing.SwingWorker;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.filechooser.FileNameExtensionFilter;

import org.apache.pdfbox.pdmodel.PDDocument;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.annotations.XYTitleAnnotation;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.SeriesRenderingOrder;
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
 * This class also includes a number of utility functions used in report
 * generation. This is done with the expectation that some (but not necessarily
 * all) such functions will be overridden by an implementing class, such as
 * functions for creating reports with multiple charts and ensuring that any
 * and all text data relevant to the panel is included in a report. In addition
 * there are methods to display on the panel the current status of the
 * experiment based on calls to a status change report in the backend.
 * 
 * @author akearns
 *
 */
public abstract class ExperimentPanel 
extends JPanel 
implements ActionListener, ChangeListener {
 
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
      cCal.setTimeInMillis(startTime);
      
      sb.append("Data start time:\n");
      sb.append( sdf.format( cCal.getTime() ) );
      sb.append('\n');
      
      cCal.setTimeInMillis(endTime);
      
      sb.append("Data end time:\n");
      sb.append( sdf.format( cCal.getTime() ) );
      sb.append('\n');
    }
    return sb.toString();
  }
  
  /**
   * Reverses an xyplot rendering order, allowing curves that would otherwise be
   * at the background are inverted and placed in the foreground instead.
   * That is, if a curve is to be rendered behind a different curve, it will be
   * rendered instead with that series in front of the other curve.
   * @param chart Chart with plot rendering order to be reversed. Must use
   * an XY plot (i.e., is an XYLineSeries chart)
   */
  public static void invertSeriesRenderingOrder(JFreeChart chart) {
    XYPlot xyp = chart.getXYPlot();
    SeriesRenderingOrder sro = xyp.getSeriesRenderingOrder();
    if ( sro.equals(SeriesRenderingOrder.FORWARD) ) {
      xyp.setSeriesRenderingOrder(SeriesRenderingOrder.REVERSE);
    } else {
      xyp.setSeriesRenderingOrder(SeriesRenderingOrder.FORWARD);
    }
    
  }
  
  protected JButton save;
  
  protected JFreeChart chart; // the chart shown in the panel
  
  protected ChartPanel chartPanel; // component used to hold the shown chart
  // (if an experiment has multiple charts to show, ideally each should be
  // selectable through some sort of menu with the active menu option deciding
  // which chart should be displayed in this panel)
  
  protected JFileChooser fc; // save image when image save button clicked
  
  public final ExperimentEnum expType; 
    // used to define experiment of each plot object (i.e., chart name)
  
  protected Experiment expResult;
    // experiment actually being run (call its 'setData' method to run backend)
    // experiments use builder pattern -- set necessary variables like
    // angle offset or x-axis units before running the experiment
  
  protected ValueAxis xAxis, yAxis;
    // default axes to use with the default chart
  
  public String[] channelType;
    // used to give details in input panel about what users needs to load where
  protected boolean set; // true if the experiment has run
  
  protected String[] plotTheseInBold; // given in the implementing function
  // this is a String because bolded names are intended to be fixed
  // (i.e., NLNM, NHNM, not dependent on user input)
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
    expResult.addChangeListener(this);
    
    chart = ChartFactory.createXYLineChart( expType.getName(), 
        "", "", null);
    chartPanel = new ChartPanel(chart);
    // chartPanel.setMouseZoomable(false);
    
    fc = new JFileChooser();
    
    save = new JButton("Save plot (PNG)");
    save.addActionListener(this);
    
    // basic layout for components (recommended to override in concrete class)
    // if specific formatting or additional components are unnecessary, the
    // implementing class can simply call super(expType) to make a panel
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
    
    if ( e.getSource() == save ) {
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
   * Clear chart, for when calculation operations have been cancelled
   */
  public void clearChart() {
    set = false;
    chart = 
        ChartFactory.createXYLineChart( 
            expType.getName(), 
            getXAxis().getLabel(),  
            getYAxis().getLabel(),  null );
    chartPanel.setChart(chart);
  }
  
  /**
   * Clear chart data and display text that it is loading new data
   */
  protected void clearChartAndSetProgressData() {
    clearChart();
    displayInfoMessage("Running calculation...");
  }
  
  /**
   * Overlay an error message in the event of an exception or other issue
   * @param errMsg Text of the message to be displayed
   */
  public void displayErrorMessage(String errMsg) {
    clearChart();
    XYPlot xyp = (XYPlot) chart.getPlot();
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
    result.setBackgroundPaint(Color.white);
    XYTitleAnnotation xyt = new XYTitleAnnotation(0.5, 0.5, result,
        RectangleAnchor.CENTER);
    xyp.clearAnnotations();
    xyp.addAnnotation(xyt);
  }
  
  /**
   * Function template for applying data taken from updateData function
   * and drawing the given charts on the screen as necessary
   */
  protected abstract void drawCharts();

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
   * Produce all data used in PDF reports as a single string that can be
   * written to a text file
   * @return Data string including all metadata and relevant infrom from
   * an experiment
   */
  public String getAllTextData() {
    StringBuilder sb = new StringBuilder( getInsetStrings() );
    if ( sb.length() > 0 ) {
      sb.append("\n\n");
    }
    String metadata = getMetadataString();
    if ( metadata.length() > 0 ) {
      sb.append(metadata);
      sb.append("\n\n");
    }
    sb.append( getTimeStampString(expResult) );
    sb.append("\n\n");
    String[] extraText = getAdditionalReportPages();
    for (String text : extraText) {
      sb.append(text);
      sb.append("\n\n");
    }
    return sb.toString();
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
   * Get index for station name for most relevant data. For example, in the
   * case of a calibration, this is usually the second input. In most cases,
   * however, the first data input will be sufficient. Used in report filename
   * generation.
   * @return Index of data used to get naem for report
   */
  protected int getIndexOfMainData() {
    return 0;
  }
  
  /**
   * Used to return any title insets as text format for saving in PDF,
   * to be overridden by any panel that uses an inset
   * @return String with any relevant parameters in it
   */
  public String getInsetStrings() {
    return "";
  }
  
  /**
   * Used to return any metadata from the experiment to be saved in PDF
   * to be overridden by panels with data that should be included in the report
   * that it would not make sense to display in the inset, such as filenames.
   * If there is no metadata besides filenames of input data, that is returned
   * by the function, but experiments may need to augment with additional data.
   * @return A string with any additional data to be included in the PDF report
   */
  public String getMetadataString() {
    List<String> names = expResult.getInputNames();
    StringBuilder sb = new StringBuilder("Input filenames, ");
    sb.append(" with SEED and RESP files paired as appropriate:\n");
    for (String name : names) {
      sb.append(name);
      sb.append('\n');
    }
    return sb.toString();
  }
  
  /**
   * Produce the filename of the report generated from this experiment.
   * Has the format TEST_STATION_YEAR.DAY unless overridden
   * @return String that will be default filename of PDF generated from data
   */
  public String getPDFFilename() {
    
    SimpleDateFormat sdf = new SimpleDateFormat("YYYY.DDD");
    sdf.setTimeZone( TimeZone.getTimeZone("UTC") );
    
    String date;
    long time = expResult.getStart();
    if (time > 0) {
      date = sdf.format(time);
    } else {
      Calendar cCal = Calendar.getInstance( sdf.getTimeZone() );
      date = sdf.format( cCal.getTime() );
    }
    
    // turn spaces into underscores
    String test = expType.getName().replace(' ', '_'); // name of experiment
    // make sure parentheses in filenames aren't causing issues
    // (i.e., better than dealing with system-specific setups to escape them)
    test = test.replace('(','_');
    test = test.replace(')','_');
    
    int idx = getIndexOfMainData();
    String name = expResult.getInputNames().get(idx); // name of input data
    
    StringBuilder sb = new StringBuilder();
    sb.append(test);
    sb.append('_');
    sb.append(name);
    sb.append('_');
    sb.append(date);
    sb.append(".pdf");
    return sb.toString();
    
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
   * Function to be overridden by implementing class that will add an extra
   * page to PDF reports including charts with less-essential data, such as
   * the plots of residual values over a range for parameter-fitting
   * @return List of charts to show on a second page of PDF reports
   */
  public JFreeChart[] getSecondPageCharts() {
    return new JFreeChart[]{};
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
   * Function to call to run experiment backend on specific data, using the
   * given swingworker
   * @param ds Data to evaluate the backend on
   * @param worker Worker thread to run the backend in, presumably the 
   * worker object originating in the main class for the suite
   */
  public SwingWorker<Boolean, Void> 
  runExperiment(final DataStore ds, SwingWorker<Boolean, Void> worker) {
    
    
    worker.execute();
    
    return worker;
    
  }

  /**
   * Takes a PDF and adds a page dedicated to string data related to the
   * data that has been passed into the experiment backend,
   * including any text that might be included in chart title insets and
   * the input data start and end timestamps
   * @param pdf PDF document to append data to
   */
  public void saveInsetDataText(PDDocument pdf) {
    
    StringBuilder sb = new StringBuilder( getInsetStrings() );
    if ( sb.length() > 0 ) {
      sb.append("\n \n");
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
    JFreeChart[] charts = getCharts();
    
    ReportingUtils.chartsToPDFPage(width, height, pdf, charts);
    JFreeChart[] page2 = getSecondPageCharts();
    if (page2.length > 0) {
      ReportingUtils.chartsToPDFPage(width, height, pdf, page2);
    }
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
   * Used to identify completion of an experiment to a containing thread
   */
  public void setDone() {
    firePropertyChange("Backend completed", false, set);
    drawCharts();
  }

  @Override
  /**
   * Used to print out status of the experiment backend onto the chart when
   * the backend status changes
   */
  public void stateChanged(ChangeEvent e) {
    if ( e.getSource() == expResult ) {
      String info = expResult.getStatus();
      displayInfoMessage(info);
    }
  }
  
  /**
   * Function template for sending input to a backend fucntion and collecting
   * the corresponding data
   * @param ds DataStore object containing seed and resp files
   */
  protected abstract void updateData(final DataStore ds);
  // details of how to run updateData are left up to the implementing panel
  // however, the boolean "set" should be set to true to enable PDF saving
  
}
