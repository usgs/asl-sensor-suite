package asl.sensor.gui;

import static asl.utils.ReportingUtils.COLORS;
import static asl.utils.ReportingUtils.chartsToPDFPage;
import static asl.utils.ReportingUtils.textListToPDFPages;
import static asl.utils.ReportingUtils.textToPDFPage;

import asl.sensor.ExperimentFactory;
import asl.sensor.experiment.Experiment;
import asl.sensor.input.Configuration;
import asl.sensor.input.DataStore;
import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Component;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.Arrays;
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
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.filechooser.FileNameExtensionFilter;
import org.apache.pdfbox.pdmodel.PDDocument;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.ChartUtils;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.annotations.XYTitleAnnotation;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.SeriesRenderingOrder;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.chart.title.TextTitle;
import org.jfree.chart.ui.RectangleAnchor;
import org.jfree.data.xy.XYSeriesCollection;

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
 */
public abstract class ExperimentPanel
    extends JPanel
    implements ActionListener, ChangeListener {

  public static final ThreadLocal<SimpleDateFormat> DATE_TIME_FORMAT =
      ThreadLocal.withInitial(() -> {
        SimpleDateFormat format = new SimpleDateFormat("yyyy.DDD.HH:mm:ss.SSS");
        format.setTimeZone(TimeZone.getTimeZone("UTC"));
        return format;
      });
  static final ThreadLocal<SimpleDateFormat> DATE_FORMAT =
      ThreadLocal.withInitial(() -> {
        SimpleDateFormat format = new SimpleDateFormat("yyyy.DDD");
        format.setTimeZone(TimeZone.getTimeZone("UTC"));
        return format;
      });
  private static final long serialVersionUID = -5591522915365766604L;
  final JButton save; // easy access to saving output as png
  final ChartPanel chartPanel; // component used to hold the shown chart
  final ExperimentFactory expType;
  final String[] channelType;
  final Set<String> seriesDashedSet;
  // colorblind-friendly red, blue, and green approximates in that order
  // (if an experiment has multiple charts to show, ideally each should be
  // selectable through some sort of menu with the active menu option used to control
  // which chart should be displayed in this panel)
  private final JFileChooser fileChooser; // save image when image save button clicked
  JFreeChart chart; // the chart shown in the panel
  // used to define experiment of each plot object (i.e., chart name)
  ValueAxis xAxis;
  ValueAxis yAxis;
  // experiment actually being run (call its 'setData' method to run backend)
  // experiments use builder pattern -- set necessary variables like
  // angle offset or x-axis units before running the experiment
  // used to give details in input panel about what users needs to load where
  boolean set; // true if the experiment has run
  // default axes to use with the default chart
  Experiment expResult;
  String[] plotTheseInBold; // given in the implementing function
  // this is a String because bolded names are intended to be fixed
  // (i.e., NLNM, NHNM, not dependent on user input)
  Map<String, Color> seriesColorMap;
  /**
   * Construct a new panel, using a backend defined by the passed-in enum
   *
   * @param experiment Experiment enum with corresponding backend for factory
   * instantiation
   */
  ExperimentPanel(ExperimentFactory experiment) {
    set = false;

    channelType = new String[DataStore.FILE_COUNT];

    // default initialization for channel type string
    Arrays.fill(channelType, "NOT USED");

    seriesColorMap = new HashMap<>();
    seriesDashedSet = new HashSet<>();
    plotTheseInBold = new String[]{};

    expType = experiment;
    expResult = experiment.createExperiment();
    expResult.addChangeListener(this);

    chart = ChartFactory.createXYLineChart(expType.getName(),
        "", "", null);
    chartPanel = new ChartPanel(chart);

    fileChooser = new JFileChooser();

    save = new JButton("Save plot (PNG)");
    save.addActionListener(this);

    // basic layout for components (recommended to override in concrete class)
    // if specific formatting or additional components are unnecessary, the
    // implementing class can simply call super(expType) to make a panel
    this.setLayout(new BoxLayout(this, BoxLayout.PAGE_AXIS));
    this.add(chartPanel);
    this.add(save);
    save.setAlignmentX(Component.CENTER_ALIGNMENT);
  }

  static Color getColor(int idx) {
    if (Configuration.getInstance().useColorblindColors()) {
      return COLORS[idx % COLORS.length];
    }
    switch (idx % 3) {
      case 0:
        return Color.RED;
      case 1:
        return Color.BLUE;
      default:
        return Color.GREEN;
    }
  }

  /**
   * Append text to a chart's title (used for distinguishing random cal types).
   *
   * @param chart Chart whose title will be modified
   * @param appendText Text to append to chart's current title
   */
  static void appendChartTitle(JFreeChart chart, String appendText) {
    String titleText = chart.getTitle().getText();
    chart.getTitle().setText(titleText + appendText);
  }
  // these are map/set because they are based on the data read in, not fixed

  static TextTitle getDefaultTextTitle() {
    TextTitle result = new TextTitle();
    Font font = result.getFont();
    font = font.deriveFont(font.getSize() + 2f);
    result.setFont(font);
    result.setBackgroundPaint(Color.WHITE);
    return result;
  }

  /**
   * Reverses an xyplot rendering order, allowing curves that would otherwise be
   * at the background are inverted and placed in the foreground instead.
   * That is, if a curve is to be rendered behind a different curve, it will be
   * rendered instead with that series in front of the other curve.
   *
   * @param chart Chart with plot rendering order to be reversed. Must use
   * an XY plot (i.e., is an XYLineSeries chart)
   */
  public static void invertSeriesRenderingOrder(JFreeChart chart) {
    XYPlot plot = chart.getXYPlot();
    if (plot.getSeriesRenderingOrder().equals(SeriesRenderingOrder.FORWARD)) {
      plot.setSeriesRenderingOrder(SeriesRenderingOrder.REVERSE);
    } else {
      plot.setSeriesRenderingOrder(SeriesRenderingOrder.FORWARD);
    }

  }

  /**
   * Handle's saving this plot's chart to file (PNG image)
   * when the save button is clicked.
   */
  @Override
  public void actionPerformed(ActionEvent event) {

    if (event.getSource() == save) {
      String ext = ".png";
      fileChooser.addChoosableFileFilter(
          new FileNameExtensionFilter("PNG image (.png)", ext));
      fileChooser.setFileFilter(fileChooser.getChoosableFileFilters()[1]);
      int returnVal = fileChooser.showSaveDialog(save);
      if (returnVal == JFileChooser.APPROVE_OPTION) {
        File selFile = fileChooser.getSelectedFile();
        if (!selFile.getName().toLowerCase().endsWith(ext)) {
          selFile = new File(selFile.toString() + ext);
        }
        try {
          ChartUtils.saveChartAsPNG(selFile, chart, 640, 480);
        } catch (IOException e1) {
          e1.printStackTrace();
        }
      }
    }
  }

  /**
   * Gets the axes to be used to plot the data
   */
  void applyAxesToChart() {
    XYPlot plot = chart.getXYPlot();
    plot.setDomainAxis(getXAxis());
    plot.setRangeAxis(getYAxis());
  }

  /**
   * Function to construct a chart from the XYSeriesCollection produced
   * from this panel's backend. Any data that requires a specific plot color,
   * dashed line, or bold line have their corresponding properties applied
   *
   * @param xyDataset Data to be plotted
   * @return XY Line Chart with the corresponding data in it
   */
  public JFreeChart buildChart(XYSeriesCollection xyDataset) {
    return buildChart(xyDataset, getXAxis(), getYAxis());
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
   *
   * @param xyDataset Data to be plotted
   * @param xAxis X-axis to be applied to the chart
   * @param yAxis Y-axis to be applied to the chart
   * @return XY Line Chart with the corresponding data in it
   */
  public JFreeChart
  buildChart(XYSeriesCollection xyDataset, ValueAxis xAxis, ValueAxis yAxis) {

    JFreeChart chart = ChartFactory.createXYLineChart(
        expType.getName(),
        xAxis.getLabel(),
        yAxis.getLabel(),
        xyDataset,
        PlotOrientation.VERTICAL,
        true, // include legend
        false,
        false);

    if (xyDataset == null) {
      return chart;
    }

    // apply effects to the components that require it (i.e., NLNM time series)
    XYPlot xyPlot = chart.getXYPlot();
    XYItemRenderer renderer = xyPlot.getRenderer();

    if (seriesColorMap.size() == 0) {
      for (int seriesIndex = 0; seriesIndex < xyDataset.getSeriesCount(); ++seriesIndex) {
        renderer.setSeriesPaint(seriesIndex, getColor(seriesIndex));
      }
    }

    // now, make everything thicker!
    for (int seriesIndex = 0; seriesIndex < xyDataset.getSeriesCount(); ++seriesIndex) {
      BasicStroke stroke = (BasicStroke) renderer.getSeriesStroke(seriesIndex);
      if (stroke == null) {
        stroke = (BasicStroke) renderer.getDefaultStroke();
      }
      float widthOffset = Configuration.getInstance().getLineWidthOffset();
      float width = stroke.getLineWidth() + widthOffset;
      int join = stroke.getLineJoin();
      int cap = stroke.getEndCap();

      stroke = new BasicStroke(width, cap, join, 10f);
      renderer.setSeriesStroke(seriesIndex, stroke);
    }

    // force certain colors and whether or not a line should be dashed
    for (String series : seriesColorMap.keySet()) {
      int seriesIndex = xyDataset.getSeriesIndex(series);
      if (seriesIndex >= 0) {
        renderer.setSeriesPaint(seriesIndex, seriesColorMap.get(series));
        BasicStroke stroke = (BasicStroke) renderer.getSeriesStroke(seriesIndex);
        float width = stroke.getLineWidth();
        int join = stroke.getLineJoin();
        int cap = stroke.getEndCap();

        stroke = new BasicStroke(width, cap, join, 10f);
        renderer.setSeriesStroke(seriesIndex, stroke);
      } else {
        continue;
      }

      if (seriesDashedSet.contains(series)) {
        // applying darkening twice makes it easier to distinguish from other lines in plot
        Color darkerColor = seriesColorMap.get(series).darker().darker();
        renderer.setSeriesPaint(seriesIndex, darkerColor);

        BasicStroke stroke = (BasicStroke) renderer.getSeriesStroke(seriesIndex);
        float width = stroke.getLineWidth();
        int join = stroke.getLineJoin();
        int cap = stroke.getEndCap();

        float[] dashing = new float[]{10, 2};

        stroke = new BasicStroke(width, cap, join, 10f, dashing, 0f);
        renderer.setSeriesStroke(seriesIndex, stroke);
      }
    }

    // EXTRA THICK
    if (!(plotTheseInBold.length == 0)) {
      for (String series : plotTheseInBold) {
        int seriesIndex = xyDataset.getSeriesIndex(series);
        if (seriesIndex < 0) {
          continue;
        }

        BasicStroke stroke = (BasicStroke) renderer.getSeriesStroke(seriesIndex);
        if (stroke == null) {
          stroke = (BasicStroke) renderer.getDefaultStroke();
        }
        stroke = new BasicStroke(stroke.getLineWidth() * 2);
        renderer.setSeriesStroke(seriesIndex, stroke);
        renderer.setSeriesPaint(seriesIndex, new Color(0, 0, 0));
      }
    }
    Font bold = xAxis.getLabelFont();
    xAxis.setLabelFont(bold);
    yAxis.setLabelFont(bold);
    xyPlot.setDomainAxis(xAxis);
    xyPlot.setRangeAxis(yAxis);

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
            getYAxis().getLabel(), null);
    chartPanel.setChart(chart);
  }

  /**
   * Clear chart data and display text that it is loading new data
   */
  void clearChartAndSetProgressData() {
    clearChart();
    displayInfoMessage("Running calculation...");
  }

  /**
   * Overlay an error message in the event of an exception or other issue
   *
   * @param errMsg Text of the message to be displayed
   */
  public void displayErrorMessage(String errMsg) {
    clearChart();
    XYPlot plot = (XYPlot) chart.getPlot();
    TextTitle result = new TextTitle();
    result.setText(errMsg);
    result.setBackgroundPaint(Color.red);
    result.setPaint(Color.white);
    XYTitleAnnotation xyt = new XYTitleAnnotation(0.5, 0.5, result,
        RectangleAnchor.CENTER);
    plot.clearAnnotations();
    plot.addAnnotation(xyt);
  }

  /**
   * Overlay informational text, such as extra results and statistics for plots
   */
  void displayInfoMessage(String infoMsg) {
    XYPlot plot = (XYPlot) chartPanel.getChart().getPlot();
    TextTitle result = new TextTitle();
    result.setText(infoMsg);
    result.setBackgroundPaint(Color.white);
    XYTitleAnnotation xyt = new XYTitleAnnotation(0.5, 0.5, result,
        RectangleAnchor.CENTER);
    plot.clearAnnotations();
    plot.addAnnotation(xyt);
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
   *
   * @return Array of strings, each one to be written to a new report page
   */
  String[] getAdditionalReportPages() {
    return new String[]{};
  }


  /**
   * Produce all data used in PDF reports as a single string that can be
   * written to a text file
   *
   * @return Data string including all metadata and relevant infrom from
   * an experiment
   */
  public String getAllTextData() {
    StringBuilder sb = new StringBuilder(expResult.getReportString());
    if (sb.length() > 0) {
      sb.append("\n\n");
    }
    String metadata = getMetadataString();
    if (metadata.length() > 0) {
      sb.append(metadata);
      sb.append("\n\n");
    }
    sb.append(expResult.getFormattedDateRange());
    sb.append("\n\n");
    String[] extraText = getAdditionalReportPages();
    for (String text : extraText) {
      sb.append(text);
      sb.append("\n\n");
    }
    return sb.toString();
  }

  /**
   * Returns the identifiers of each input plot being used, such as
   * "calibration input" for the calibration tests.
   *
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
   *
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
   *
   * @return Index of data used to get naem for report
   */
  int getIndexOfMainData() {
    return 0;
  }

  /**
   * Used to return any metadata from the experiment to be saved in PDF
   * to be overridden by panels with data that should be included in the report
   * that it would not make sense to display in the inset, such as filenames.
   * If there is no metadata besides filenames of input data, that is returned
   * by the function, but experiments may need to augment with additional data.
   *
   * @return A string with any additional data to be included in the PDF report
   */
  String getMetadataString() {
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
   *
   * @return String that will be default filename of PDF generated from data
   */
  public String getPDFFilename() {

    SimpleDateFormat dateFormat = DATE_FORMAT.get();

    String date;
    long time = expResult.getStart();
    if (time > 0) {
      date = dateFormat.format(time);
    } else {
      Calendar cCal = Calendar.getInstance(dateFormat.getTimeZone());
      date = dateFormat.format(cCal.getTime());
    }

    // turn spaces into underscores
    String test = expType.getName().replace(' ', '_'); // name of experiment
    // make sure parentheses in filenames aren't causing issues
    // (i.e., better than dealing with system-specific setups to escape them)
    test = test.replace('(', '_');
    test = test.replace(')', '_');

    List<String> names = expResult.getInputNames();
    String name = "";
    if (names.size() > getIndexOfMainData()) {
      name = names.get(getIndexOfMainData()); // name of input data
    } else if (names.size() > 0) {
      name = names.get(0);
    }

    return test + '_' + name + '_' + date + ".pdf";

  }

  /**
   * For report generation, give the list of indices of response files used
   * in the plot
   *
   * @return list where each index is a relevant response file
   */
  public int[] getResponseIndices() {
    return expResult.listActiveResponseIndices();
  }

  /**
   * Function to be overridden by implementing class that will add an extra
   * page to PDF reports including charts with less-essential data, such as
   * the plots of residual values over a range for parameter-fitting
   *
   * @return List of charts to show on a second page of PDF reports
   */
  JFreeChart[] getSecondPageCharts() {
    return new JFreeChart[]{};
  }

  /**
   * Default x-axis return function.
   * Though the x-axis is a local variable, some panels may have multiple unit
   * types for the x-axis (i.e., for units of seconds vs. Hz); accessing
   * the x-axis object through this function allows for overrides allowing for
   * more flexibility.
   *
   * @return ValueAxis to be applied to chart
   */
  public ValueAxis getXAxis() {
    return xAxis;
  }

  /**
   * Default y-axis return function. As with getXAxis, designed to be overridden
   * for charts that may use multiple scales.
   *
   * @return ValueAxis to be applied to chart
   */
  public ValueAxis getYAxis() {
    return yAxis;
  }

  /**
   * Function used to query backend on whether or not a DataStore has all the
   * data that a backend needs to calculate. This is used mainly to inform
   * the main window (see SensorSuite class) that the generate result button
   * can be set active
   *
   * @param dataStore DataStore to run data check on
   * @return True if the backend can run with the data provided
   */
  public boolean hasEnoughData(final DataStore dataStore) {
    return expResult.hasEnoughData(dataStore);
  }

  /**
   * True if data has been loaded into the experiment backend yet
   *
   * @return True if there is data to process
   */
  public boolean hasRun() {
    return set;
  }

  /**
   * Get the number of panels to display to fit all data needed by the program
   *
   * @return Number of plots to show in the input panel
   */
  public abstract int panelsNeeded();

  /**
   * Number of panels to return in an output report
   *
   * @return Number of panels to include
   */
  public int plotsToShow() {
    return panelsNeeded();
  }

  /**
   * Takes a PDF and adds a page dedicated to string data related to the
   * data that has been passed into the experiment backend,
   * including any text that might be included in chart title insets and
   * the input data start and end timestamps
   *
   * @param pdf PDF document to append data to
   */
  private void saveInsetDataText(PDDocument pdf) {

    StringBuilder sb = new StringBuilder(expResult.getReportString());
    if (sb.length() > 0) {
      sb.append("\n \n");
    }
    String metadata = getMetadataString();
    if (metadata.length() > 0) {
      sb.append(metadata);
      sb.append("\n \n");
    }
    sb.append(expResult.getFormattedDateRange());
    textToPDFPage(sb.toString(), pdf);
    textListToPDFPages(pdf, getAdditionalReportPages());
  }

  /**
   * Loads in charts used in this panel and prints them out in a PDF document
   * and includes any relevant metadata / analysis results as plain text
   *
   * @param pdf Document to save data to
   */
  public void savePDFResults(PDDocument pdf) {

    int width = 1280;
    int height = 960;
    JFreeChart[] charts = getCharts();

    chartsToPDFPage(width, height, pdf, charts);
    JFreeChart[] page2 = getSecondPageCharts();
    if (page2.length > 0) {
      chartsToPDFPage(width, height, pdf, page2);
    }
    saveInsetDataText(pdf);

  }

  /**
   * Used to plot the results of a backend function from an experiment
   * using a collection of XYSeries mapped by strings. This will be set to
   * the default chart object held by the panel.
   *
   * @param xyDataset collection of XYSeries to plot
   */
  void setChart(XYSeriesCollection xyDataset) {
    chart = buildChart(xyDataset);
  }

  /**
   * Used to identify completion of an experiment to a containing thread
   */
  public void setDone() {
    firePropertyChange("Backend completed", false, set);
    drawCharts();
  }

  /**
   * Used to print out status of the experiment backend onto the chart when
   * the backend status changes
   */
  @Override
  public void stateChanged(ChangeEvent event) {
    if (event.getSource() == expResult) {
      String info = expResult.getStatus();
      displayInfoMessage(info);
    }
  }

  /**
   * Function template for sending input to a backend fucntion and collecting the corresponding data.
   * Details of how to run updateData are left up to the implementing panel however, the boolean "set" should be set to true to enable PDF saving
   *
   * @param dataStore DataStore object containing seed and resp files
   */
  protected abstract void updateData(final DataStore dataStore);


}
