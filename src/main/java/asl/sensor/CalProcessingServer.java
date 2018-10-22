package asl.sensor;

import asl.sensor.experiment.GainExperiment;
import asl.sensor.experiment.GainSixExperiment;
import asl.sensor.output.CalResult;
import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.image.BufferedImage;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.time.Instant;
import java.time.OffsetDateTime;
import java.time.format.DateTimeFormatter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Set;
import javax.imageio.ImageIO;
import org.apache.commons.math3.complex.Complex;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.DateAxis;
import org.jfree.chart.axis.LogarithmicAxis;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.plot.Marker;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.ValueMarker;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.data.xy.XYSeriesCollection;
import asl.sensor.experiment.RandomizedExperiment;
import asl.sensor.experiment.SineExperiment;
import asl.sensor.experiment.StepExperiment;
import asl.sensor.gui.ExperimentPanel;
import asl.sensor.input.DataBlock;
import asl.sensor.input.DataStore;
import asl.sensor.input.InstrumentResponse;
import asl.sensor.utils.ReportingUtils;
import asl.sensor.utils.TimeSeriesUtils;
import edu.iris.dmc.seedcodec.CodecException;
import edu.sc.seis.seisFile.mseed.SeedFormatException;
import py4j.GatewayServer;
import py4j.Py4JNetworkException;

/**
 * CalProcessingServer allows for processing calibrations in a python environment using Py4J.
 *
 * It uses the Py4J default port: 25333
 * If a process is already using that port it silently terminates.
 *
 * @author akearns - KBRWyle
 * @author jholland - USGS
 */
public class CalProcessingServer {


  public CalProcessingServer() {
  }

  /**
   * Enumerate names of embedded resp files
   * @return List of names of embedded resps (will match common sensor & digitizer gain setups)
   */
  public static String[] getEmbeddedRESPFilenames() {
    Set<String> respFilenames = InstrumentResponse.parseInstrumentList();

    List<String> names = new ArrayList<>(respFilenames);
    Collections.sort(names);
    String[] nameArray = new String[names.size()];
    for (int k = 0; k < nameArray.length; ++k) {
      String name = names.get(k);
      name = name.replace("resps/", "");
      nameArray[k] = name;
    }
    return nameArray;
  }

  public static void main(String[] args) {
    GatewayServer gatewayServer = new GatewayServer(new CalProcessingServer());
    try {
      gatewayServer.start();
    } catch (Py4JNetworkException e) {
      System.exit(0);
    }
    System.out.println("Gateway Server Started");
  }

  /**
   * Acquire data and run a gain experiment over it. Angle and gain references can be set
   * independently using the command line parameters. Each RESP file can be set as embedded
   * individually.
   *
   * @param north1FileName File name of data to be used as first north input (N1)
   * @param east1FileName File name of data to be used as first east input (E1)
   * @param vert1FileName File name of data to be used as first vertical input (Z1)
   * @param north2FileName File name of data to be used as second north input (N2)
   * @param east2FileName File name of data to be used as second east input (E2)
   * @param vert2FileName File name of data to be used as second vertical input (Z2)
   * @param north1RespName Name of data to be used as first north response (N1)
   * @param east1RespName Name of data to be used as first east response (E1)
   * @param vert1RespName Name of data to be used as first vertical response (Z1)
   * @param north2RespName Name of data to be used as second north response (N2)
   * @param east2RespName Name of data to be used as second east response (E2)
   * @param vert2RespName Name of data to be used as second vertical response (Z2)
   * @param north1RespEmbedded True if N1 resp is an embedded resp file
   * @param east1RespEmbedded True if E1 resp is an embedded resp file
   * @param vert1RespEmbedded True if Z1 resp is an embedded resp file
   * @param north2RespEmbedded True if N2 resp is an embedded resp file
   * @param east2RespEmbedded True if E2 resp is an embedded resp file
   * @param vert2RespEmbedded True if Z2 resp is an embedded resp file
   * @param startDate ISO-861 formatted datetime string with timezone offset; start of data window
   * @param endDate ISO-861 formatted datetime string with timezone offset; end of data window
   * @param useFirstDataAsAngleRef True if N-E1 data will be used for rotation reference
   * @param useFirstDataAsGainRef True if N-E-Z1 data will be used for rotation reference
   * @return Data from running the experiment (plots and gain statistics)
   * @throws IOException If a string does not refer to a valid accessible file
   * @throws SeedFormatException If a data file cannot be parsed as a seed file
   * @throws CodecException If there is an issue with the compression of the seed files
   */
  public CalResult runGain(String north1FileName, String east1FileName, String vert1FileName,
      String north2FileName, String east2FileName, String vert2FileName, String north1RespName,
      String east1RespName, String vert1RespName, String north2RespName, String east2RespName,
      String vert2RespName, boolean north1RespEmbedded, boolean east1RespEmbedded,
      boolean vert1RespEmbedded, boolean north2RespEmbedded, boolean east2RespEmbedded,
      boolean vert2RespEmbedded, String startDate, String endDate, boolean useFirstDataAsAngleRef,
      boolean useFirstDataAsGainRef) throws IOException, CodecException, SeedFormatException {
    DateTimeFormatter dtf = DateTimeFormatter.ISO_OFFSET_DATE_TIME;
    OffsetDateTime startDateTime = OffsetDateTime.parse(startDate, dtf);
    OffsetDateTime endDateTime = OffsetDateTime.parse(endDate, dtf);
    long start = startDateTime.toInstant().toEpochMilli();
    long end = endDateTime.toInstant().toEpochMilli();

    String[] seedFileNames = new String[]{
        north1FileName, east1FileName, vert1FileName, north2FileName, east2FileName, vert2FileName};
    String[] respFileNames = new String[]{
        north1RespName, east1RespName, vert1RespName, north2RespName, east2RespName, vert2RespName};
    boolean[] embedResps = new boolean[]{
        north1RespEmbedded, east1RespEmbedded, vert1RespEmbedded,
        north2RespEmbedded, east2RespEmbedded, vert2RespEmbedded};

    DataStore ds = new DataStore();
    for (int i = 0; i < seedFileNames.length; ++i) {
      DataBlock db = TimeSeriesUtils.getFirstTimeSeries(seedFileNames[i]);
      ds.setBlock(i, db);
      InstrumentResponse ir;
      if(embedResps[i]) {
        ir = InstrumentResponse.loadEmbeddedResponse(respFileNames[i]);
      } else {
        ir = new InstrumentResponse(respFileNames[i]);
      }
      ds.setResponse(i, ir);
    }
    ds.trim(start, end);

    return runExpGetDataGain(ds, useFirstDataAsAngleRef, useFirstDataAsGainRef);
  }

  private CalResult runExpGetDataGain(DataStore ds, boolean firstAngleRef, boolean firstGainRef)
      throws IOException {

    // input is sorted such that the gain reference is always plotted first;
    // so if this is false, this means second set of data in plot set is reference angle
    boolean refDataMatches = (firstAngleRef == firstGainRef);

    GainSixExperiment gainSix = new GainSixExperiment();

    if (firstAngleRef) {
      gainSix.setFirstDataAsAngleReference();
    } else {
      gainSix.setSecondDataAsAngleReference();
    }

    int gainRef = firstGainRef ? 0 : 1;
    gainSix.setReferenceIndex(gainRef);

    gainSix.setRangeForStatistics(GainExperiment.DEFAULT_LOW_BOUND,
        GainExperiment.DEFAULT_UP_BOUND);

    gainSix.runExperimentOnData(ds);

    List<XYSeriesCollection> results = gainSix.getData();

    // plot has 3 components: source, destination, NLNM line plot
    String[] orientation = new String[]{"North", "East", "Vertical"};
    JFreeChart[] charts = new JFreeChart[3];
    for (int i = 0; i < charts.length; ++i) {
      XYSeriesCollection timeseriesIn = results.get(i);
      XYSeriesCollection timeseries = new XYSeriesCollection();
      timeseries.addSeries(timeseriesIn.getSeries(gainRef));
      timeseries.addSeries(timeseriesIn.getSeries((gainRef + 1) % 2));
      timeseries.addSeries(timeseriesIn.getSeries("NLNM"));

      charts[i] = ChartFactory.createXYLineChart(
          "Gain Experiment -- " + orientation[i],
          "Period (s)",
          "Power (rel. 1 (m/s^2)^2/Hz)",
          timeseries,
          PlotOrientation.VERTICAL,
          true, // include legend
          false,
          false);

      // add vertical lines to plot over rage of data for statistics (3 to 9s by default)
      XYPlot plot = charts[i].getXYPlot();
      Marker startMarker = new ValueMarker(GainExperiment.DEFAULT_LOW_BOUND);
      startMarker.setStroke(new BasicStroke((float) 1.5));
      Marker endMarker = new ValueMarker(GainExperiment.DEFAULT_UP_BOUND);
      endMarker.setStroke(new BasicStroke((float) 1.5));
      plot.addDomainMarker(startMarker);
      plot.addDomainMarker(endMarker);

      // ensure that NLNM lines are bolder, colored black
      XYItemRenderer renderer = plot.getRenderer();
      // series index 0 - ref data; series index 1 - other data; series index 2 - NLNM plot
      BasicStroke stroke = (BasicStroke) renderer.getSeriesStroke(2);
      if (stroke == null) {
        stroke = (BasicStroke) renderer.getBaseStroke();
      }
      stroke = new BasicStroke(stroke.getLineWidth() * 2);
      renderer.setSeriesStroke(2, stroke);
      renderer.setSeriesPaint(2, new Color(0, 0, 0));
    }

    double northAzimuth = gainSix.getNorthAzimuthDegrees();
    double eastAzimuth = gainSix.getEastAzimuthDegrees();
    double[][] statistics = gainSix.getStatistics();

    BufferedImage[] images = ReportingUtils.chartsToImageList(1, 1280, 960, charts);
    byte[][] pngByteArrays = new byte[images.length][];
    for (int i = 0; i < images.length; ++i) {
      ByteArrayOutputStream out = new ByteArrayOutputStream();
      ImageIO.write(images[i], "png", out);
      pngByteArrays[i] = out.toByteArray();
    }

    return CalResult.buildSixGainData(refDataMatches, northAzimuth, eastAzimuth,
        statistics, pngByteArrays);

  }

  /**
   * Acquire data and run randomized calibration solver over it.
   * Returns the experiment (all data kept locally to maintain thread safety)
   *
   * @param calFileName Filename of calibration signal
   * @param outFileName Filename of sensor output
   * @param respName Filename of response to load in
   * @param useEmbeddedResp True if response is an embedded response in program
   * @param startDate ISO-861 formatted datetime string with timezone offset; start of data window
   * @param endDate ISO-861 formatted datetime string with timezone offset; end of data window
   * @param lowFreq True if a low-freq cal should be run
   * @return Data from running the experiment (plots and fit pole/zero values)
   * @throws IOException If a string does not refer to a valid accessible file
   * @throws SeedFormatException If a data file cannot be parsed as a seed file
   * @throws CodecException If there is an issue with the compression of the seed files
   */
  public CalResult runRand(String calFileName, String outFileName,
      String respName, boolean useEmbeddedResp, String startDate, String endDate, boolean lowFreq)
      throws IOException, SeedFormatException, CodecException {
    DateTimeFormatter dtf = DateTimeFormatter.ISO_OFFSET_DATE_TIME;
    OffsetDateTime startDateTime = OffsetDateTime.parse(startDate, dtf);
    OffsetDateTime endDateTime = OffsetDateTime.parse(endDate, dtf);
    long start = startDateTime.toInstant().toEpochMilli();
    long end = endDateTime.toInstant().toEpochMilli();

    DataStore ds = new DataStore();
    DataBlock calBlock = TimeSeriesUtils.getFirstTimeSeries(calFileName);
    DataBlock outBlock = TimeSeriesUtils.getFirstTimeSeries(outFileName);
    InstrumentResponse ir;
    if (useEmbeddedResp) {
      ir = InstrumentResponse.loadEmbeddedResponse(respName);
    } else {
      Instant epoch = InstrumentResponse.getRespFileClosestEpoch(respName, start, end);
      ir = new InstrumentResponse(respName, epoch);
    }

    ds.setBlock(0, calBlock);
    ds.setBlock(1, outBlock);
    ds.setResponse(1, ir);
    ds.trim(start, end);
    if (lowFreq) {
      ds.resample(10.); // more than 5 Hz should be unnecessary for low-frequency curve fitting
    }

    return runExpGetDataRand(ds, lowFreq);
  }

  /**
   * Acquire data and run randomized calibration solver over it.
   * Used to handle calibrations that cross day boundaries.
   * Returns the experiment (all data kept locally to maintain thread safety)
   *
   * @param calFileNameD1 Filename of calibration signal (day 1)
   * @param calFileNameD2 Filename of calibration signal (day 2)
   * @param outFileNameD1 Filename of sensor output (day 1)
   * @param outFileNameD2 Filename of sensor output (day 2)
   * @param respName Filename of response to load in
   * @param useEmbeddedResp True if response is an embedded response in program
   * @param startDate ISO-861 formatted datetime string with timezone offset; start of data window
   * @param endDate ISO-861 formatted datetime string with timezone offset; end of data window
   * @param lowFreq True if a low-freq cal should be run
   * @return Data from running the experiment (plots and fit pole/zero values)
   * @throws IOException If a string does not refer to a valid accessible file
   * @throws SeedFormatException If a data file cannot be parsed as a seed file
   * @throws CodecException If there is an issue with the compression of the seed files
   */
  public CalResult runRand(String calFileNameD1, String calFileNameD2,
      String outFileNameD1, String outFileNameD2, String respName, boolean useEmbeddedResp,
      String startDate, String endDate, boolean lowFreq)
      throws IOException, SeedFormatException, CodecException {
    DateTimeFormatter dtf = DateTimeFormatter.ISO_OFFSET_DATE_TIME;
    OffsetDateTime startDateTime = OffsetDateTime.parse(startDate, dtf);
    OffsetDateTime endDateTime = OffsetDateTime.parse(endDate, dtf);

    long start = startDateTime.toInstant().toEpochMilli();
    long end = endDateTime.toInstant().toEpochMilli();

    DataStore ds = new DataStore();
    String[] calFileName = new String[]{calFileNameD1, calFileNameD2};
    String[] outFileName = new String[]{outFileNameD1, outFileNameD2};
    DataBlock calBlock = TimeSeriesUtils.getFirstTimeSeries(calFileName);
    DataBlock outBlock = TimeSeriesUtils.getFirstTimeSeries(outFileName);
    InstrumentResponse ir;
    if (useEmbeddedResp) {
      ir = InstrumentResponse.loadEmbeddedResponse(respName);
    } else {
      Instant epoch = InstrumentResponse.getRespFileClosestEpoch(respName, start, end);
      ir = new InstrumentResponse(respName, epoch);
    }

    ds.setBlock(0, calBlock);
    ds.setBlock(1, outBlock);
    ds.setResponse(1, ir);
    ds.trim(start, end);

    return runExpGetDataRand(ds, lowFreq);
  }

  /**
   * Acquire data and run step calibration solver over it.
   * Used to handle calibrations that cross day boundaries.
   * Returns the experiment (all data kept locally to maintain thread safety)
   *
   * @param calFileNameD1 Filename of calibration signal (day 1)
   * @param calFileNameD2 Filename of calibration signal (day 2)
   * @param outFileNameD1 Filename of sensor output (day 1)
   * @param outFileNameD2 Filename of sensor output (day 2)
   * @param respName Filename of response to load in
   * @param useEmbeddedResp True if response is an embedded response in program
   * @param startDate ISO-861 formatted datetime string with timezone offset; start of data window
   * @param endDate ISO-861 formatted datetime string with timezone offset; end of data window
   * @return Data from running the experiment (plots and fit corner/damping values)
   * @throws IOException If a string does not refer to a valid accessible file
   * @throws SeedFormatException If a data file cannot be parsed as a seed file
   * @throws CodecException If there is an issue with the compression of the seed files
   */
  public CalResult runStep(String calFileNameD1, String calFileNameD2, String outFileNameD1,
      String outFileNameD2, String respName, boolean useEmbeddedResp, String startDate,
      String endDate) throws SeedFormatException, CodecException, IOException {
    DateTimeFormatter dtf = DateTimeFormatter.ISO_OFFSET_DATE_TIME;
    OffsetDateTime startDateTime = OffsetDateTime.parse(startDate, dtf);
    OffsetDateTime endDateTime = OffsetDateTime.parse(endDate, dtf);

    long start = startDateTime.toInstant().toEpochMilli();
    long end = endDateTime.toInstant().toEpochMilli();

    DataStore ds = new DataStore();
    String[] calFileName = new String[]{calFileNameD1, calFileNameD2};
    String[] outFileName = new String[]{outFileNameD1, outFileNameD2};
    DataBlock calBlock = TimeSeriesUtils.getFirstTimeSeries(calFileName);
    DataBlock outBlock = TimeSeriesUtils.getFirstTimeSeries(outFileName);
    InstrumentResponse ir;
    if (useEmbeddedResp) {
      ir = InstrumentResponse.loadEmbeddedResponse(respName);
    } else {
      Instant epoch = InstrumentResponse.getRespFileClosestEpoch(respName, start, end);
      ir = new InstrumentResponse(respName, epoch);
    }

    ds.setBlock(0, calBlock);
    ds.setBlock(1, outBlock);
    ds.setResponse(1, ir);
    ds.trim(start, end);

    return runExpGetDataStep(ds);

  }

  /**
   * Acquire data and run step calibration solver over it.
   * Returns the experiment (all data kept locally to maintain thread safety)
   *
   * @param calFileName Filename of calibration signal
   * @param outFileName Filename of sensor output
   * @param respName Filename of response to load in
   * @param useEmbeddedResp True if response is an embedded response in program
   * @param startDate ISO-861 formatted datetime string with timezone offset; start of data window
   * @param endDate ISO-861 formatted datetime string with timezone offset; end of data window
   * @return Data from running the experiment (plots and fit corner/damping values)
   * @throws IOException If a string does not refer to a valid accessible file
   * @throws SeedFormatException If a data file cannot be parsed as a seed file
   * @throws CodecException If there is an issue with the compression of the seed files
   */
  public CalResult runStep(String calFileName, String outFileName, String respName,
      boolean useEmbeddedResp, String startDate,String endDate)
      throws SeedFormatException, CodecException, IOException {
    DateTimeFormatter dtf = DateTimeFormatter.ISO_OFFSET_DATE_TIME;
    OffsetDateTime startDateTime = OffsetDateTime.parse(startDate, dtf);
    OffsetDateTime endDateTime = OffsetDateTime.parse(endDate, dtf);

    long start = startDateTime.toInstant().toEpochMilli();
    long end = endDateTime.toInstant().toEpochMilli();

    DataStore ds = new DataStore();
    DataBlock calBlock = TimeSeriesUtils.getFirstTimeSeries(calFileName);
    DataBlock outBlock = TimeSeriesUtils.getFirstTimeSeries(outFileName);
    InstrumentResponse ir;
    if (useEmbeddedResp) {
      ir = InstrumentResponse.loadEmbeddedResponse(respName);
    } else {
      Instant epoch = InstrumentResponse.getRespFileClosestEpoch(respName, start, end);
      ir = new InstrumentResponse(respName, epoch);
    }

    ds.setBlock(0, calBlock);
    ds.setBlock(1, outBlock);
    ds.setResponse(1, ir);
    ds.trim(start, end);

    return runExpGetDataStep(ds);

  }

  /**
   * Acquire data and run sine calibration solver over it.
   * Returns the experiment (all data kept locally to maintain thread safety)
   *
   * @param calFileName Filename of calibration signal
   * @param outFileName Filename of sensor output
   * @param startDate ISO-861 formatted datetime string with timezone offset; start of data window
   * @param endDate ISO-861 formatted datetime string with timezone offset; end of data window
   * @return Data from running the experiment (plots and amplitude estimations)
   * @throws IOException If a string does not refer to a valid accessible file
   * @throws SeedFormatException If a data file cannot be parsed as a seed file
   * @throws CodecException If there is an issue with the compression of the seed files
   */
  public CalResult runSine(String calFileName, String outFileName, String startDate,
      String endDate) throws SeedFormatException, CodecException, IOException {
    DateTimeFormatter dtf = DateTimeFormatter.ISO_OFFSET_DATE_TIME;
    OffsetDateTime startDateTime = OffsetDateTime.parse(startDate, dtf);
    OffsetDateTime endDateTime = OffsetDateTime.parse(endDate, dtf);

    long start = startDateTime.toInstant().toEpochMilli();
    long end = endDateTime.toInstant().toEpochMilli();

    DataStore ds = new DataStore();
    DataBlock calBlock = TimeSeriesUtils.getFirstTimeSeries(calFileName);
    DataBlock outBlock = TimeSeriesUtils.getFirstTimeSeries(outFileName);

    ds.setBlock(0, calBlock);
    ds.setBlock(1, outBlock);
    ds.trim(start, end);

    return runExpGetDataSine(ds);
  }

  /**
   * Acquire data and run sine calibration solver over it.
   * Used to handle calibrations that cross day boundaries.
   * Returns the experiment (all data kept locally to maintain thread safety)
   *
   * @param calFileNameD1 Filename of calibration signal (day 1)
   * @param calFileNameD2 Filename of calibration signal (day 2)
   * @param outFileNameD1 Filename of sensor output (day 1)
   * @param outFileNameD2 Filename of sensor output (day 2)
   * @param startDate ISO-861 formatted datetime string with timezone offset; start of data window
   * @param endDate ISO-861 formatted datetime string with timezone offset; end of data window
   * @return Data from running the experiment (plots and amplitude estimations)
   * @throws IOException If a string does not refer to a valid accessible file
   * @throws SeedFormatException If a data file cannot be parsed as a seed file
   * @throws CodecException If there is an issue with the compression of the seed files
   */
  public CalResult runSine(String calFileNameD1, String calFileNameD2, String outFileNameD1,
      String outFileNameD2, String startDate, String endDate)
      throws SeedFormatException, CodecException, IOException {
    DateTimeFormatter dtf = DateTimeFormatter.ISO_OFFSET_DATE_TIME;
    OffsetDateTime startDateTime = OffsetDateTime.parse(startDate, dtf);
    OffsetDateTime endDateTime = OffsetDateTime.parse(endDate, dtf);

    long start = startDateTime.toInstant().toEpochMilli();
    long end = endDateTime.toInstant().toEpochMilli();

    DataStore ds = new DataStore();
    String[] calFileName = new String[]{calFileNameD1, calFileNameD2};
    String[] outFileName = new String[]{outFileNameD1, outFileNameD2};
    DataBlock calBlock = TimeSeriesUtils.getFirstTimeSeries(calFileName);
    DataBlock outBlock = TimeSeriesUtils.getFirstTimeSeries(outFileName);

    ds.setBlock(0, calBlock);
    ds.setBlock(1, outBlock);
    ds.trim(start, end);

    return runExpGetDataSine(ds);
  }

  private CalResult runExpGetDataSine(DataStore ds) throws IOException {
    SineExperiment sine = new SineExperiment();
    sine.runExperimentOnData(ds);
    List<XYSeriesCollection> plots = sine.getData();
    double calAmplitude = sine.getCalAmplitude();
    double outAmplitude = sine.getOutAmplitude();
    double estFreq = sine.getEstSineFreq();
    double ratio = calAmplitude / outAmplitude;

    DateAxis timeAxis = new DateAxis();
    timeAxis.setDateFormatOverride(ExperimentPanel.DATE_TIME_FORMAT.get());
    Font bold = timeAxis.getLabelFont().deriveFont(Font.BOLD);
    timeAxis.setLabelFont(bold);

    JFreeChart sineChart = ChartFactory.createXYLineChart(
        "Sine Calibration",
        "Time of sample (Julian date)",
        "Normalized calibration signals (counts)",
        plots.get(0));
    sineChart.getXYPlot().setDomainAxis(timeAxis);

    JFreeChart linearityChart = ChartFactory.createXYLineChart(
        "Sine cal. Linearity",
        "Value of sampled cal data (counts)",
        "Value of sampled sensor output (counts)",
        plots.get(1));

    JFreeChart[] charts = {sineChart, linearityChart};

    BufferedImage[] images = ReportingUtils.chartsToImageList(1, 1280, 960, charts);
    byte[][] pngByteArrays = new byte[images.length][];
    for (int i = 0; i < images.length; ++i) {
      ByteArrayOutputStream out = new ByteArrayOutputStream();
      ImageIO.write(images[i], "png", out);
      pngByteArrays[i] = out.toByteArray();
    }
    return CalResult.buildSineCalData(pngByteArrays, calAmplitude, outAmplitude, estFreq, ratio);
  }

  private CalResult runExpGetDataStep(DataStore ds) throws IOException {
    StepExperiment step = new StepExperiment();
    step.runExperimentOnData(ds);
    double[] fitParams = step.getFitParams();
    double[] initParams = step.getInitParams();
    List<XYSeriesCollection> plots = step.getData();
    // order of plots -- step function, resp amplitudes, resp phases

    NumberAxis stepAxis = new NumberAxis("Step counts");
    DateAxis timeAxis = new DateAxis("Time of sample (Julian date)");
    timeAxis.setDateFormatOverride(ExperimentPanel.DATE_TIME_FORMAT.get());
    NumberAxis ampAxis = new NumberAxis("RESP Amplitude [10 * log10(RESP(f))]");
    NumberAxis phaseAxis = new NumberAxis("RESP Phase (deg.)");
    LogarithmicAxis freqAxis = new LogarithmicAxis("Frequency (f)");
    Font bold = stepAxis.getLabelFont().deriveFont(Font.BOLD);
    stepAxis.setLabelFont(bold);
    timeAxis.setLabelFont(bold);
    ampAxis.setLabelFont(bold);
    phaseAxis.setLabelFont(bold);
    freqAxis.setLabelFont(bold);

    JFreeChart stepChart = ChartFactory.createXYLineChart(
        "Step Calibration",
        timeAxis.getLabel(),
        stepAxis.getLabel(),
        plots.get(0));
    stepChart.getXYPlot().setDomainAxis(timeAxis);
    stepChart.getXYPlot().setRangeAxis(stepAxis);

    JFreeChart respAmpChart = ChartFactory.createXYLineChart(
        "Step Calibration - Resp Amplitude Comparison",
        freqAxis.getLabel(),
        ampAxis.getLabel(),
        plots.get(1));
    respAmpChart.getXYPlot().setDomainAxis(freqAxis);
    respAmpChart.getXYPlot().setRangeAxis(ampAxis);

    JFreeChart respPhaseChart = ChartFactory.createXYLineChart(
        "Step Calibration - Resp Phase Comparison",
        freqAxis.getLabel(),
        phaseAxis.getLabel(),
        plots.get(2));
    respPhaseChart.getXYPlot().setDomainAxis(freqAxis);
    respPhaseChart.getXYPlot().setRangeAxis(phaseAxis);
    JFreeChart[] charts = {stepChart, respAmpChart, respPhaseChart};

    BufferedImage[] images = ReportingUtils.chartsToImageList(1, 1280, 960, charts);
    byte[][] pngByteArrays = new byte[images.length][];
    for (int i = 0; i < images.length; ++i) {
      ByteArrayOutputStream out = new ByteArrayOutputStream();
      ImageIO.write(images[i], "png", out);
      pngByteArrays[i] = out.toByteArray();
    }

    return CalResult.buildStepCalData(pngByteArrays, fitParams, initParams);
  }

  private CalResult runExpGetDataRand(DataStore dataStore, boolean isLowFrequency)
      throws IOException {

    RandomizedExperiment randomExperiment = new RandomizedExperiment();

    randomExperiment.setLowFrequencyCalibration(isLowFrequency);
    randomExperiment.runExperimentOnData(dataStore);

    Complex[] fitZerosComplex = randomExperiment.getFitResponse().getZeros()
        .toArray(new Complex[]{});
    Complex[] fitPolesComplex = randomExperiment.getFitResponse().getPoles()
        .toArray(new Complex[]{});
    Complex[] initialZerosComplex = dataStore.getResponse(1).getZeros().toArray(new Complex[]{});
    Complex[] initialPolesComplex = dataStore.getResponse(1).getPoles().toArray(new Complex[]{});

    double[] fitZeros = new double[2 * fitZerosComplex.length];
    double[] initialZeros = new double[fitZeros.length];
    double[] fitPoles = new double[2 * fitPolesComplex.length];
    double[] initialPoles = new double[fitPoles.length];
    for (int i = 0; i < fitZerosComplex.length; ++i) {
      int realIndex = 2 * i;
      int imaginaryIndex = realIndex + 1;
      fitZeros[realIndex] = fitZerosComplex[i].getReal();
      fitZeros[imaginaryIndex] = fitZerosComplex[i].getImaginary();
      initialZeros[realIndex] = initialZerosComplex[i].getReal();
      initialZeros[imaginaryIndex] = initialZerosComplex[i].getImaginary();
    }
    for (int i = 0; i < fitPolesComplex.length; ++i) {
      int realIndex = 2 * i;
      int imaginaryIndex = realIndex + 1;
      fitPoles[realIndex] = fitPolesComplex[i].getReal();
      fitPoles[imaginaryIndex] = fitPolesComplex[i].getImaginary();
      initialPoles[realIndex] = initialPolesComplex[i].getReal();
      initialPoles[imaginaryIndex] = initialPolesComplex[i].getImaginary();
    }

    List<XYSeriesCollection> xySeriesCollections = randomExperiment.getData();
    JFreeChart[] charts = new JFreeChart[xySeriesCollections.size()];

    String xAxisTitle = "Frequency (Hz)";
    String amplitudeAxisTitle = "20 * log10( RESP(f) )";
    String phaseAxisTitle = "phi(RESP(f))";

    ValueAxis xAxis = new LogarithmicAxis(xAxisTitle);
    ValueAxis residualXAxis = new LogarithmicAxis(xAxisTitle);
    NumberAxis amplitudeAxis = new NumberAxis(amplitudeAxisTitle);
    amplitudeAxis.setAutoRange(true);
    amplitudeAxis.setAutoRangeIncludesZero(false);
    ValueAxis phaseAxis = new NumberAxis(phaseAxisTitle);
    phaseAxis.setAutoRange(true);
    ValueAxis residualPhaseAxis = new NumberAxis("Phase error (degrees)");
    ValueAxis residualAmplitudeAxis = new NumberAxis("Amplitude error (percentage)");
    Font bold = xAxis.getLabelFont().deriveFont(Font.BOLD);
    xAxis.setLabelFont(bold);
    amplitudeAxis.setLabelFont(bold);
    phaseAxis.setLabelFont(bold);
    residualXAxis.setLabelFont(bold);
    residualPhaseAxis.setLabelFont(bold);
    residualAmplitudeAxis.setLabelFont(bold);
    XYPlot xyPlot;

    String title;
    if (isLowFrequency) {
      title = "Low-freq random cal";
    } else {
      title = "High-freq random cal";
    }

    charts[0] = ChartFactory.createXYLineChart(
        title + " Amplitude",
        xAxis.getLabel(),
        amplitudeAxis.getLabel(),
        xySeriesCollections.get(0),
        PlotOrientation.VERTICAL,
        true, // include legend
        false,
        false);
    xyPlot = charts[0].getXYPlot();
    xyPlot.setDomainAxis(xAxis);
    xyPlot.setRangeAxis(amplitudeAxis);
    ExperimentPanel.invertSeriesRenderingOrder(charts[0]);

    charts[1] = ChartFactory.createXYLineChart(
        title + " Phase",
        xAxis.getLabel(),
        phaseAxis.getLabel(),
        xySeriesCollections.get(1),
        PlotOrientation.VERTICAL,
        true, // include legend
        false,
        false);
    xyPlot = charts[1].getXYPlot();
    xyPlot.setDomainAxis(xAxis);
    xyPlot.setRangeAxis(phaseAxis);
    ExperimentPanel.invertSeriesRenderingOrder(charts[1]);

    charts[2] = ChartFactory.createXYLineChart(
        title + " Amplitude Error",
        residualXAxis.getLabel(),
        residualAmplitudeAxis.getLabel(),
        xySeriesCollections.get(2),
        PlotOrientation.VERTICAL,
        true, // include legend
        false,
        false);
    xyPlot = charts[2].getXYPlot();
    xyPlot.setDomainAxis(residualXAxis);
    xyPlot.setRangeAxis(residualAmplitudeAxis);
    xyPlot.getRenderer().setSeriesPaint(1, Color.GREEN);
    ExperimentPanel.invertSeriesRenderingOrder(charts[2]);

    charts[3] = ChartFactory.createXYLineChart(
        title + "Phase Error",
        residualXAxis.getLabel(),
        residualPhaseAxis.getLabel(),
        xySeriesCollections.get(3),
        PlotOrientation.VERTICAL,
        true, // include legend
        false,
        false);
    xyPlot = charts[3].getXYPlot();
    xyPlot.setDomainAxis(residualXAxis);
    xyPlot.setRangeAxis(residualPhaseAxis);
    xyPlot.getRenderer().setSeriesPaint(1, Color.GREEN);
    ExperimentPanel.invertSeriesRenderingOrder(charts[3]);

    if (!isLowFrequency) {
      Marker maxFitMarker = new ValueMarker(randomExperiment.getMaxFitFrequency());
      maxFitMarker.setStroke(new BasicStroke((float) 1.5));
      charts[0].getXYPlot().addDomainMarker(maxFitMarker);
      charts[1].getXYPlot().addDomainMarker(maxFitMarker);
    }

    BufferedImage[] images = ReportingUtils.chartsToImageList(1, 1280, 960, charts);
    byte[][] pngByteArrays = new byte[images.length][];
    for (int i = 0; i < images.length; ++i) {
      ByteArrayOutputStream out = new ByteArrayOutputStream();
      ImageIO.write(images[i], "png", out);
      pngByteArrays[i] = out.toByteArray();
    }

    return CalResult.buildRandomCalData(fitPoles, fitZeros, initialPoles, initialZeros,
        pngByteArrays);

  }

}
