package asl.sensor;

import java.awt.BasicStroke;
import java.awt.Font;
import java.awt.image.BufferedImage;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.time.OffsetDateTime;
import java.time.format.DateTimeFormatter;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TimeZone;
import javax.imageio.ImageIO;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.util.Pair;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.LogarithmicAxis;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.plot.Marker;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.ValueMarker;
import org.jfree.chart.plot.XYPlot;
import org.jfree.data.xy.XYSeriesCollection;
import asl.sensor.experiment.RandomizedExperiment;
import asl.sensor.experiment.SineExperiment;
import asl.sensor.experiment.StepExperiment;
import asl.sensor.gui.ExperimentPanel;
import asl.sensor.gui.RandomizedPanel;
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


  @SuppressWarnings({"WeakerAccess", "unused"})
  public abstract class CalResult {
    protected Map<String, double[]> numerMap;
    protected Map<String, byte[]> imageMap;

    public CalResult() {
      numerMap = new HashMap<String, double[]>();
      imageMap = new HashMap<String, byte[]>();
    }

    public Map<String, byte[]> getImageMap() {
      return imageMap;
    }

    public Map<String, double[]> getNumerMap() {
      return numerMap;
    }
  }

  public class StepData extends CalResult {
    // constuctor to be used with step calibrations
    StepData(byte[][] images, double[] initParams, double[] fitParams) {
      super();
      double fitCorner = fitParams[0];
      double fitDamping = fitParams[1];
      double fitResid = fitParams[2];
      double initCorner = initParams[0];
      double initDamping = initParams[1];
      double initResid = initParams[2];
      numerMap.put("Fit-corner", new double[]{fitCorner});
      numerMap.put("Fit-damping", new double[]{fitDamping});
      numerMap.put("Fit-residual", new double[]{fitResid});
      numerMap.put("Initial-corner", new double[]{initCorner});
      numerMap.put("Initial-damping", new double[]{initDamping});
      numerMap.put("Initial-residual", new double[]{initResid});
      imageMap.put("Step-plot", images[0]);
      imageMap.put("Resp-amplitudes", images[1]);
      imageMap.put("Resp-phases", images[2]);
    }
  }

  public class SineData extends CalResult {
    SineData(byte[][] images, double calAmp, double outAmp, double freq, double ratio) {
      super();
      numerMap.put("Calibration-amplitude", new double[]{calAmp});
      numerMap.put("Output-signal-amplitude", new double[]{outAmp});
      numerMap.put("Estimated signal frequency", new double[]{freq});
      numerMap.put("Calibration-to-output ratio", new double[]{ratio});
      imageMap.put("Sines-plot", images[0]);
      imageMap.put("Linearity", images[1]);
    }
  }

  //RandData is used externally and requires Public access on getters
  @SuppressWarnings({"WeakerAccess", "unused"})
  public class RandData {

    private double[] initialPoles;
    private double[] initialZeros;
    private double[] fitPoles;
    private double[] fitZeros;
    private byte[][] images;
    private String[] gapNameIdentifiers;
    private Date[][] gapStarts;
    private Date[][] gapEnds;

    RandData(double[] fitPoles, double[] fitZeros, double[] initialPoles, double[] initialZeros, byte[][] images,
        String[] gapNames, Date[][] gapStartTimes, Date[][] gapEndTimes) {
      this.fitPoles = fitPoles;
      this.fitZeros = fitZeros;
      this.initialPoles = initialPoles;
      this.initialZeros = initialZeros;
      this.images = images;
      gapNameIdentifiers = gapNames;
      gapStarts = gapStartTimes;
      gapEnds = gapEndTimes;
    }

    public Map<String, byte[]> getImageMap() {
      Map<String, byte[]> imageMap = new HashMap<String, byte[]>();
      imageMap.put("amplitude", images[0]);
      imageMap.put("phase", images[1]);
      imageMap.put("amplitude-error", images[2]);
      imageMap.put("phase-error", images[3]);
      return imageMap;
    }

    public Map<String, double[]> getNumericMap() {
      Map<String, double[]> numerMap = new HashMap<String, double[]>();
      numerMap.put("fit-poles", fitPoles);
      numerMap.put("fit-zeros", fitZeros);
      numerMap.put("initial-poles", initialPoles);
      numerMap.put("initial-zeros", initialZeros);
      return numerMap;
    }

    public byte[] getAmpErrorImage() {
      return images[2];
    }

    public byte[] getAmpImage() {
      return images[0];
    }

    public double[] getFitPoles() {
      return fitPoles;
    }

    public double[] getFitZeros() {
      return fitZeros;
    }

    public Date[][] getGapEndDates() {
      return gapEnds;
    }

    public String[] getGapIdentifiers() {
      return gapNameIdentifiers;
    }

    public String getGapInfoAsString() {
      SimpleDateFormat sdf = new SimpleDateFormat("DD.HH:m:s");
      sdf.setTimeZone(TimeZone.getTimeZone("UTC"));
      return getGapInfoAsString(sdf);
    }

    public String getGapInfoAsString(DateFormat df) {
      StringBuilder sb = new StringBuilder();
      for (int j = 0; j < gapNameIdentifiers.length; ++j) {
        sb.append(gapNameIdentifiers[j]);
        sb.append(":\n");
        for (int i = 0; i < gapStarts[j].length; ++i) {
          sb.append("\t");
          Date start = gapStarts[j][i];
          Date end = gapEnds[j][i];
          sb.append(df.format(start));
          sb.append("\t");
          sb.append(df.format(end));
          sb.append("\n");
        }
        sb.append("\n");
      }
      return sb.toString();
    }

    public Date[][] getGapStartDates() {
      return gapStarts;
    }

    public double[] getInitPoles() {
      return initialPoles;
    }

    public double[] getInitZeros() {
      return initialZeros;
    }

    public byte[] getPhaseErrorImage() {
      return images[3];
    }

    public byte[] getPhaseImage() {
      return images[1];
    }
  }

  /**
   * Get all metadata from the function in a single file
   *
   * @return text representation of data from experiment
   */
  @SuppressWarnings("unused")
  public static String getMetadataFromExp(RandomizedExperiment exp) {
    String[] data = RandomizedPanel.getInsetString(exp);
    StringBuilder sb = new StringBuilder();
    for (String dataPart : data) {
      sb.append(dataPart);
      sb.append('\n');
    }
    return sb.toString();
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

  public CalProcessingServer() {
  }

  /**
   * Acquire data and run calibration over it.
   * Returns the experiment (all data kept locally to maintain thread safety)
   *
   * @param calFileName Filename of calibration signal
   * @param outFileName Filename of sensor output
   * @param respName Filename of response to load in
   * @param useEmbeddedResp True if response is an embedded response in program
   * @param startDate Long representing ms-since-epoch of data start time
   * @param endDate Long representing ms-since-epoch of data end time
   * @param lowFreq True if a low-freq cal should be run
   * @return Data from running the experiment (plots and fit pole/zero values)
   * @throws IOException If a string does not refer to a valid accessible file
   */
  public RandData runRand(String calFileName, String outFileName,
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
      ir = new InstrumentResponse(respName);
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
   * Acquire data and run calibration over it.
   * Returns the experiment (all data kept locally to maintain thread safety)
   *
   * @param calFileName Filename of calibration signal
   * @param outFileName Filename of sensor output
   * @param respName Filename of response to load in
   * @param useEmbeddedResp True if response is an embedded response in program
   * @param startDate Long representing ms-since-epoch of data start time
   * @param endDate Long representing ms-since-epoch of data end time
   * @param lowFreq True if a low-freq cal should be run
   * @return Data from running the experiment (plots and fit pole/zero values)
   * @throws IOException If a string does not refer to a valid accessible file
   */
  public RandData populateDataAndRun(String calFileName, String outFileName,
      String respName, boolean useEmbeddedResp, String startDate, String endDate, boolean lowFreq)
      throws IOException, SeedFormatException, CodecException {
    return
        runRand(calFileName, outFileName, respName, useEmbeddedResp, startDate, endDate, lowFreq);
  }

  /**
   * Acquire data and run calibration over it. Used to handle calibrations that cross day boundaries
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
   */
  public RandData runRand(String calFileNameD1, String calFileNameD2,
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
      ir = new InstrumentResponse(respName);
    }

    ds.setBlock(0, calBlock);
    ds.setBlock(1, outBlock);
    ds.setResponse(1, ir);
    ds.trim(start, end);

    return runExpGetDataRand(ds, lowFreq);
  }

  public CalResult runStep(String calFileName, String outFileName,
      String respName, boolean useEmbeddedResp, String startDate, String endDate)
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
      ir = new InstrumentResponse(respName);
    }

    ds.setBlock(0, calBlock);
    ds.setBlock(1, outBlock);
    ds.setResponse(1, ir);
    ds.trim(start, end);

    return runExpGetDataStep(ds);

  }

  public CalResult runSine(String calFileName, String outFileName,
      String respName, boolean useEmbeddedResp, String startDate, String endDate)
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
      ir = new InstrumentResponse(respName);
    }

    ds.setBlock(0, calBlock);
    ds.setBlock(1, outBlock);
    ds.setResponse(1, ir);
    ds.trim(start, end);

    return runExpGetDataSine(ds);

  }

  @SuppressWarnings("unused")
  private CalResult runExpGetDataSine(DataStore ds) throws IOException {
    SineExperiment sine = new SineExperiment();
    sine.runExperimentOnData(ds);
    List<XYSeriesCollection> plots = sine.getData();
    double calAmplitude = sine.getCalAmplitude();
    double outAmplitude = sine.getOutAmplitude();
    double estFreq = sine.getEstSineFreq();
    double ratio = calAmplitude / outAmplitude;

    JFreeChart sineChart = ChartFactory.createXYLineChart(
        "Sine Calibration",
        "Time from data start (s)",
        "Normalized calibration signals (counts)",
        plots.get(0));

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
   return new SineData(pngByteArrays, calAmplitude, outAmplitude, estFreq, ratio);
  }

  @SuppressWarnings("unused")
  private CalResult runExpGetDataStep(DataStore ds) throws IOException {
    StepExperiment step = new StepExperiment();
    step.runExperimentOnData(ds);
    double[] fitParams = step.getFitParams();
    double[] initParams = step.getInitParams();
    List<XYSeriesCollection> plots = step.getData();
    // order of plots -- step function, resp amplitudes, resp phases

    NumberAxis stepAxis = new NumberAxis("Step counts");
    NumberAxis timeAxis = new NumberAxis("Time from data start (s)");
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
    return new StepData(pngByteArrays, fitParams, initParams);
  }

  /**
   * Acquire data and run calibration over it. Used to handle calibrations that cross day boundaries
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
   */
  @SuppressWarnings("unused")
  public RandData populateDataAndRun(String calFileNameD1, String calFileNameD2,
      String outFileNameD1, String outFileNameD2, String respName, boolean useEmbeddedResp,
      String startDate, String endDate, boolean lowFreq)
      throws IOException, SeedFormatException, CodecException {

    return runRand(calFileNameD1, calFileNameD2, outFileNameD1, outFileNameD2, respName,
        useEmbeddedResp, startDate, endDate, lowFreq);

  }

  private RandData runExpGetDataRand(DataStore dataStore, boolean isLowFrequency) throws IOException {

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
    ValueAxis amplitudeAxis = new NumberAxis(amplitudeAxisTitle);
    amplitudeAxis.setAutoRange(true);
    ((NumberAxis) amplitudeAxis).setAutoRangeIncludesZero(false);
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

    Map<String, List<Pair<Date, Date>>> gaps = randomExperiment.getGapRegions();
    String[] names = gaps.keySet().toArray(new String[]{});
    Date[][] gapStarts = new Date[names.length][];
    Date[][] gapEnds = new Date[names.length][];
    for (int j = 0; j < names.length; ++j) {
      String name = names[j];
      List<Pair<Date, Date>> dates = gaps.get(name);
      Date[] starts = new Date[dates.size()];
      Date[] ends = new Date[dates.size()];
      for (int i = 0; i < dates.size(); ++i) {
        starts[i] = dates.get(i).getFirst();
        ends[i] = dates.get(i).getSecond();
      }
      gapStarts[j] = starts;
      gapEnds[j] = ends;
    }

    return new RandData(fitPoles, fitZeros, initialPoles, initialZeros, pngByteArrays,
        names, gapStarts, gapEnds);

  }

}
