package asl.sensor;

import asl.sensor.experiment.RandomizedExperiment;
import asl.sensor.gui.ExperimentPanel;
import asl.sensor.gui.RandomizedPanel;
import asl.sensor.input.DataBlock;
import asl.sensor.input.DataStore;
import asl.sensor.input.InstrumentResponse;
import asl.sensor.utils.ReportingUtils;
import asl.sensor.utils.TimeSeriesUtils;
import java.awt.BasicStroke;
import java.awt.Font;
import java.awt.image.BufferedImage;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;
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
import py4j.GatewayServer;
import py4j.Py4JNetworkException;

/**
 * @author akearns
 */
public class CalProcessingServer {

  public CalProcessingServer() {
  }
  
  /**
   * Acquire data and run calibration over it. Used to handle calibrations that cross day boundaries
   * Returns the experiment (all data kept locally to maintain thread safety)
   * @param calFileNameD1 Filename of calibration signal (day 1)
   * @param calFileNameD2 Filename of calibration signal (day 2)
   * @param outFileNameD1 Filename of sensor output (day 1)
   * @param outFileNameD2 Filename of sensor output (day 2)
   * @param respName Filename of response to load in
   * @param respEmbd True if response is an embedded response in program
   * @param startTime Long representing ms-since-epoch of data start time
   * @param endTime Long representing ms-since-epoch of data end time
   * @param lowFreq True if a low-freq cal should be run
   * @return Data from running the experiment (plots and fit pole/zero values)
   * @throws IOException If a string does not refer to a valid accessible file
   */
  public RandData populateDataAndRun(String calFileNameD1, String calFileNameD2, 
      String outFileNameD1, String outFileNameD2, String respName, boolean respEmbd, long startTime,
      long endTime, boolean lowFreq) throws IOException {
    
      DataStore ds = new DataStore();
      String[] calFileName = new String[]{calFileNameD1, calFileNameD2};
      String[] outFileName = new String[]{outFileNameD1, outFileNameD2};
      DataBlock calBlock = TimeSeriesUtils.getFirstTimeSeries(calFileName);
      DataBlock outBlock = TimeSeriesUtils.getFirstTimeSeries(outFileName);
      InstrumentResponse ir;
      if (respEmbd) {
        ir = InstrumentResponse.loadEmbeddedResponse(respName);
      } else{
        ir = new InstrumentResponse(respName);
      }

      ds.setBlock(0, calBlock);
      ds.setBlock(1, outBlock);
      ds.setResponse(1, ir);
      ds.trim(startTime, endTime);
      
      return runExpGetData(ds, lowFreq);
    
  }

  /**
   * Acquire data and run calibration over it.
   * Returns the experiment (all data kept locally to maintain thread safety)
   * @param calFileName Filename of calibration signal
   * @param outFileName Filename of sensor output 
   * @param respName Filename of response to load in
   * @param respEmbd True if response is an embedded response in program
   * @param startTime Long representing ms-since-epoch of data start time
   * @param endTime Long representing ms-since-epoch of data end time
   * @param lowFreq True if a low-freq cal should be run
   * @return Data from running the experiment (plots and fit pole/zero values)
   * @throws IOException If a string does not refer to a valid accessible file
   */
  public RandData populateDataAndRun(String calFileName, String outFileName, 
      String respName, boolean respEmbd, long startTime, long endTime, boolean lowFreq) 
      throws IOException {
    
      DataStore ds = new DataStore();
      DataBlock calBlock = TimeSeriesUtils.getFirstTimeSeries(calFileName);
      DataBlock outBlock = TimeSeriesUtils.getFirstTimeSeries(outFileName);
      InstrumentResponse ir;
      if (respEmbd) {
        ir = InstrumentResponse.loadEmbeddedResponse(respName);
      } else{
        ir = new InstrumentResponse(respName);
      }

      ds.setBlock(0, calBlock);
      ds.setBlock(1, outBlock);
      ds.setResponse(1, ir);
      ds.trim(startTime, endTime);

      return runExpGetData(ds, lowFreq);

  }
  
  private RandData runExpGetData(DataStore ds, boolean lowFreq) throws IOException {
    
    RandomizedExperiment re = new RandomizedExperiment();
    
    re.setLowFreq(lowFreq);
    re.runExperimentOnData(ds);
    
    Complex[] fitZerosCpx = re.getFitResponse().getZeros().toArray(new Complex[]{});
    Complex[] fitPolesCpx = re.getFitResponse().getPoles().toArray(new Complex[]{});
    Complex[] initZerosCpx = ds.getResponse(1).getZeros().toArray(new Complex[]{});
    Complex[] initPolesCpx = ds.getResponse(1).getPoles().toArray(new Complex[]{});
    
    double[] zeros = new double[2 * fitZerosCpx.length];
    double[] initZeros = new double[zeros.length];
    double[] poles = new double[2 * fitPolesCpx.length];
    double[] initPoles = new double[poles.length];
    for (int i = 0; i < fitZerosCpx.length; ++i) {
      int reIdx = 2 * i; int imIdx = reIdx + 1;
      zeros[reIdx] = fitZerosCpx[i].getReal();
      zeros[imIdx] = fitZerosCpx[i].getImaginary();
      initZeros[reIdx] = initZerosCpx[i].getReal();
      initZeros[imIdx] = initZerosCpx[i].getImaginary();
    }
    for(int i = 0; i < fitPolesCpx.length; ++i) {
      int reIdx = 2 * i; int imIdx = reIdx + 1;
      poles[reIdx] = fitPolesCpx[i].getReal();
      poles[imIdx] = fitPolesCpx[i].getImaginary();
      initPoles[reIdx] = initPolesCpx[i].getReal();
      initPoles[imIdx] = initPolesCpx[i].getImaginary();
    }
    
    List<XYSeriesCollection> xyscList = re.getData();
    JFreeChart[] chartArr = new JFreeChart[xyscList.size()];
    
    String xAxisTitle = "Frequency (Hz)";
    String ampAxisTitle = "10 * log10( RESP(f) )";
    String phaseAxisTitle = "phi(RESP(f))";
    
    ValueAxis xAxis = new LogarithmicAxis(xAxisTitle);
    ValueAxis residXAxis = new LogarithmicAxis(xAxisTitle);
    ValueAxis ampAxis = new NumberAxis(ampAxisTitle);
    ampAxis.setAutoRange(true);
    ( (NumberAxis) ampAxis).setAutoRangeIncludesZero(false);
    ValueAxis phaseAxis = new NumberAxis(phaseAxisTitle);
    phaseAxis.setAutoRange(true);
    ValueAxis residPhaseAxis = new NumberAxis("Phase error (degrees)");
    ValueAxis residAmpAxis = new NumberAxis("Amplitude error (percentage)");
    Font bold = xAxis.getLabelFont().deriveFont(Font.BOLD);
    xAxis.setLabelFont(bold);
    ampAxis.setLabelFont(bold);
    phaseAxis.setLabelFont(bold);
    residXAxis.setLabelFont(bold);
    residPhaseAxis.setLabelFont(bold);
    residAmpAxis.setLabelFont(bold);
    XYPlot xyp;
    
    StringBuilder title = new StringBuilder();
    if (lowFreq) {
      title.append("Low-freq random cal");
    } else {
      title.append("High-freq random cal");
    }
    
    chartArr[0] = ChartFactory.createXYLineChart(
        title.toString(),
        xAxis.getLabel(),
        ampAxis.getLabel(),
        xyscList.get(0),
        PlotOrientation.VERTICAL,
        true, // include legend
        false, 
        false);
    xyp = chartArr[0].getXYPlot();
    xyp.setDomainAxis(xAxis);
    xyp.setRangeAxis(ampAxis);
    ExperimentPanel.invertSeriesRenderingOrder(chartArr[0]);
    
    chartArr[1] = ChartFactory.createXYLineChart(
        title.toString(),
        xAxis.getLabel(),
        phaseAxis.getLabel(),
        xyscList.get(1),
        PlotOrientation.VERTICAL,
        true, // include legend
        false, 
        false);
    xyp = chartArr[1].getXYPlot();
    xyp.setDomainAxis(xAxis);
    xyp.setRangeAxis(phaseAxis);
    ExperimentPanel.invertSeriesRenderingOrder(chartArr[1]);
    
    chartArr[2] = ChartFactory.createXYLineChart(
        title.toString(),
        residXAxis.getLabel(),
        residAmpAxis.getLabel(),
        xyscList.get(2),
        PlotOrientation.VERTICAL,
        true, // include legend
        false, 
        false);
    xyp = chartArr[2].getXYPlot();
    xyp.setDomainAxis(residXAxis);
    xyp.setRangeAxis(residAmpAxis);
    ExperimentPanel.invertSeriesRenderingOrder(chartArr[2]);
    
    chartArr[3] = ChartFactory.createXYLineChart(
        title.toString(),
        residXAxis.getLabel(),
        residPhaseAxis.getLabel(),
        xyscList.get(3),
        PlotOrientation.VERTICAL,
        true, // include legend
        false, 
        false);
    xyp = chartArr[3].getXYPlot();
    xyp.setDomainAxis(residXAxis);
    xyp.setRangeAxis(residPhaseAxis);
    ExperimentPanel.invertSeriesRenderingOrder(chartArr[3]);
    
    if (!lowFreq) {
      Marker maxFitMarker = new ValueMarker( re.getMaxFitFrequency() );
      maxFitMarker.setStroke( new BasicStroke( (float) 1.5 ) );
      chartArr[0].getXYPlot().addDomainMarker(maxFitMarker);
      chartArr[1].getXYPlot().addDomainMarker(maxFitMarker);
    }
    
    BufferedImage[] bi = ReportingUtils.chartsToImageList(1, 1280, 960, chartArr);
    byte[][] pngByteArrays = new byte[bi.length][];
    for (int i = 0; i < bi.length; ++i) {
      ByteArrayOutputStream out = new ByteArrayOutputStream();
      ImageIO.write(bi[i], "png", out);
      pngByteArrays[i] = out.toByteArray();
    }
    
    Map<String, List<Pair<Date, Date>>> gaps = re.getGapRegions();
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
    
    return new RandData(poles, zeros, initPoles, initZeros, pngByteArrays, 
        names, gapStarts, gapEnds);
    
  }
  
  /**
   * get all metadata from the function in a single file
   * @param exp
   * @return
   */
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
    } catch (Py4JNetworkException e){
      System.out.println("Already Running: Closing process");
      System.exit(0);
    }
    System.out.println("Gateway Server Started");
  }

  public class RandData {
    
    private double[] initPoles;
    private double[] initZeros;
    private double[] fitPoles;
    private double[] fitZeros;
    private byte[][] pngs;
    private String[] gapNameIdentifiers;
    private Date[][] gapStarts;
    private Date[][] gapEnds;
    
    public RandData(double[] fp, double[] fz, double[] ip, double[] iz, byte[][] im, 
        String[] nm, Date[][] gpa, Date[][] gpb) {
      fitPoles = fp;
      fitZeros = fz;
      initPoles = ip;
      initZeros = iz;
      pngs = im;
      gapNameIdentifiers = nm;
      gapStarts = gpa;
      gapEnds = gpb;
    }
    
    public double[] getInitPoles() {
      return initPoles;
    }
    
    public double[] getInitZeros() {
      return initZeros;
    }
    
    public double[] getFitPoles() {
      return fitPoles;
    }
    
    public double[] getFitZeros() {
      return fitZeros;
    }
    
    public byte[] getAmpImage() {
      return pngs[0];
    }
    
    public byte[] getPhaseImage() {
      return pngs[1];
    }
    
    public byte[] getAmpErrorImage() {
      return pngs[2];
    }
    
    public byte[] getPhaseErrorImage() {
      return pngs[3];
    }
    
    public String[] getGapIdentifiers() {
      return gapNameIdentifiers;
    }
    
    public Date[][] getGapStartDates() {
      return gapStarts;
    }
    
    public Date[][] getGapEndDates() {
      return gapEnds;
    }
    
    public String getGapInfoAsString() {
      SimpleDateFormat sdf = new SimpleDateFormat("DD.HH:m:s");
      sdf.setTimeZone( TimeZone.getTimeZone("UTC") );
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
          sb.append( df.format(start) );
          sb.append("\t");
          sb.append( df.format(end) );
          sb.append("\n");
        }
        sb.append("\n");
      }
      return sb.toString();
    }
  }
  
}
