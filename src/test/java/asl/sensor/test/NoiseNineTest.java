package asl.sensor.test;

import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.awt.Font;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.List;
import java.util.TimeZone;

import org.apache.pdfbox.pdmodel.PDDocument;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.LogarithmicAxis;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.data.xy.XYSeriesCollection;
import org.junit.Test;

import asl.sensor.experiment.ExperimentEnum;
import asl.sensor.experiment.NoiseNineExperiment;
import asl.sensor.gui.InputPanel;
import asl.sensor.gui.NoiseNinePanel;
import asl.sensor.input.DataStore;
import asl.sensor.utils.ReportingUtils;
import asl.sensor.utils.TimeSeriesUtils;

public class NoiseNineTest {

  String currentDir = System.getProperty("user.dir");

  @Test
  public void canRunAndPlotTest1() {

    String folder = currentDir + "/data/noisenine/";
    String[] types = new String[]{"00","10","60"};;
    String freqName = "_BH";
    String[] components = new String[]{"1","2","Z"};
    String ending = ".512.seed";
    String respName = "STS-1_Q330HR_BH_20";
    
    DataStore ds = 
        setUpTest(folder, types, freqName, components, ending, respName, true);

    SimpleDateFormat sdf = InputPanel.SDF;
    sdf.setTimeZone( TimeZone.getTimeZone("UTC") );
    // sdf.setLenient(false);

    Calendar cCal = Calendar.getInstance( sdf.getTimeZone() );
    cCal.setTimeInMillis( ds.getBlock(0).getStartTime() );
    cCal.set(Calendar.HOUR, 7);
    System.out.println( "start: " + sdf.format( cCal.getTime() ) );
    long start = cCal.getTime().getTime();
    cCal.set(Calendar.HOUR, 8);
    cCal.set(Calendar.MINUTE, 00);
    System.out.println( "end: " + sdf.format( cCal.getTime() ) );
    long end = cCal.getTime().getTime();

    ds.trim(start, end);
    // ds.matchIntervals();

    // trim down internal representation of the data to avoid memory errors
    /*
    for (int i = 0; i < 9; ++i) {
      DataBlock db = ds.getBlock(i);
      db.setData( db.getData() );
      System.out.println( db.size() );
    }
    */
    
    NoiseNineExperiment nne = new NoiseNineExperiment();
    assertTrue( nne.hasEnoughData(ds) );
    nne.runExperimentOnData(ds);

    List<XYSeriesCollection> xysc = nne.getData();
    JFreeChart[] jfcl = new JFreeChart[xysc.size()];

    String xAxisTitle = "Frequency (Hz)";
    NumberAxis xAxis = new LogarithmicAxis(xAxisTitle);
    Font bold = xAxis.getLabelFont().deriveFont(Font.BOLD);
    xAxis.setLabelFont(bold);
    String yAxisTitle = "Power (rel. 1 (m/s^2)^2/Hz)";
    String[] orientations = new String[]{" (North)", " (East)", " (Vertical)"};

    for (int i = 0; i < xysc.size(); ++i) {
      jfcl[i] = ChartFactory.createXYLineChart(
          ExperimentEnum.RANDM.getName() + orientations[i],
          xAxisTitle,
          yAxisTitle,
          xysc.get(i),
          PlotOrientation.VERTICAL,
          true,
          false,
          false);

      XYPlot xyp = jfcl[i].getXYPlot();

      //xyp.clearAnnotations();
      //xyp.addAnnotation(xyt);

      xyp.setDomainAxis( xAxis );
    }

    StringBuilder sb = new StringBuilder();
    String insets = NoiseNinePanel.getInsetString(nne);
    sb.append('\n');
    sb.append( NoiseNinePanel.getTimeStampString(nne) );
    sb.append('\n');
    sb.append("INPUTTED FILES:");
    sb.append('\n');
    
    List<String> names = nne.getInputNames();
    
    for (String name : names) {
      sb.append( name );
      sb.append('\n');
    }
    

    int width = 1280; int height = 960;

    PDDocument pdf = new PDDocument();
    ReportingUtils.chartsToPDFPage(width, height, pdf, jfcl);
    ReportingUtils.textListToPDFPages( pdf, insets, sb.toString() );

    String testResultFolder = currentDir + "/testResultImages/";
    File dir = new File(testResultFolder);
    if ( !dir.exists() ) {
      dir.mkdir();
    }

    String testResult = testResultFolder + "Nine-Noise-Test.pdf";
    try {
      pdf.save( new File(testResult) );
      pdf.close();
    } catch (IOException e) {
      e.printStackTrace();
      fail();
    }
    
    System.out.println("Output result has been written");

  }
  
  @Test
  public void canRunAndPlotTest2() {
    String folder = currentDir + "/data/noise_test/";
    String[] types = new String[]{"00","10","30"};
    String freqName = "_BH";
    String[] components = new String[]{"1","2","Z"};
    String ending = ".512.seed";
    String respName = "responses/RESP.XX.MOFO.00.BHZ";
    
    DataStore ds = 
        setUpTest(folder, types, freqName, components, ending, respName, false);

    SimpleDateFormat sdf = InputPanel.SDF;
    sdf.setTimeZone( TimeZone.getTimeZone("UTC") );
    // sdf.setLenient(false);

    Calendar cCal = Calendar.getInstance( sdf.getTimeZone() );
    cCal.setTimeInMillis( ds.getBlock(0).getStartTime() );
    cCal.set(Calendar.HOUR_OF_DAY, 12);
    cCal.set(Calendar.MINUTE, 0);
    cCal.set(Calendar.SECOND, 0);
    System.out.println( "start: " + sdf.format( cCal.getTime() ) );
    long start = cCal.getTime().getTime();
    cCal.set(Calendar.HOUR_OF_DAY, 13);
    cCal.set(Calendar.MINUTE, 00);
    cCal.set(Calendar.SECOND, 0);
    System.out.println( "end: " + sdf.format( cCal.getTime() ) );
    long end = cCal.getTime().getTime();
    
    ds.trim(start, end, DataStore.FILE_COUNT);
    // ds.matchIntervals();

    // trim down internal representation of the data to avoid memory errors
    /*
    for (int i = 0; i < 9; ++i) {
      DataBlock db = ds.getBlock(i);
      db.setData( db.getData() );
    }
    */
    
    NoiseNineExperiment nne = new NoiseNineExperiment();
    assertTrue( nne.hasEnoughData(ds) );
    nne.runExperimentOnData(ds);

    List<XYSeriesCollection> xysc = nne.getData();
    JFreeChart[] jfcl = new JFreeChart[xysc.size()];

    String xAxisTitle = "Frequency (Hz)";
    NumberAxis xAxis = new LogarithmicAxis(xAxisTitle);
    Font bold = xAxis.getLabelFont().deriveFont(Font.BOLD);
    xAxis.setLabelFont(bold);
    String yAxisTitle = "Power (rel. 1 (m/s^2)^2/Hz)";
    String[] orientations = new String[]{" (North)", " (East)", " (Vertical)"};
    
    for (int i = 0; i < xysc.size(); ++i) {
      
      // make sure each series is populated with data from correct orientation
      // i.e., BH1s for north, BH2s for south, etc.
      String check = freqName + components[i];
      XYSeriesCollection coll = xysc.get(i);
      // has 3 components, so do 3 checks
      // nlnm, nhnm should be 4th and 5th entries in data
      for (int j = 0; j < 3; ++j) {
        String key = (String) coll.getSeriesKey(j);
        assertTrue(key.contains(check));
      }
      
      jfcl[i] = ChartFactory.createXYLineChart(
          ExperimentEnum.RANDM.getName() + orientations[i],
          xAxisTitle,
          yAxisTitle,
          coll,
          PlotOrientation.VERTICAL,
          true,
          false,
          false);

      XYPlot xyp = jfcl[i].getXYPlot();

      //xyp.clearAnnotations();
      //xyp.addAnnotation(xyt);

      xyp.setDomainAxis( xAxis );
    }

    String insets = NoiseNinePanel.getInsetString(nne);
    StringBuilder sb = new StringBuilder();
    // sb.append( NoiseNinePanel.getInsetString(nne) );
    // sb.append('\n');
    sb.append( NoiseNinePanel.getTimeStampString(nne) );
    sb.append('\n');
    sb.append("INPUTTED FILES:");
    sb.append('\n');
    
    List<String> names = nne.getInputNames();
    
    for (String name : names) {
      sb.append( name );
      sb.append('\n');
    }
    

    int width = 1280; int height = 960;

    PDDocument pdf = new PDDocument();
    ReportingUtils.chartsToPDFPage(width, height, pdf, jfcl);
    ReportingUtils.textListToPDFPages( pdf, insets, sb.toString() );

    String testResultFolder = currentDir + "/testResultImages/";
    File dir = new File(testResultFolder);
    if ( !dir.exists() ) {
      dir.mkdir();
    }

    String testResult = testResultFolder + "Nine-Noise-Test-2.pdf";
    try {
      pdf.save( new File(testResult) );
      pdf.close();
    } catch (IOException e) {
      e.printStackTrace();
      fail();
    }
    
    System.out.println("Output result has been written");

  }
  
  public DataStore setUpTest(String folder, String[] types, String freqName, 
      String[] components, String ending, String respName, boolean isEmbed) {
    
    DataStore ds = new DataStore();

    for (int i = 0; i < types.length; ++i) {
      for (int j = 0; j < components.length; ++j) {
        int indexInStore = i * types.length + j;

        String fName = folder + types[i] + freqName + components[j] + ending;

        System.out.println(fName);

        List<String> dataNames = new ArrayList<String>();
        try {
          dataNames = new ArrayList<String>( 
              TimeSeriesUtils.getMplexNameSet(fName) );
        } catch (FileNotFoundException e) {
          fail();
          e.printStackTrace();
        }
        ds.setBlock(indexInStore, fName, dataNames.get(0) );
        if (isEmbed) {
          ds.setEmbedResponse(indexInStore, respName);
        } else {
          ds.setResponse(indexInStore, respName);
        }
      }
    }
    
    return ds;
    
  }

}
