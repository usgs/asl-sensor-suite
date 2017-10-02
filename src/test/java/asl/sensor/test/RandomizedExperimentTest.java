package asl.sensor.test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.awt.Font;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.List;
import java.util.TimeZone;

import javax.imageio.ImageIO;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.linear.RealVector;
import org.apache.pdfbox.pdmodel.PDDocument;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.LogarithmicAxis;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.junit.Test;

import asl.sensor.experiment.ExperimentEnum;
import asl.sensor.experiment.ExperimentFactory;
import asl.sensor.experiment.RandomizedExperiment;
import asl.sensor.gui.InputPanel;
import asl.sensor.gui.RandomizedPanel;
import asl.sensor.input.DataBlock;
import asl.sensor.input.DataStore;
import asl.sensor.input.InstrumentResponse;
import asl.sensor.utils.FFTResult;
import asl.sensor.utils.NumericUtils;
import asl.sensor.utils.ReportingUtils;
import asl.sensor.utils.TimeSeriesUtils;

public class RandomizedExperimentTest {

  
  public DataStore getFromList(List<String> setUpFilenames) throws IOException {
    
    String respName = setUpFilenames.get(0);
    String calName = setUpFilenames.get(1);
    String sensOutName = setUpFilenames.get(2);
    String metaName;
    
    System.out.println(respName);
    System.out.println(calName);
    System.out.println(sensOutName);
    
    metaName = TimeSeriesUtils.getMplexNameList(calName).get(0);
    DataBlock calib = TimeSeriesUtils.getTimeSeries(calName, metaName);
    
    InstrumentResponse ir = new InstrumentResponse(respName);
    
    metaName = TimeSeriesUtils.getMplexNameList(sensOutName).get(0);
    DataBlock sensor = TimeSeriesUtils.getTimeSeries(sensOutName, metaName);

    DataStore ds = new DataStore();
    ds.setBlock(0, calib);
    ds.setBlock(1, sensor);
    ds.setResponse(1, ir);
    
    return ds;
    
  }
  
  public Calendar getStartCalendar(DataStore ds) {
    SimpleDateFormat sdf = InputPanel.SDF;
    sdf.setTimeZone( TimeZone.getTimeZone("UTC") );
    Calendar cCal = Calendar.getInstance( sdf.getTimeZone() );
    
    cCal.setTimeInMillis( ds.getBlock(0).getStartTime() );
    return cCal;
  }
  
  @Test
  public void ResponseCorrectConvertedToVectorHighFreq() {
    String fname = "responses/TST5_response.txt";
    boolean lowFreq = false;
    InstrumentResponse ir;
    try {
      
      ir = new InstrumentResponse(fname);
      List<Complex> poles = new ArrayList<Complex>( ir.getPoles() );
      // using an unnecessarily high nyquist rate here
      RealVector high = ir.polesToVector(lowFreq, 1E8);
      
      int complexIndex = 2; // start at second pole
      int vectorIndex = 0;
      
      while ( vectorIndex < high.getDimension() ) {
        // return current index
        double real = high.getEntry(vectorIndex++);
        double imag = high.getEntry(vectorIndex++);
        
        double poleImag = poles.get(complexIndex).getImaginary();
        
        assertEquals( real, poles.get(complexIndex).getReal(), 0.0 );
        assertEquals( imag, poleImag, 0.0 );
        
        if (poleImag != 0) {
          // complex conjugate case
          ++complexIndex;
          assertEquals( real, poles.get(complexIndex).getReal(), 0.0 );
          assertEquals( imag, -poles.get(complexIndex).getImaginary(), 0.0 );
        }
        
        ++complexIndex;
        
      }
      
    } catch (IOException e) {
      fail();
      e.printStackTrace();
    }
  }
  
  @Test
  public void ResponseCorrectlyConvertedToVectorLowFreq() {
    String fname = "responses/TST5_response.txt";
    boolean lowFreq = true;
    InstrumentResponse ir;
    try {
      
      ir = new InstrumentResponse(fname);
      List<Complex> poles = new ArrayList<Complex>( ir.getPoles() );
      // again, use a very high nyquist rate
      RealVector low = ir.polesToVector(lowFreq, 1E8);
      
      // only test lower two poless
      assertEquals( low.getEntry(0), poles.get(0).getReal(), 0.0 );
      assertEquals( low.getEntry(1), poles.get(0).getImaginary(), 0.0 );
      
      assertEquals( low.getEntry(0), poles.get(1).getReal(), 0.0 );
      assertEquals( low.getEntry(1), -poles.get(1).getImaginary(), 0.0 );
      
    } catch (IOException e) {
      fail();
      e.printStackTrace();
    }
  }
  
  @Test
  public void ResponseSetCorrectlyHighFreq() {
    
    String fname = "responses/TST5_response.txt";
    
    InstrumentResponse ir;

      try {
        ir = new InstrumentResponse(fname);
        boolean lowFreq = false;
        
        List<Complex> poles = new ArrayList<Complex>( ir.getPoles() );
        List<Complex> replacements = new ArrayList<Complex>();
        
        int start = 2;
        if ( poles.get(0).getImaginary() == 0 ) {
          start = 1;
        }
        
        for (int i = start; i < poles.size(); ++i) {
          if ( poles.get(i).getImaginary() == 0 ) {
            Complex c = poles.get(i);
            replacements.add(c.subtract(1));
            int next = i+1;
            while (next < poles.size() && poles.get(next).equals(c)) {
              ++next; // skip duplicates
            }
          } else {
            Complex c = poles.get(i);
            c = c.subtract( new Complex(1, 1) );
            replacements.add(c);
            ++i;
          }
        }
        
        //System.out.println(poles);
        //System.out.println(replacements);
        
        double[] newPoles = new double[replacements.size() * 2];
        for (int i = 0; i < newPoles.length; i += 2) {
          int poleIdx = i / 2;
          Complex c = replacements.get(poleIdx);
          newPoles[i] = c.getReal();
          newPoles[i + 1] = c.getImaginary();
        }
        
        InstrumentResponse ir2 = 
            ir.buildResponseFromFitVector(newPoles, lowFreq, 0);
        
        List<Complex> testList = ir2.getPoles();
        //System.out.println(testList);
        int offsetIdx = 0;
        for (int i = 0; i < poles.size(); ++i) {
          if (i < start) {
            assertTrue( poles.get(i).equals( testList.get(i) ) );
          } else {
            Complex c = replacements.get(offsetIdx);
            assertTrue( testList.get(i).equals(c) );
            if ( poles.get(i).getImaginary() != 0 ) {
              Complex c1 = new Complex(1, 1);
              assertTrue(poles.get(i).equals( c.add(c1) ));
              ++i;
              Complex c2 = new Complex(1, -1);
              assertTrue( testList.get(i).equals( c.conjugate() ) );
              assertTrue( poles.get(i).equals( c.conjugate().add(c2) ) );
            } else {
              assertTrue( poles.get(i).equals(c.add(1)) );
            }
            ++offsetIdx;
          }
        }
        
      } catch (IOException e) {
        e.printStackTrace();
      }
    
  }
  
  @Test
  public void testSingleSideFFTValues() {
    // just to test to see how crosspower compares to single-side fft
    String folder = "data/highfrq-majo/";
    String calInFile = "CB_BC1.512.seed";
    String sensorOutFile = "10_EHZ.512.seed";
    try {
      DataBlock cal = TimeSeriesUtils.getFirstTimeSeries(folder + calInFile);
      DataBlock out = TimeSeriesUtils.getFirstTimeSeries(folder + sensorOutFile);
      
      Calendar start = cal.getTrimmedStartCalendar();
      Calendar end = cal.getTrimmedStartCalendar();
      start.set(Calendar.HOUR_OF_DAY, 18);
      start.set(Calendar.MINUTE, 20);
      end.set(Calendar.HOUR_OF_DAY, 18);
      end.set(Calendar.MINUTE, 35);
      cal.trim(start, end);
      out.trim(start, end);
      
      FFTResult outputSingleSide = FFTResult.singleSidedFFT(out, false);
      FFTResult calSingleSide = FFTResult.singleSidedFFT(cal, false);
      
      Complex[] outputFFT = outputSingleSide.getFFT();
      Complex[] calFFT = calSingleSide.getFFT();
      double[] freqs = outputSingleSide.getFreqs();
      
      double minFreq = 30;
      
      XYSeries numXYS = new XYSeries("Nom FFT amp [out x cal]");
      XYSeries denXYS = new XYSeries("Dnm FFT amp [cal x cal]");
      double nyquist = cal.getSampleRate() / 2;
      for (int i = 0; freqs[i] < .8 * nyquist; ++i) {
        if (freqs[i] < minFreq) {
          continue;
        }
        // Complex scaleFreq = new Complex(0., NumericUtils.TAU * freqs[i]);
        double dBNumer = 10 * Math.log10(outputFFT[i].abs());
        double dBDenom = 10 * Math.log10(calFFT[i].abs());
        numXYS.add(freqs[i], dBNumer);
        denXYS.add(freqs[i], dBDenom);
      }
      
      XYSeriesCollection xysc = new XYSeriesCollection();
      xysc.addSeries(numXYS);
      xysc.addSeries(denXYS);
      JFreeChart chart = ChartFactory.createXYLineChart(
          "High-freq. Cal Deconvolution verification",
          "Frequency",
          "Spectrum magnitude (log scale)",
          xysc,
          PlotOrientation.VERTICAL,
          true,
          false,
          false);
      LogarithmicAxis la = new LogarithmicAxis("Freq. (Hz)");
      chart.getXYPlot().setDomainAxis(la);
      
      String pngFname = "testResultImages/MAJO-plot-singleside-rcal.png";
      int width = 640; int height = 480;
      BufferedImage bi = 
          ReportingUtils.chartsToImage(width, height, chart);
      File file = new File(pngFname);
      ImageIO.write(bi, "png", file);
      
      
    } catch (FileNotFoundException e) {
      e.printStackTrace();
      fail();
    } catch (IOException e) {
      e.printStackTrace();
      fail();
    }
    
  }
  
  @Test
  public void testCrossPowerCalculations() {
    String folder = "data/highfrq-majo/";
    String calInFile = "CB_BC1.512.seed";
    String sensorOutFile = "10_EHZ.512.seed";
    try {
      DataBlock cal = TimeSeriesUtils.getFirstTimeSeries(folder + calInFile);
      DataBlock out = TimeSeriesUtils.getFirstTimeSeries(folder + sensorOutFile);
      
      Calendar start = cal.getTrimmedStartCalendar();
      Calendar end = cal.getTrimmedStartCalendar();
      start.set(Calendar.HOUR_OF_DAY, 18);
      start.set(Calendar.MINUTE, 20);
      end.set(Calendar.HOUR_OF_DAY, 18);
      end.set(Calendar.MINUTE, 35);
      cal.trim(start, end);
      out.trim(start, end);
      
      FFTResult numer = FFTResult.spectralCalc(out, cal);
      FFTResult denom = FFTResult.spectralCalc(cal, cal);
      
      Complex[] numerFFT = numer.getFFT();
      Complex[] denomFFT = denom.getFFT();
      double[] freqs = numer.getFreqs();
      
      double minFreq = 30;
      
      XYSeries numXYS = new XYSeries("Nom FFT amp [out x cal]");
      XYSeries denXYS = new XYSeries("Dnm FFT amp [cal x cal]");
      XYSeries dcvXYS = new XYSeries("Deconvolution (nom/dnm), scaled to vel.");
      double nyquist = cal.getSampleRate() / 2;
      for (int i = 0; freqs[i] < .8 * nyquist; ++i) {
        if (freqs[i] < minFreq) {
          continue;
        }
        Complex scaleFreq = new Complex(0., NumericUtils.TAU * freqs[i]);
        double dBNumer = 10 * Math.log10(numerFFT[i].abs());
        double dBDenom = 10 * Math.log10(denomFFT[i].abs());
        Complex div = numerFFT[i].divide(denomFFT[i].abs());
        Complex todB = div.multiply(scaleFreq);
        double dBDiv = 10 * Math.log10(todB.abs());
        numXYS.add(freqs[i], dBNumer);
        denXYS.add(freqs[i], dBDenom);
        dcvXYS.add(freqs[i], dBDiv);
      }
      
      XYSeriesCollection xysc = new XYSeriesCollection();
      xysc.addSeries(numXYS);
      xysc.addSeries(denXYS);
      JFreeChart chart = ChartFactory.createXYLineChart(
          "High-freq. Cal Deconvolution verification",
          "Frequency",
          "Spectrum magnitude (log scale)",
          xysc,
          PlotOrientation.VERTICAL,
          true,
          false,
          false);
      LogarithmicAxis la = new LogarithmicAxis("Freq. (Hz)");
      chart.getXYPlot().setDomainAxis(la);
      
      xysc = new XYSeriesCollection();
      xysc.addSeries(dcvXYS);
      JFreeChart chart2 = ChartFactory.createXYLineChart(
          "High-freq. Cal Deconvolution verification",
          "Frequency",
          "Spectrum magnitude (log scale)",
          xysc,
          PlotOrientation.VERTICAL,
          true,
          false,
          false);
      chart2.getXYPlot().setDomainAxis(la);      
      
      String pngFname = "testResultImages/MAJO-plot-highfreq-spectral.png";
      int width = 1280; int height = 480;
      BufferedImage bi = 
          ReportingUtils.chartsToImage(width, height, chart, chart2);
      File file = new File(pngFname);
      ImageIO.write(bi, "png", file);
      
      
    } catch (FileNotFoundException e) {
      e.printStackTrace();
      fail();
    } catch (IOException e) {
      e.printStackTrace();
      fail();
    }
    
  }
  
  @Test
  public void responseSetCorrectlyLowFreq() {
    
    String fname = "responses/TST5_response.txt";
    
      InstrumentResponse ir;
      try {
        ir = new InstrumentResponse(fname);
        boolean lowFreq = true;
        List<Complex> poles = new ArrayList<Complex>( ir.getPoles() );

        double[] newPoles = new double[2];
        newPoles[0] = 0.;
        newPoles[1] = 1.;
        
        Complex c = new Complex( newPoles[0], newPoles[1] );
        
        InstrumentResponse ir2 = 
            ir.buildResponseFromFitVector(newPoles, lowFreq, 0);
        List<Complex> poles2 = ir2.getPoles();

        List<Complex> testList = new ArrayList<Complex>(poles);
        testList.set(0, c);
        testList.set( 1, c.conjugate() );
        
        // System.out.println(testList);
        // System.out.println(poles);
        // System.out.println(poles2);
        
        for (int i = 0; i < poles.size(); ++i) {
          if (i < 2) {
            assertFalse( poles.get(i).equals( poles2.get(i) ) );
            assertTrue( poles2.get(i).equals( testList.get(i) ) );
          }
        }
        
        
      } catch (IOException e) {
        fail();
        e.printStackTrace();
      }
    
  }
  
  public DataStore setUpTest1() throws IOException {
    
    List<String> fileList = new ArrayList<String>();
    String respName = "responses/RESP.XX.NS088..BHZ.STS1.360.2400";
    String dataFolderName = "data/random_cal/"; 
    String calName =  dataFolderName + "_EC0.512.seed";
    String sensOutName = dataFolderName + "00_EHZ.512.seed";
    
    fileList.add(respName);
    fileList.add(calName);
    fileList.add(sensOutName);
    
    DataStore ds = getFromList(fileList);   
    Calendar cCal = getStartCalendar(ds);
    
    cCal.set(Calendar.MINUTE, 36);
    cCal.set(Calendar.SECOND, 0);
    long start = cCal.getTime().getTime();
    
    cCal.set(Calendar.MINUTE, 41);
    // System.out.println( "end: " + sdf.format( cCal.getTime() ) );
    long end = cCal.getTime().getTime();
    
    ds.trim(start, end);
    
    return ds;
  }

  public DataStore setUpTest2() throws IOException {
    List<String> fileList = new ArrayList<String>();
    String respName = "responses/RESP.XX.NS088..BHZ.STS1.360.2400";
    String dataFolderName = "data/random_cal_2/"; 
    String calName =  dataFolderName + "CB_BC0.512.seed";
    String sensOutName = dataFolderName + "00_EHZ.512.seed";
    
    fileList.add(respName);
    fileList.add(calName);
    fileList.add(sensOutName);
    
    DataStore ds = getFromList(fileList);
    
    // response we want is embedded
    InstrumentResponse ir;
    ir = InstrumentResponse.loadEmbeddedResponse("T-compact_Q330HR_BH_40");
    ds.setResponse(1, ir);
    
    Calendar cCal = getStartCalendar(ds);
    
    cCal.set(Calendar.MINUTE, 52);
    cCal.set(Calendar.SECOND, 0);
    long start = cCal.getTime().getTime();
    
    int hour = cCal.get(Calendar.HOUR);
    cCal.set(Calendar.HOUR, hour + 1);
    cCal.set(Calendar.MINUTE, 12);
    
    // System.out.println( "end: " + sdf.format( cCal.getTime() ) );
    long end = cCal.getTime().getTime();
    
    ds.trim(start, end);
    
    return ds;
  }
  
  public DataStore setUpTest3() throws IOException {
    List<String> fileList = new ArrayList<String>();
    String respName = "responses/RESP.XX.NS088..BHZ.STS1.360.2400";
    String dataFolderName = "data/random_cal_3/"; 
    String calName =  dataFolderName + "BC0.512.seed";
    String sensOutName = dataFolderName + "00_BHZ.512.seed";
    
    fileList.add(respName);
    fileList.add(calName);
    fileList.add(sensOutName);
    
    DataStore ds = getFromList(fileList);
    
    // response we want is embedded
    InstrumentResponse ir;
    ir = InstrumentResponse.loadEmbeddedResponse("KS54000_Q330HR_BH_40");
    ds.setResponse(1, ir);
    
    Calendar cCal = getStartCalendar(ds);
    
    cCal.set(Calendar.HOUR_OF_DAY, 21);
    cCal.set(Calendar.MINUTE, 24);
    cCal.set(Calendar.SECOND, 0);
    long start = cCal.getTime().getTime();
    
    // commented out -- calibration ends when the data does
    //int hour = cCal.get(Calendar.HOUR);
    /*
    cCal.set(Calendar.DAY_OF_YEAR, 4);
    cCal.set(Calendar.HOUR_OF_DAY, 0);
    cCal.set(Calendar.MINUTE, 0);
    */
    
    // System.out.println( "end: " + sdf.format( cCal.getTime() ) );
    long end = ds.getBlock(0).getEndTime();
    
    ds.trim(start, end);
    
    return ds;
  }
  
  public DataStore setUpTest4() throws IOException {
    
    List<String> fileList = new ArrayList<String>();
    String respName = "responses/RESP.XX.NS088..BHZ.STS1.360.2400";
    String dataFolderName = "data/random_cal_4/"; 
    String calName =  dataFolderName + "CB_BC0.512.seed";
    String sensOutName = dataFolderName + "00_EHZ.512.seed";
    
    fileList.add(respName);
    fileList.add(calName);
    fileList.add(sensOutName);
    
    DataStore ds = getFromList(fileList);
    
    // response we want is embedded
    InstrumentResponse ir;
    ir = InstrumentResponse.loadEmbeddedResponse("KS54000_Q330HR_BH_40");
    ds.setResponse(1, ir);
    
    Calendar cCal = getStartCalendar(ds);
    
    cCal.set(Calendar.HOUR_OF_DAY, 20);
    cCal.set(Calendar.MINUTE, 16);
    cCal.set(Calendar.SECOND, 0);
    long start = cCal.getTime().getTime();
    
    cCal.set(Calendar.MINUTE, 26);
    // System.out.println( "end: " + sdf.format( cCal.getTime() ) );
    long end = cCal.getTime().getTime();
    
    ds.trim(start, end);
    
    return ds;
  }

  @Test
  public void subtestLowFreqPSDs() {
    try {
      DataStore ds = setUpTest3();
      DataBlock cal = ds.getBlock(0);
      DataBlock out = ds.getBlock(1);
      int maxLen = Math.max( out.size(), cal.size() );
      int windowSize = 2;
      while (windowSize <= maxLen) {
        windowSize *= 2;
      }
      windowSize *= 2;
      FFTResult fft1 = 
          FFTResult.spectralCalcMultitaper(cal, out);
      FFTResult fft2 =
          FFTResult.spectralCalcMultitaper(out, out);
      Complex[] calSpec = fft1.getFFT();
      Complex[] outSpec = fft2.getFFT();
      double[] freqs = fft1.getFreqs();
      XYSeries calXYS = new XYSeries("Calibration spectrum amplitude");
      XYSeries outXYS = new XYSeries("Sensor out spectrum amplitude");
      XYSeries divXYS = new XYSeries("Power-spectral division");
      XYSeries subXYS = new XYSeries("Power-spectral dB subtraction");
      for (int i = 0; i < freqs.length; ++i) {
        if ( freqs[i] != 0 && freqs[i] <= .05 && freqs[i] >= .0001) {
          double calDB = 10 * Math.log10( calSpec[i].abs() );
          double outDB = 10 * Math.log10( outSpec[i].abs() );
          calXYS.add( 1/freqs[i], calDB );
          outXYS.add( 1/freqs[i], outDB );
          divXYS.add( 1/freqs[i], 
              10 * Math.log10( outSpec[i].divide(calSpec[i].abs()).abs() ) );
          subXYS.add( 1/freqs[i], outDB - calDB);
        }
      }
      
      XYSeriesCollection xysc = new XYSeriesCollection();
      xysc.addSeries(calXYS);
      xysc.addSeries(outXYS);
      JFreeChart chart = ChartFactory.createXYLineChart(
          "Low-Freq Cal Input verification",
          "Frequency",
          "Spectrum magnitude (log scale)",
          xysc,
          PlotOrientation.VERTICAL,
          true,
          false,
          false);
      
      xysc = new XYSeriesCollection(divXYS);
      xysc.addSeries(subXYS);
      JFreeChart chart2 = ChartFactory.createXYLineChart(
          "Low-Freq Cal Deconvolution verification",
          "Frequency",
          "Spectrum magnitude (log scale)",
          xysc,
          PlotOrientation.VERTICAL,
          true,
          false,
          false);
      LogarithmicAxis la = new LogarithmicAxis("Period (s)");
      
      chart.getXYPlot().setDomainAxis(la);
      chart2.getXYPlot().setDomainAxis(la);
      
      int width = 1280; int height = 960;
      
      PDDocument pdf = new PDDocument();
      ReportingUtils.chartsToPDFPage(width, height, pdf, chart, chart2);
      
      String currentDir = System.getProperty("user.dir");
      String testResultFolder = currentDir + "/testResultImages/";
      File dir = new File(testResultFolder);
      if ( !dir.exists() ) {
        dir.mkdir();
      }
      String testResult = 
          testResultFolder + "Low-frequency-spectral-data.pdf";
      pdf.save( new File(testResult) );
      pdf.close();
      System.out.println("Output result has been written");
      
    } catch (IOException e) {
      fail();
      e.printStackTrace();
    }
  }
  
  @Test
  public void subtestLowFreqPSDsMAJO() {
    try {
      String calname = "./data/random_cal_lowfrq/BC0.512.seed";
      String outname = "./data/random_cal_lowfrq/BHZ.512.seed";
      String calMplex = TimeSeriesUtils.getMplexNameList(calname).get(0);
      String outMplex = TimeSeriesUtils.getMplexNameList(outname).get(0);
      DataBlock cal = TimeSeriesUtils.getTimeSeries(calname, calMplex);
      DataBlock out = TimeSeriesUtils.getTimeSeries(outname, outMplex);
      long start = Math.max(cal.getStartTime(), out.getStartTime());
      long end = Math.min(cal.getEndTime(), out.getEndTime());
      cal.trim(start, end);
      out.trim(start, end);
      int maxLen = Math.max( out.size(), cal.size() );
      int windowSize = 2;
      while (windowSize <= maxLen) {
        windowSize *= 2;
      }
      windowSize *= 2;
      FFTResult fft1 = 
          FFTResult.spectralCalcMultitaper(cal, out);
      FFTResult fft2 =
          FFTResult.spectralCalcMultitaper(out, out);
      
      FFTResult fft3 =
          FFTResult.spectralCalc(cal, out);
      FFTResult fft4 =
          FFTResult.spectralCalc(out, out);
      Complex[] calSpec = fft1.getFFT();
      Complex[] outSpec = fft2.getFFT();
      double[] freqs1 = fft1.getFreqs();
      
      Complex[] calSpec2 = fft3.getFFT();
      Complex[] outSpec2 = fft4.getFFT();
      double[] freqs2 = fft3.getFreqs();
      
      XYSeries calXYS = new XYSeries("Cal vs out cross-amplitude MT");
      XYSeries outXYS = new XYSeries("Out vs out cross-amplitude MT");
      XYSeries calXYS2 = new XYSeries("Cal vs out cross-amplitude CT");
      XYSeries outXYS2 = new XYSeries("Out vs out cross-amplitude CT");
      XYSeries calPXYS = new XYSeries("Cal vs out cross-phase MT");
      XYSeries outPXYS = new XYSeries("Out vs out cross-phase MT");
      XYSeries calPXYS2 = new XYSeries("Cal vs out cross-phase CT");
      XYSeries outPXYS2 = new XYSeries("Out vs out cross-phase CT");
      for (int i = 0; i < freqs1.length; ++i) {
        if ( freqs1[i] != 0 && freqs1[i] <= .05 && freqs1[i] >= .0001) {
          double calDB = 10 * Math.log10( calSpec[i].abs() );
          double outDB = 10 * Math.log10( outSpec[i].abs() );
          double calP = NumericUtils.atanc(calSpec[i]);
          double outP = NumericUtils.atanc(outSpec[i]);
          calXYS.add( 1/freqs1[i], calDB );
          outXYS.add( 1/freqs1[i], outDB );
          calPXYS.add( 1/freqs1[i], calP );
          outPXYS.add( 1/freqs1[i], outP);
        }
      }
      
      for (int i = 0; i < freqs2.length; ++i) {
        if ( freqs2[i] != 0 && freqs2[i] <= .05 && freqs2[i] >= .0001) {
          double calDB = 10 * Math.log10( calSpec2[i].abs() );
          double outDB = 10 * Math.log10( outSpec2[i].abs() );
          double calP = NumericUtils.atanc(calSpec2[i]);
          double outP = NumericUtils.atanc(outSpec2[i]);
          calXYS2.add( 1/freqs2[i], calDB );
          outXYS2.add( 1/freqs2[i], outDB );
          calPXYS2.add( 1/freqs2[i], calP );
          outPXYS2.add( 1/freqs2[i], outP);
        }
      }
      
      XYSeriesCollection xysc = new XYSeriesCollection();
      xysc.addSeries(calXYS);
      xysc.addSeries(outXYS);
      xysc.addSeries(calXYS2);
      xysc.addSeries(outXYS2);
      JFreeChart chart = ChartFactory.createXYLineChart(
          "Low-Freq Cal Input verification",
          "Frequency",
          "Spectrum magnitude (log scale)",
          xysc,
          PlotOrientation.VERTICAL,
          true,
          false,
          false);
      
      xysc = new XYSeriesCollection();
      xysc.addSeries(calPXYS);
      xysc.addSeries(outPXYS);
      xysc.addSeries(calPXYS2);
      xysc.addSeries(outPXYS2);
      JFreeChart chart2 = ChartFactory.createXYLineChart(
          "Low-Freq Cal Deconvolution verification",
          "Frequency",
          "Spectrum magnitude (log scale)",
          xysc,
          PlotOrientation.VERTICAL,
          true,
          false,
          false);
      LogarithmicAxis la = new LogarithmicAxis("Period (s)");
      
      chart.getXYPlot().setDomainAxis(la);
      chart2.getXYPlot().setDomainAxis(la);
      
      int width = 1280; int height = 960;
      
      PDDocument pdf = new PDDocument();
      ReportingUtils.chartsToPDFPage(width, height, pdf, chart, chart2);
      
      String currentDir = System.getProperty("user.dir");
      String testResultFolder = currentDir + "/testResultImages/";
      File dir = new File(testResultFolder);
      if ( !dir.exists() ) {
        dir.mkdir();
      }
      String testResult = 
          testResultFolder + "Low-frequency-spectral-data_MAJO.pdf";
      pdf.save( new File(testResult) );
      pdf.close();
      System.out.println("Output result has been written");
      
    } catch (IOException e) {
      fail();
      e.printStackTrace();
    }
  }
  
  @Test
  public void testCalculationResult1() {
    
    String currentDir = System.getProperty("user.dir");
    // int testNumber = 3; // use to switch automated report data
    boolean lowFreq = false;
    
    try {
      
      DataStore ds = setUpTest1();
      InstrumentResponse ir = ds.getResponse(1);
      
      double nyq = ds.getBlock(0).getSampleRate() / 2.;
      System.out.println("NYQUIST RATE: " + nyq);
      
      RandomizedExperiment rCal = (RandomizedExperiment)
          ExperimentFactory.createExperiment(ExperimentEnum.RANDM);
      
      rCal.setLowFreq(lowFreq);
      
      assertTrue( rCal.hasEnoughData(ds) );
      rCal.runExperimentOnData(ds);
      
      double bestResid = rCal.getFitResidual();
      
      int width = 1280;
      int height = 960;
      
      List<XYSeriesCollection> xysc = rCal.getData();
      String[] yAxisTitles = new String[]{"Resp(f), dB", "Angle / TAU"};
      JFreeChart[] jfcl = new JFreeChart[yAxisTitles.length];
      
      String xAxisTitle = "Frequency (Hz)";
      NumberAxis xAxis = new LogarithmicAxis(xAxisTitle);
      Font bold = xAxis.getLabelFont().deriveFont(Font.BOLD);
      xAxis.setLabelFont(bold);
      
      StringBuilder sb = new StringBuilder();
      
      String[] resultString = RandomizedPanel.getInsetString(rCal);
      for (String resultPart : resultString) {
        sb.append( resultPart );
        sb.append('\n');
      }
      sb.append( RandomizedPanel.getTimeStampString(rCal) );
      sb.append('\n');
      sb.append("Input files:\n");
      sb.append( ds.getBlock(0).getName() );
      sb.append(" (calibration)\n");
      sb.append( ds.getBlock(1).getName() );
      sb.append(" (sensor output)\n");
      sb.append("Response file used:\n");
      sb.append( ds.getResponse(1).getName() );
      sb.append("\n \n");
      
      String page1 = sb.toString();
      
      String[] addtlPages = ( RandomizedPanel.getAdditionalReportPages(rCal) );
      // technically 'page 2' but really second part of first dataset report
      // and I'm too lazy to rename everything to reflect that
      String page1Part2 = addtlPages[0];
      
      sb = new StringBuilder();
      
      // expected best fit params, for debugging
      sb.append("BELOW RESULTS FOR EXPECTED BEST FIT (YELLOW CURVE)\n");
      double[] expectedParams = new double[]{-3.580104E+1, +7.122400E+1};
      ir = ir.buildResponseFromFitVector(expectedParams, lowFreq, 0);
      ir.setName("Best-fit params");
      ds.setResponse(1, ir);
      rCal.runExperimentOnData(ds);
      
      // residual from other code's best-fit parameters
      // compare to best-fit residual and assume difference is < 5%
      double expectedResid = rCal.getInitResidual();
      
      double pctDiff = 
          Math.abs( 100 * (bestResid - expectedResid) / bestResid );
      
      if (pctDiff > 15) {
        System.out.println(rCal.getFitPoles());
        System.out.println(rCal.getInitialPoles());
        System.out.println(bestResid + ", " + expectedResid);
      }
      assertTrue("PCT DIFF EXPECTED <15%, GOT " + pctDiff, pctDiff < 15);

      // add initial curve from expected fit params to report
      XYSeries expectedInitialCurve = rCal.getData().get(0).getSeries(0);
      xysc.get(0).addSeries(expectedInitialCurve);
      XYSeries expectedInitialAngle = rCal.getData().get(1).getSeries(0);
      xysc.get(1).addSeries(expectedInitialAngle);

      resultString = RandomizedPanel.getInsetString(rCal);
      for (String resultPart : resultString) {
        sb.append( resultPart );
        sb.append('\n');
      }
      
      for (int i = 0; i < jfcl.length; ++i) {
        
        jfcl[i] = ChartFactory.createXYLineChart(
            ExperimentEnum.RANDM.getName(),
            xAxisTitle,
            yAxisTitles[i],
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
      
      String page2 = sb.toString();
      
      PDDocument pdf = new PDDocument();
      ReportingUtils.chartsToPDFPage(width, height, pdf, jfcl);
      ReportingUtils.textListToPDFPages(pdf, page1, page1Part2, page2);
      
      String testResultFolder = currentDir + "/testResultImages/";
      File dir = new File(testResultFolder);
      if ( !dir.exists() ) {
        dir.mkdir();
      }
      
      String testResult = 
          testResultFolder + "Random-Calib-Test-1.pdf";
      pdf.save( new File(testResult) );
      pdf.close();
      System.out.println("Output result has been written");
      
    } catch (IOException e) {
      e.printStackTrace();
      fail();
    }
  }
  
  // @Test
  public void testCalculationResult2() {
    String currentDir = System.getProperty("user.dir");
    boolean lowFreq = false;
    
    try {
      
      DataStore ds = setUpTest2();
      // InstrumentResponse ir = ds.getResponse(1);
      
      RandomizedExperiment rCal = (RandomizedExperiment)
          ExperimentFactory.createExperiment(ExperimentEnum.RANDM);
      
      rCal.setLowFreq(lowFreq);
      
      assertTrue( rCal.hasEnoughData(ds) );
      rCal.runExperimentOnData(ds);
      
      int width = 1280;
      int height = 960;
      
      List<XYSeriesCollection> xysc = rCal.getData();
      String[] yAxisTitles = new String[]{"Resp(f), dB", "Angle / TAU"};
      JFreeChart[] jfcl = new JFreeChart[yAxisTitles.length];
      
      String xAxisTitle = "Frequency (Hz)";
      NumberAxis xAxis = new LogarithmicAxis(xAxisTitle);
      Font bold = xAxis.getLabelFont().deriveFont(Font.BOLD);
      xAxis.setLabelFont(bold);
      
      StringBuilder sb = new StringBuilder();
      String[] resultString = RandomizedPanel.getInsetString(rCal);
      for (String resultPart : resultString) {
        sb.append( resultPart );
        sb.append('\n');
      }
      sb.append('\n');
      sb.append( RandomizedPanel.getTimeStampString(rCal) );
      sb.append('\n');
      sb.append("Input files:\n");
      sb.append( ds.getBlock(0).getName() );
      sb.append(" (calibration)\n");
      sb.append( ds.getBlock(1).getName() );
      sb.append(" (sensor output)\n");
      sb.append("Response file used:\n");
      sb.append( ds.getResponse(1).getName() );
      sb.append("\n \n");
      
      String page1 = sb.toString();
      
      sb = new StringBuilder();
      
      for (int i = 0; i < jfcl.length; ++i) {
        
        jfcl[i] = ChartFactory.createXYLineChart(
            ExperimentEnum.RANDM.getName(),
            xAxisTitle,
            yAxisTitles[i],
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
      
      String page2 = sb.toString();
      
      PDDocument pdf = new PDDocument();
      ReportingUtils.chartsToPDFPage(width, height, pdf, jfcl);
      ReportingUtils.textListToPDFPages(pdf, page1, page2);
      
      String testResultFolder = currentDir + "/testResultImages/";
      File dir = new File(testResultFolder);
      if ( !dir.exists() ) {
        dir.mkdir();
      }
      
      String testResult = 
          testResultFolder + "Random-Calib-Test-2.pdf";
      pdf.save( new File(testResult) );
      pdf.close();
      System.out.println("Output result has been written");
      
    } catch (IOException e) {
      e.printStackTrace();
      fail();
    }
  }
  
  @Test
  public void testCalculationResult3() {
    String currentDir = System.getProperty("user.dir");
    boolean lowFreq = true;
    
    try {
      
      DataStore ds = setUpTest3();
      
      // InstrumentResponse ir = ds.getResponse(1);
      
      RandomizedExperiment rCal = (RandomizedExperiment)
          ExperimentFactory.createExperiment(ExperimentEnum.RANDM);
      
      rCal.setLowFreq(lowFreq);
      
      assertTrue( rCal.hasEnoughData(ds) );
      rCal.runExperimentOnData(ds);
      
      int width = 1280;
      int height = 960;
      
      List<XYSeriesCollection> xysc = rCal.getData();
      String[] yAxisTitles = new String[]{"Resp(f), dB", "Angle / TAU"};
      JFreeChart[] jfcl = new JFreeChart[yAxisTitles.length];
      
      String xAxisTitle = "Frequency (Hz)";
      NumberAxis xAxis = new LogarithmicAxis(xAxisTitle);
      Font bold = xAxis.getLabelFont().deriveFont(Font.BOLD);
      xAxis.setLabelFont(bold);
      
      StringBuilder sb = new StringBuilder();
      String[] resultString = RandomizedPanel.getInsetString(rCal);
      for (String resultPart : resultString) {
        sb.append( resultPart );
        sb.append('\n');
      }
      sb.append('\n');
      sb.append( RandomizedPanel.getTimeStampString(rCal) );
      sb.append('\n');
      sb.append("Input files:\n");
      sb.append( ds.getBlock(0).getName() );
      sb.append(" (calibration)\n");
      sb.append( ds.getBlock(1).getName() );
      sb.append(" (sensor output)\n");
      sb.append("Response file used:\n");
      sb.append( ds.getResponse(1).getName() );
      sb.append("\n \n");
      
      String page1 = sb.toString();
      
      sb = new StringBuilder();
      
      for (int i = 0; i < jfcl.length; ++i) {
        
        jfcl[i] = ChartFactory.createXYLineChart(
            ExperimentEnum.RANDM.getName(),
            xAxisTitle,
            yAxisTitles[i],
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
      
      String page2 = sb.toString();
      
      PDDocument pdf = new PDDocument();
      ReportingUtils.chartsToPDFPage(width, height, pdf, jfcl);
      ReportingUtils.textListToPDFPages(pdf, page1, page2);
      
      String testResultFolder = currentDir + "/testResultImages/";
      File dir = new File(testResultFolder);
      if ( !dir.exists() ) {
        dir.mkdir();
      }
      
      String testResult = 
          testResultFolder + "Random-Calib-Test-3.pdf";
      pdf.save( new File(testResult) );
      pdf.close();
      System.out.println("Output result has been written");
      
    } catch (IOException e) {
      e.printStackTrace();
      fail();
    }
  }
  
  @Test
  public void testCalculationResult4() {
    String currentDir = System.getProperty("user.dir");
    boolean lowFreq = false;
    
    try {
      
      DataStore ds = setUpTest4();
      
      // InstrumentResponse ir = ds.getResponse(1);
      
      RandomizedExperiment rCal = (RandomizedExperiment)
          ExperimentFactory.createExperiment(ExperimentEnum.RANDM);
      
      rCal.setLowFreq(lowFreq);
      
      assertTrue( rCal.hasEnoughData(ds) );
      rCal.runExperimentOnData(ds);
      
      int width = 1280;
      int height = 960;
      
      List<XYSeriesCollection> xysc = rCal.getData();
      String[] yAxisTitles = new String[]{"Resp(f), dB", "Angle / TAU"};
      JFreeChart[] jfcl = new JFreeChart[yAxisTitles.length];
      
      String xAxisTitle = "Frequency (Hz)";
      NumberAxis xAxis = new LogarithmicAxis(xAxisTitle);
      Font bold = xAxis.getLabelFont().deriveFont(Font.BOLD);
      xAxis.setLabelFont(bold);
      
      StringBuilder sb = new StringBuilder();
      String[] resultString = RandomizedPanel.getInsetString(rCal);
      for (String resultPart : resultString) {
        sb.append( resultPart );
        sb.append('\n');
      }
      sb.append('\n');
      sb.append( RandomizedPanel.getTimeStampString(rCal) );
      sb.append('\n');
      sb.append("Input files:\n");
      sb.append( ds.getBlock(0).getName() );
      sb.append(" (calibration)\n");
      sb.append( ds.getBlock(1).getName() );
      sb.append(" (sensor output)\n");
      sb.append("Response file used:\n");
      sb.append( ds.getResponse(1).getName() );
      sb.append("\n \n");
      
      String page1 = sb.toString();
      
      sb = new StringBuilder();
      
      for (int i = 0; i < jfcl.length; ++i) {
        
        jfcl[i] = ChartFactory.createXYLineChart(
            ExperimentEnum.RANDM.getName(),
            xAxisTitle,
            yAxisTitles[i],
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
      
      String page2 = sb.toString();
      
      PDDocument pdf = new PDDocument();
      ReportingUtils.chartsToPDFPage(width, height, pdf, jfcl);
      ReportingUtils.textListToPDFPages(pdf, page1, page2);
      
      String testResultFolder = currentDir + "/testResultImages/";
      File dir = new File(testResultFolder);
      if ( !dir.exists() ) {
        dir.mkdir();
      }
      
      String testResult = 
          testResultFolder + "Random-Calib-Test-4.pdf";
      pdf.save( new File(testResult) );
      pdf.close();
      System.out.println("Output result has been written");
      
    } catch (IOException e) {
      e.printStackTrace();
      fail();
    }
  }
}
