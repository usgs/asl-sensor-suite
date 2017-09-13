package asl.sensor.test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Calendar;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.TimeZone;

import javax.imageio.ImageIO;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.complex.ComplexFormat;
import org.apache.commons.math3.transform.DftNormalization;
import org.apache.commons.math3.transform.FastFourierTransformer;
import org.apache.commons.math3.transform.TransformType;
import org.apache.commons.math3.util.Pair;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.LogarithmicAxis;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.plot.SeriesRenderingOrder;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.junit.Test;

import asl.sensor.gui.InputPanel;
import asl.sensor.input.DataBlock;
import asl.sensor.utils.FFTResult;
import asl.sensor.utils.ReportingUtils;
import asl.sensor.utils.TimeSeriesUtils;

public class FFTResultTest {

  
  //@Test
  public void testComplexDivision() {
    // test is commented out because it is slow and unlikely to regress
    // because i and j would be the most intractable possible iteration vars
    for (int a = 1; a < 100000; ++a) {
      for (int b = a; b < 100000; ++b) {
        Complex numer = new Complex(a, b);
        Complex denom = new Complex(b, -a);
        Complex unit = numer.divide(numer); // expect 1
        Complex imag = numer.divide(denom); // expect i
        assertEquals(unit.getReal(), 1., 1E-15);
        assertEquals(unit.getImaginary(), 0., 1E-15);
        assertEquals(imag.getReal(), 0., 1E-15);
        assertEquals(imag.getImaginary(), 1., 1E-15);
      }
    }
  }
  
  @Test
  public void ohNoMorePrintFunctions() {
    
    // DecimalFormat df = new DecimalFormat(); // may need to tweak precision
    ComplexFormat cf = new ComplexFormat("j");
    
    // intended to be not-quite line-by-line replication of the PSD operations
    // it is said that unit tests are meant to be atomic. this ain't that
    String location = "testResultImages/psdOutputFiles/";
    File folder = new File(location);
    if ( !folder.exists() ) {
      folder.mkdir();
    }
    PrintWriter out;
    
    // get data first of all
    String inName = "data/noise_1/00_BH0.512.seed";
    DataBlock db;
    try {
      db = TimeSeriesUtils.getFirstTimeSeries(inName);
      Calendar startCal = db.getStartCalendar();
      startCal.set(Calendar.HOUR_OF_DAY, 0);
      startCal.set(Calendar.MINUTE, 59);
      startCal.set(Calendar.SECOND, 59);
      startCal.set(Calendar.MILLISECOND, 994);
      Calendar endCal = (Calendar) startCal.clone();
      endCal.set(Calendar.HOUR_OF_DAY, 7);
      endCal.set(Calendar.MINUTE, 0);
      endCal.set(Calendar.SECOND, 0);
      endCal.set(Calendar.MILLISECOND, 25);
      db.trim(startCal, endCal);
      
      double[] list1 = db.getData();
      String name = db.getName();
      String pref = location + name;
      out = new PrintWriter(pref + "-rawData.txt");
      out.write(Arrays.toString(list1));
      out.close();
      long interval = db.getInterval();
      // divide into windows of 1/4, moving up 1/16 of the data at a time
      
      final double TAPER_WIDTH = .1;
      int range = list1.length/4;
      int slider = range/4;
      
      // period is 1/sample rate in seconds
      // since the interval data is just that multiplied by a large number
      // let's divide it by that large number to get our period
      
      // shouldn't need to worry about a cast here
      double period = 1.0 / TimeSeriesUtils.ONE_HZ_INTERVAL;
      period *= interval;
      
      int padding = 2;
      while (padding < range) {
        padding *= 2;
      }
      
      int singleSide = padding / 2 + 1;
      double deltaFreq = 1. / (padding * period);
      
      Complex[] powSpectDens = new Complex[singleSide];
      double wss = 0;
      
      int segsProcessed = 0;
      int rangeStart = 0;
      int rangeEnd = range;
      
      for (int i = 0; i < powSpectDens.length; ++i) {
        powSpectDens[i] = Complex.ZERO;
      }
      
      while ( rangeEnd <= list1.length ) {
        
        String initValues, detrended, demeaned, tapered; 
        StringBuilder fftOutput = new StringBuilder();
        StringBuilder psdBins = new StringBuilder();
        
        Complex[] fftResult1 = new Complex[singleSide]; // first half of FFT reslt

        // give us a new list we can modify to get the data of
        double[] data1Range = 
            Arrays.copyOfRange(list1, rangeStart, rangeEnd);
         
        initValues = Arrays.toString(data1Range);
        
        // double arrays initialized with zeros, set as a power of two for FFT
        // (i.e., effectively pre-padded on initialization)
        double[] toFFT1 = new double[padding];
        
        // demean and detrend work in-place on the list
        TimeSeriesUtils.detrend(data1Range);
        detrended = Arrays.toString(data1Range);
        TimeSeriesUtils.demeanInPlace(data1Range);
        demeaned = Arrays.toString(data1Range);
        wss = FFTResult.cosineTaper(data1Range, TAPER_WIDTH);
        tapered = Arrays.toString(data1Range);
        // presumably we only need the last value of wss
        
        
        for (int i = 0; i < data1Range.length; ++i) { 
          toFFT1[i] = data1Range[i];
        }
        
        FastFourierTransformer fft = 
            new FastFourierTransformer(DftNormalization.STANDARD);

        Complex[] frqDomn1 = fft.transform(toFFT1, TransformType.FORWARD);
        fftOutput.append('[');
        for (int i = 0; i < frqDomn1.length; ++i) {
          fftOutput.append(cf.format(frqDomn1[i]));
          if (i + 1 < frqDomn1.length) {
            fftOutput.append(", ");
          }
        }
        fftOutput.append(']');
        // use arraycopy now (as it's fast) to get the first half of the fft
        System.arraycopy(frqDomn1, 0, fftResult1, 0, fftResult1.length);
        Complex[] psdBin = new Complex[singleSide];
        psdBins.append('[');
        for (int i = 0; i < singleSide; ++i) {
          Complex val1 = fftResult1[i];
          psdBin[i] = val1.multiply( val1.conjugate() );

          psdBins.append( psdBin[i].getReal() ); // imag part should be 0
          assertEquals(psdBin[i].getImaginary(), 0., 1E-15);
          if (i + 1 < singleSide) {
            psdBins.append(", ");
          }
          
          powSpectDens[i] = powSpectDens[i].add(psdBin[i]);
        }
        psdBins.append(']');
        ++segsProcessed;
        rangeStart  += slider;
        rangeEnd    += slider;
        
        out = 
            new PrintWriter(pref + "-psdSteps_" + segsProcessed + ".txt");
        StringBuilder sb = new StringBuilder();
        // initValues, detrended, demeaned, tapered, fftOutput, psdBins;
        sb.append("THIS IS THE DATA FOR WINDOW " + segsProcessed);
        sb.append(" WITH THE FOLLOWING FORMAT:\n");
        sb.append("1. RAW DATA GOING INTO PSD OF WINDOW SIZE\n");
        sb.append("2. DETRENDED DATA\n");
        sb.append("3. DATA WITH COSINE TAPER APPLIED\n");
        sb.append("4. OUTPUT OF FFT FOR MODIFIED DATA\n");
        sb.append("5. BINNED PSD DATA (SHOULD BE REAL-VALUED ONLY)\n\n");
        sb.append(initValues);
        sb.append("\n\n");
        sb.append(detrended);
        sb.append("\n\n");
        sb.append(demeaned);
        sb.append("\n\n");
        sb.append(tapered);
        sb.append("\n\n");
        sb.append(fftOutput);
        sb.append("\n\n");
        sb.append(psdBins);
        sb.append("\n\n");
        out.write(sb.toString());
        out.close();
      }
      
      // normalization time!
      
      double psdNormalization = 2.0 * period / padding;
      double windowCorrection = wss / (double) range;
      // it only uses the last value of wss, but that was how the original
      // code was
      
      psdNormalization /= windowCorrection;
      psdNormalization /= segsProcessed; // NOTE: divisor here should be 13
      
      double[] frequencies = new double[singleSide];
      
      for (int i = 0; i < singleSide; ++i) {
        powSpectDens[i] = powSpectDens[i].multiply(psdNormalization);
        frequencies[i] = i * deltaFreq;
      }
      
      out = new PrintWriter(pref + "-outputPSD.txt");
      out.write( Arrays.toString(powSpectDens) );
      out.close();
      // do smoothing over neighboring frequencies; values taken from 
      // asl.timeseries' PSD function
      int nSmooth = 11, nHalf = 5;
      Complex[] psdCFSmooth = new Complex[singleSide];
      
      int iw = 0;
      
      for (iw = 0; iw < nHalf; ++iw) {
        psdCFSmooth[iw] = powSpectDens[iw];
      }
      
      // iw should be icenter of nsmooth point window
      for (; iw < singleSide - nHalf; ++iw){
        int k1 = iw - nHalf;
        int k2 = iw + nHalf;
        
        Complex sumC = Complex.ZERO;
        for (int k = k1; k < k2; ++k) {
          sumC = sumC.add(powSpectDens[k]);
        }
        psdCFSmooth[iw] = sumC.divide((double) nSmooth);
      }
      
      // copy remaining into smoothed array
      for (; iw < singleSide; ++iw) {
        psdCFSmooth[iw] = powSpectDens[iw];
      }

      out = new PrintWriter(pref + "-smoothedPSD.txt");
      out.write( Arrays.toString(psdCFSmooth) );
      out.close();
      
    } catch (FileNotFoundException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
      fail();
    }
    
  }
  
  @Test
  public void cosineTaperTest() {
    double[] x = { 5, 5, 5, 5, 5 };
    double[] toTaper = x.clone();
    double[] tapered = { 0d, 4.5d, 5d, 4.5d, 0d };

    double power = FFTResult.cosineTaper(toTaper, 0.25);
    
    assertEquals(new Double(Math.round(power)), new Double(4));
    
    for (int i = 0; i < x.length; i++) {
      // precision to nearest tenth?
      assertEquals(toTaper[i], tapered[i], 0.1);
    }
  }
  
  @Test
  public void fftInversionTest() {
    double[] timeSeries = {10, 11, 12, 11, 10, 11, 12, 11, 10, 11, 12};
    
    int padSize = 2;
    while (padSize < timeSeries.length) {
      padSize *= 2;
    }
    
    double[] paddedTS = new double[padSize];
    for (int i = 0; i < timeSeries.length; ++i) {
      paddedTS[i] = timeSeries[i];
    }
    
    // System.out.println(paddedTS.length);
    
    FastFourierTransformer fft = 
        new FastFourierTransformer(DftNormalization.UNITARY);
    
    Complex[] frqDomn = fft.transform(paddedTS, TransformType.FORWARD);
    
    padSize = frqDomn.length/2 + 1;
    // System.out.println(padSize);
    
    Complex[] trim = new Complex[padSize];
    
    for (int i = 0; i < trim.length; ++i) {
      trim[i] = frqDomn[i];
    }
    
    padSize = (trim.length - 1) * 2;
    
    // System.out.println(padSize);
    
    Complex[] frqDomn2 = new Complex[padSize];
    
    for (int i = 0; i < padSize; ++i) {
      if (i < trim.length) {
        frqDomn2[i] = trim[i];
      } else {
        int idx = padSize - i;
        frqDomn2[i] = trim[idx].conjugate();
      }
      
      // System.out.println(frqDomn[i]+"|"+frqDomn2[i]);
      
    }
    
    Complex[] inverseFrqDomn = fft.transform(frqDomn2, TransformType.INVERSE);
    double[] result = new double[timeSeries.length];
   
    for (int i = 0; i < timeSeries.length; ++i) {
      result[i] = Math.round( inverseFrqDomn[i].getReal() );
      // System.out.println( result[i] + "," + inverseFrqDomn[i].getReal() );
      assertEquals(timeSeries[i], result[i], 0.1);
    }
    
  }
  
  @Test
  public void fftZerosTestMultitaper() {
    long interval = TimeSeriesUtils.ONE_HZ_INTERVAL;
    double[] data = new double[1000];
    // likely unnecessary loop, double arrays initialized at 0
    for (int i = 0; i < data.length; ++i) {
      data[i] = 0.;
    }
    FFTResult fftr = FFTResult.spectralCalcMultitaper(data, data, interval);
    Complex[] values = fftr.getFFT();
    for (Complex c : values) {
      assertTrue(c.equals(Complex.ZERO));
    }
  }
  
  @Test
  public void fftZerosTestWelch() {
    long interval = TimeSeriesUtils.ONE_HZ_INTERVAL;
    double[] data = new double[1000];
    // likely unnecessary loop, double arrays initialized at 0
    for (int i = 0; i < data.length; ++i) {
      data[i] = 0.;
    }
    FFTResult fftr = FFTResult.spectralCalc(data, data, interval);
    Complex[] values = fftr.getFFT();
    for (Complex c : values) {
      assertTrue(c.equals(Complex.ZERO));
    }
  }
  
  @Test
  public void lowPassFilterTest() {
    double[] timeSeries = new double[400];
    
    for (int i = 0; i < timeSeries.length; ++i) {
      if (i % 2 == 0) {
        timeSeries[i] = -10;
      } else
        timeSeries[i] = 10;
    }
    
    double sps = 40.;
    
    double[] lowPassed = FFTResult.bandFilter(timeSeries, sps, 0., 1.);
    
    for (int i = 1; i < (lowPassed.length - 1); ++i) {
      assertTrue( Math.abs( lowPassed[i] ) < 1. );
    }
    
  }
  
  //@Test
  public void multitaperSmootherData() {
    String name = "data/random_cal_lowfrq/BHZ.512.seed";
    try {
      DataBlock db = TimeSeriesUtils.getFirstTimeSeries(name);

      Calendar cal = db.getStartCalendar();
      Calendar cal2 = db.getStartCalendar();
      cal2.set(Calendar.HOUR_OF_DAY, 7);
      cal2.set(Calendar.MINUTE, 30);
      db.trim(cal, cal2);
      
      FastFourierTransformer fft = 
          new FastFourierTransformer(DftNormalization.STANDARD);
      
      int tCount = 12;
      double[] data = db.getData();
      double[][] tapers = FFTResult.getMultitaperSeries(data.length, tCount);
      
      // sum of smoothness of data
      double[] smoothness = new double[tCount];
      
      int padding = 1;
      while (padding < data.length) {
        padding *= 2;
      }
      
      for (int j = 0; j < tapers.length; ++j) {
        double[] taper = tapers[j];
        double[] taperedData = new double[padding];
        double smoothSum = 0;
        double taperSum = 0;
        for (int i = 0; i < data.length; ++i) {
          taperedData[i] = data[i] * taper[i];
          taperSum += Math.abs(taper[i]);
        }
        Complex[] fftData = fft.transform(taperedData, TransformType.FORWARD);
        for (int i = 1; i < fftData.length; ++i) {
          fftData[i] = fftData[i].divide(taperSum);
          Complex temp = fftData[i].subtract(fftData[i-1]);
          smoothSum += temp.abs();
        }
        
        smoothness[j] = smoothSum;
      }
      // smoother curves have lower value (less change between points)
      boolean improving = true;
      for (int i = 1; i < tCount; ++i) {
        System.out.println(smoothness[i] + ", prev " + smoothness[i-1]);
        improving &= smoothness[i] < smoothness[i-1];
      }
      assertTrue(improving);
      
    } catch (FileNotFoundException e) {
      e.printStackTrace();
      fail();
    }
  }

  
  @Test
  public void rangeCopyTest() {
  
    Number[] numbers = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
    
    int low = 5;
    int high = 9;
    
    List<Number> numList = Arrays.asList(numbers);
    List<Number> subseq = new ArrayList<Number>(numList.subList(low, high));
    
    for (int i = 0; i < subseq.size(); ++i) {
      int fullListIdx = i + low;
      assertEquals( numList.get(fullListIdx), subseq.get(i) );
    }
    
    for (int i = 0; i < subseq.size(); ++i) {
      Number temp = subseq.get(i);
      temp = 2000;
      assertNotEquals(subseq.get(i),temp); // can't try to change "in-place"
      subseq.set(i, 100);
    }
    
    for (int i = 0; i < subseq.size(); ++i) {
      int fullListIdx = i + low;
      assertNotEquals( numList.get(fullListIdx), subseq.get(i) );
    }
    
  }
  
  @Test
  public void showMultitaperPlot() {
    final int TAPERS = 12;
    StringBuilder sb = new StringBuilder();
    double[][] taper = FFTResult.getMultitaperSeries(512, TAPERS);
    XYSeriesCollection xysc = new XYSeriesCollection();
    for (int j = 0; j < taper.length; ++j) {
      XYSeries xys = new XYSeries("Taper " + j);
      double[] taperLine = taper[j];
      //System.out.println("TAPER LINE LEN: " + taperLine.length);
      //System.out.println("LAST VALUE: " + taperLine[taperLine.length - 1]);
      for (int i = 0; i < taperLine.length; ++i) {
        xys.add(i, taperLine[i]);
        sb.append(taperLine[i]);
        if ( i + 1 < taperLine.length) {
          sb.append(", ");
        }
      }
      sb.append("\n");
      xysc.addSeries(xys);
    }
    
    JFreeChart chart = 
        ChartFactory.createXYLineChart("MULTITAPER", "taper series index", 
            "taper value", xysc);
    
    BufferedImage bi = ReportingUtils.chartsToImage(1280, 960, chart);
    try {
      
      String testResultFolder = "testResultImages/";
      File folder = new File(testResultFolder);
      if ( !folder.exists() ) {
        System.out.println("Writing directory " + testResultFolder);
        folder.mkdirs();
      }
      
      File file = new File(testResultFolder + "multitaper plot.png");
      
      String testResult = 
          testResultFolder + "multitaper_ascii.csv";
      PrintWriter out = new PrintWriter(testResult);
      out.println( sb.toString() );
      out.close();
      
      ImageIO.write( bi, "png", file );
    } catch (IOException e) {
      e.printStackTrace();
      fail();
    }
  }

  
  @Test
  public void spectrumTest() {
    
    long interval = TimeSeriesUtils.ONE_HZ_INTERVAL;
    long timeStart = 0L;
    String name1 = "XX_FAKE_LH1_00";
    String name2 = "XX_FAKE_LH2_00";
    int len = 10000;
    double[] timeSeries = new double[len];
    double[] secondSeries = new double[len];
    
    for (int i = 0; i < len; ++i) {
      timeSeries[i] = Math.sin(i);
      secondSeries[i] = Math.sin(i) + Math.sin(2 * i);
    }
    
    DataBlock db = new DataBlock(timeSeries, interval, name1, timeStart);
    DataBlock db2 = new DataBlock(secondSeries, interval, name2, timeStart);
    FFTResult fft = FFTResult.spectralCalc(db, db2);
    
    XYSeriesCollection xysc = new XYSeriesCollection();
    String name = name1+"_"+name2;
    XYSeries xysr = new XYSeries(name + " spectrum (Real part)");
    XYSeries xysi = new XYSeries(name + " spectrum (Imag part)");
    for (int i = 0; i < fft.size(); ++i) {
      
      double freq = fft.getFreq(i);
      if ( freq <= 0. ) {
        continue;
      }
      
      xysr.add( freq, fft.getFFT(i).getReal() );
      xysi.add( freq, fft.getFFT(i).getImaginary() );
    }
    // xysc.addSeries(xysr);
    xysc.addSeries(xysi);
    
    JFreeChart jfc = ChartFactory.createXYLineChart(
        "SPECTRUM TEST CHART",
        "frequency",
        "value of spectrum",
        xysc);
    
    ValueAxis x = new LogarithmicAxis("frequency");
    jfc.getXYPlot().setDomainAxis(x);
    
    BufferedImage bi = ReportingUtils.chartsToImage(640, 480, jfc);
    String currentDir = System.getProperty("user.dir");
    String testResultFolder = currentDir + "/testResultImages/";
    File dir = new File(testResultFolder);
    if ( !dir.exists() ) {
      dir.mkdir();
    }
    
    String testResult = 
        testResultFolder + "spectrum.png";
    File file = new File(testResult);
    try {
      ImageIO.write(bi, "png", file);
    } catch (IOException e) {
      fail();
      e.printStackTrace();
    }
    
  }
  //@Test commented out because saves a lot of text data
  public void testAutomateRingler() {
    String name = "data/random_cal_lowfrq/BHZ.512.seed";
    try {
      DataBlock db = TimeSeriesUtils.getFirstTimeSeries(name);

      Calendar cal = db.getStartCalendar();
      Calendar cal2 = db.getStartCalendar();
      cal2.set(Calendar.HOUR_OF_DAY, 7);
      cal2.set(Calendar.MINUTE, 30);
      
      name = "data/random_cal_lowfrq/BC0.512.seed";
      
      db.trim(cal, cal2);
      
      String dir = "testResultImages/ringlerTaperResults/";
      File folder = new File(dir);
      if ( !folder.exists() ) {
        System.out.println("Writing directory " + dir);
        folder.mkdirs();
      }
      
      double[] data = db.getData();
      int padding = 1;
      while (padding < data.length) {
        padding *= 2;
      }
      
      double[] frequencies = new double[padding];
      double deltaFreq = db.getSampleRate() / padding;
      for (int i = 0; i < frequencies.length; ++i) {
        frequencies[i] = deltaFreq * i;
      }
      
      FastFourierTransformer fft = 
          new FastFourierTransformer(DftNormalization.STANDARD);
      
      StringBuilder sb = new StringBuilder();
      PrintWriter out;
      
      for (int i = 0; i < data.length; ++i) {
        sb.append(data[i]);
        if (i + 1 < data.length) {
          sb.append("\n");
        }
      }
      String dataOut = "raw_data_MAJO_BHZ.txt";
      out = new PrintWriter(dir + dataOut);
      out.println(sb.toString());
      out.close();
      
      double[] cloned = data.clone();
      TimeSeriesUtils.detrend(cloned);
      TimeSeriesUtils.demeanInPlace(cloned);
      int taperCt = 12;
      double[][] tapers = FFTResult.getMultitaperSeries(data.length, taperCt);
      
      int firstHalf = padding / 2 + 1;
      Complex[] fftFinal = new Complex[firstHalf];
      for (int i = 0; i < firstHalf; ++i) {
        fftFinal[i] = Complex.ZERO;
      }
      
      for (int j = 0; j < tapers.length; ++j) {
        double[] convert = new double[padding];
        double[] taper = tapers[j];
        double taperSum = 0;
        sb = new StringBuilder();
        for (int i = 0; i < data.length; ++i) {
          taperSum += Math.abs(taper[i]);
          convert[i] = cloned[i] * taper[i];
          
          sb.append(cloned[i]);
          sb.append(", ");
          sb.append(taper[i]);
          sb.append(", ");
          sb.append(convert[i]);
          if (i + 1 < data.length) {
            sb.append("\n");
          }
          
        }
        
        String file;
        
        file = "taper_" + j + ".csv";
        out = new PrintWriter(dir + file);
        out.println(sb.toString());
        out.close();
        
        
        Complex[] freqSpace = fft.transform(convert, TransformType.FORWARD);
        file = "fft_" + j + ".csv";
        sb = new StringBuilder();
        
        for (int i = 0; i < firstHalf; ++i) {
          fftFinal[i] = fftFinal[i].add(freqSpace[i].divide(taperSum));
          sb.append(frequencies[i]);
          sb.append(", ");
          sb.append(freqSpace[i].getReal());
          sb.append(", ");
          sb.append(freqSpace[i].getImaginary());
          if (i + 1 < freqSpace.length) {
            sb.append("\n");
          }
        }
        out = new PrintWriter(dir + file);
        out.println(sb.toString());
        out.close();
      }
      
      XYSeries xys = new XYSeries("PSD of MAJO LFQ - MULTI");
      XYSeries xysTest = new XYSeries("PSD of MAJO LFQ - WELCH");
      FFTResult fftr = FFTResult.spectralCalc(db, db);
      
      sb = new StringBuilder();
      Complex[] psdTest = fftr.getFFT();
      double[] freqList = fftr.getFreqs();
      for (int i = 1; freqList[i] < 0.1; ++i) {
        sb.append(freqList[i]);
        sb.append(", ");
        sb.append(psdTest[i].getReal());
        sb.append(", ");
        sb.append(psdTest[i].getImaginary());
        if (freqList[i+1] < 0.1) {
          sb.append("\n");
        }
        xysTest.add(freqList[i], 10 * Math.log10(psdTest[i].abs()));
      }
      out = new PrintWriter(dir + "BHZ-welch-fft.csv");
      out.println(sb.toString());
      out.close();
      
      sb = new StringBuilder();
      String outputFinal = "FFT_result.csv";
      for (int i = 0; i < firstHalf; ++i) {
        double fq = frequencies[i];
        fftFinal[i] = fftFinal[i].divide(taperCt);
        Complex cross = fftFinal[i].multiply( fftFinal[i].conjugate() );
        sb.append(fftFinal[i].getReal());
        sb.append(", ");
        sb.append(fftFinal[i].getImaginary());
        sb.append(", ");
        sb.append(cross.getReal());
        sb.append(", ");
        sb.append(cross.getImaginary());
        if (i + 1 < firstHalf) {
          sb.append('\n');
        }
        
        if (i > 0 && fq < 0.1) {
          xys.add( fq, 10 * Math.log10( cross.abs() ) );
        }
        
      }
      
      XYSeriesCollection xysc = new XYSeriesCollection();
      xysc.addSeries(xys);
      xysc.addSeries(xysTest);
      
      JFreeChart chart = 
          ChartFactory.createXYLineChart("psd", "x", "y", xysc);
      chart.getXYPlot().setDomainAxis(new LogarithmicAxis("frq"));
      
      BufferedImage bi = ReportingUtils.chartsToImage(640, 480, chart);
      String currentDir = System.getProperty("user.dir");
      String testResultFolder = currentDir + "/testResultImages/";
      String testResult = 
          testResultFolder + "PSD-MAJO.png";
      File file = new File(testResult);
      ImageIO.write(bi, "png", file);
      
      out = new PrintWriter(dir + outputFinal);
      out.write(sb.toString());
      out.close();
      
    } catch (FileNotFoundException e) {
      e.printStackTrace();
      fail();
    } catch (IOException e) {
      e.printStackTrace();
    }
  }
  
  //@Test commented out because saves a lot of text data
  public void testAutomateRingler2() {
    String name = "data/random_cal_lowfrq/BHZ.512.seed";
    try {
      DataBlock db = TimeSeriesUtils.getFirstTimeSeries(name);

      Calendar cal = db.getStartCalendar();
      Calendar cal2 = db.getStartCalendar();
      cal2.set(Calendar.HOUR_OF_DAY, 7);
      cal2.set(Calendar.MINUTE, 30);
      
      name = "data/random_cal_lowfrq/BC0.512.seed";
      db = TimeSeriesUtils.getFirstTimeSeries(name);
      
      db.trim(cal, cal2);
      
      String dir = "testResultImages/ringlerTaperResults2/";
      File folder = new File(dir);
      if ( !folder.exists() ) {
        System.out.println("Writing directory " + dir);
        folder.mkdirs();
      }
      
      double[] data = db.getData();
      int padding = 1;
      while (padding < data.length) {
        padding *= 2;
      }
      
      double[] frequencies = new double[padding];
      double deltaFreq = db.getSampleRate() / padding;
      for (int i = 0; i < frequencies.length; ++i) {
        frequencies[i] = deltaFreq * i;
      }
      
      FastFourierTransformer fft = 
          new FastFourierTransformer(DftNormalization.STANDARD);
      
      StringBuilder sb = new StringBuilder();
      PrintWriter out;
      
      for (int i = 0; i < data.length; ++i) {
        sb.append(data[i]);
        if (i + 1 < data.length) {
          sb.append("\n");
        }
      }
      String dataOut = "raw_data_MAJO_BC0.txt";
      out = new PrintWriter(dir + dataOut);
      out.println(sb.toString());
      out.close();
      
      double[] cloned = data.clone();
      TimeSeriesUtils.detrend(cloned);
      TimeSeriesUtils.demeanInPlace(cloned);
      int taperCt = 12;
      double[][] tapers = FFTResult.getMultitaperSeries(data.length, taperCt);
      
      int firstHalf = padding / 2 + 1;
      Complex[] fftFinal = new Complex[firstHalf];
      for (int i = 0; i < firstHalf; ++i) {
        fftFinal[i] = Complex.ZERO;
      }
      
      for (int j = 0; j < tapers.length; ++j) {
        double[] convert = new double[padding];
        double[] taper = tapers[j];
        double taperSum = 0;
        sb = new StringBuilder();
        for (int i = 0; i < data.length; ++i) {
          taperSum += Math.abs(taper[i]);
          convert[i] = cloned[i] * taper[i];
          
          sb.append(cloned[i]);
          sb.append(", ");
          sb.append(taper[i]);
          sb.append(", ");
          sb.append(convert[i]);
          if (i + 1 < data.length) {
            sb.append("\n");
          }
          
        }
        
        String file;
        
        file = "taper_" + j + ".csv";
        out = new PrintWriter(dir + file);
        out.println(sb.toString());
        out.close();
        
        
        Complex[] freqSpace = fft.transform(convert, TransformType.FORWARD);
        file = "fft_" + j + ".csv";
        sb = new StringBuilder();
        
        for (int i = 0; i < firstHalf; ++i) {
          fftFinal[i] = fftFinal[i].add(freqSpace[i].divide(taperSum));
          sb.append(frequencies[i]);
          sb.append(", ");
          sb.append(freqSpace[i].getReal());
          sb.append(", ");
          sb.append(freqSpace[i].getImaginary());
          if (i + 1 < freqSpace.length) {
            sb.append("\n");
          }
        }
        out = new PrintWriter(dir + file);
        out.println(sb.toString());
        out.close();
      }
      
      XYSeries xys = new XYSeries("PSD of MAJO LFQ - MULTI");
      XYSeries xysTest = new XYSeries("PSD of MAJO LFQ - WELCH");
      FFTResult fftr = FFTResult.spectralCalc(db, db);
      
      sb = new StringBuilder();
      Complex[] psdTest = fftr.getFFT();
      double[] freqList = fftr.getFreqs();
      for (int i = 1; freqList[i] < 0.1; ++i) {
        sb.append(freqList[i]);
        sb.append(", ");
        sb.append(psdTest[i].getReal());
        sb.append(", ");
        sb.append(psdTest[i].getImaginary());
        if (freqList[i+1] < 0.1) {
          sb.append("\n");
        }
        xysTest.add(freqList[i], 10 * Math.log10(psdTest[i].abs()));
      }
      out = new PrintWriter(dir + "BC0-welch-fft.csv");
      out.println(sb.toString());
      out.close();
      
      sb = new StringBuilder();
      String outputFinal = "FFT_result.csv";
      for (int i = 0; i < firstHalf; ++i) {
        double fq = frequencies[i];
        fftFinal[i] = fftFinal[i].divide(taperCt);
        Complex cross = fftFinal[i].multiply( fftFinal[i].conjugate() );
        sb.append(fftFinal[i].getReal());
        sb.append(", ");
        sb.append(fftFinal[i].getImaginary());
        sb.append(", ");
        sb.append(cross.getReal());
        sb.append(", ");
        sb.append(cross.getImaginary());
        if (i + 1 < firstHalf) {
          sb.append('\n');
        }
        
        if (i > 0 && fq < 0.1) {
          xys.add( fq, 10 * Math.log10( fftFinal[i].abs() ) );
          // xys.add( fq, 10 * Math.log10( cross.getReal() ) );
        }
        
      }
      
      XYSeriesCollection xysc = new XYSeriesCollection();
      xysc.addSeries(xys);
      xysc.addSeries(xysTest);
      
      JFreeChart chart = 
          ChartFactory.createXYLineChart("psd", "x", "y", xysc);
      chart.getXYPlot().setDomainAxis(new LogarithmicAxis("frq"));
      
      BufferedImage bi = ReportingUtils.chartsToImage(640, 480, chart);
      String currentDir = System.getProperty("user.dir");
      String testResultFolder = currentDir + "/testResultImages/";
      String testResult = 
          testResultFolder + "PSD-MAJO-CAL.png";
      File file = new File(testResult);
      ImageIO.write(bi, "png", file);
      
      out = new PrintWriter(dir + outputFinal);
      out.write(sb.toString());
      out.close();
      
    } catch (FileNotFoundException e) {
      e.printStackTrace();
      fail();
    } catch (IOException e) {
      e.printStackTrace();
    }
  }
  
  //@Test
  public void testDemeaning() {
    // temporarily commented out while I deal with this thing refactoring
    String dataFolderName = "data/random_cal/"; 
    String sensOutName = dataFolderName + "00_EHZ.512.seed";
    
    String metaName;
    try {
      metaName = TimeSeriesUtils.getMplexNameList(sensOutName).get(0);
      Pair<Long, Map<Long, double[]>> dataMap = 
          TimeSeriesUtils.getTimeSeriesMap(sensOutName, metaName);
      DataBlock sensor = TimeSeriesUtils.mapToTimeSeries(dataMap, metaName);
      
      XYSeriesCollection xysc = new XYSeriesCollection();
      XYSeries meaned = new XYSeries(metaName + "FFT, mean kept");
      XYSeries demeaned = new XYSeries(metaName + "FFT, mean removed");
      
      Map<Long, double[]> map = dataMap.getSecond();
      
      int padding = 2;
      while (padding < sensor.size()) {
        padding *= 2;
      }
      
      // System.out.println(map.size()+", "+sensor.size());
      
      double[] meanedTimeSeries = new double[padding];
      List<Long> time = new ArrayList<Long>( map.keySet() );
      Collections.sort(time);
      long interval = sensor.getInterval();
      
      int arrIdx = 0;
      for (int i = 0; i < time.size(); ++i) {
        long timeNow = time.get(i);
        // do a demean here so that we can add 0-values to empty points
        // without those values affecting the removal of a DC offset later
        for ( Number sample : map.get(timeNow) ) {
          meanedTimeSeries[arrIdx] = sample.doubleValue();
          ++arrIdx;
        }
        

        if ( (i + 1) < time.size() ) {
          long timeNext = time.get(i + 1);
          // is there a discrepancy, and is it big enough for us to care?
          if (timeNext - timeNow != interval) {
            
            if (timeNext - timeNow < interval / 2) {
              continue;
            }
            
            if (timeNext - timeNow < interval * 1.75) {
              continue; // assume just a time rounding error between records
            }
            
            // long gap = timeNext - timeNow;
            // System.out.println("FOUND GAP: " + timeNow + ", " + timeNext);
            // System.out.println("(Itvl: " + interval + "; gap: " + gap + ")");
            while (timeNext - timeNow > interval) {
              meanedTimeSeries[arrIdx] = 0.;
              ++arrIdx;
              timeNow += interval;
            }
          }
        }
      }
      
      assertEquals(arrIdx, sensor.size());
      
      double[] demeanedTimeSeries = sensor.getData();
      /*
      for (int i = 0; i < sensor.size(); ++i) {
        demeanedTimeSeries[i] = sensor.getData()[i];
      }
      */
      
      Complex[] meanedFFT = FFTResult.simpleFFT(meanedTimeSeries);
      Complex[] demeanedFFT = FFTResult.simpleFFT(demeanedTimeSeries);
      
      int numPoints = meanedFFT.length / 2 + 1;
      double nyquist = sensor.getSampleRate() / 2;
      double deltaFrq = nyquist / (numPoints - 1);
      
      double deltaFrqB = sensor.getSampleRate() / meanedFFT.length;
      
      // System.out.println("PEAK FREQ: " + deltaFrq * (numPoints - 1) );
      assertEquals(deltaFrq, deltaFrqB, 1E-7);
      
      for (int i = 1; i < numPoints; ++i) {
        double frq = i * deltaFrq;
        meaned.add(frq, 10 * Math.log10(meanedFFT[i].abs()));
        demeaned.add(frq, 10 * Math.log10(demeanedFFT[i].abs()));
      }
      
      xysc.addSeries(meaned);
      xysc.addSeries(demeaned);
      
      JFreeChart chart = ChartFactory.createXYLineChart(
          "Test demeaning operation",
          "frequency (hz)",
          "10 * log10 of FFT amplitude",
          xysc);
      chart.getXYPlot().setSeriesRenderingOrder(SeriesRenderingOrder.FORWARD);
      ValueAxis va = new LogarithmicAxis ("freq[hz]");
      chart.getXYPlot().setDomainAxis(va);
      
      
      String folderName = "testResultImages";
      File folder = new File(folderName);
      if ( !folder.exists() ) {
        System.out.println("Writing directory " + folderName);
        folder.mkdirs();
      }

      BufferedImage bi = ReportingUtils.chartsToImage(1280, 960, chart);
      File file = new File("testResultImages/demeaning-FFT-test.png");
      ImageIO.write( bi, "png", file );
      
    } catch (FileNotFoundException e) {
      fail();
      e.printStackTrace();
    } catch (IOException e) {
      fail();
      e.printStackTrace();
    }

  }
  
  @Test
  public void testMultitaper() {
    int size = 2000;
    List<Double> timeSeries = new ArrayList<Double>();
    for (int i = 0; i < size; ++i) {
      if (i % 2 == 0) {
        timeSeries.add(-500.);
      } else {
        timeSeries.add(500.);
      }
    }
    
    final int TAPERS = 12;
    double[][] taper = FFTResult.getMultitaperSeries(size, TAPERS);
    for (int j = 0; j < taper.length; ++j) {  
      double[] toFFT = new double[size];
      int l = toFFT.length-1; // last point
      double[] taperCurve = taper[j];
      double taperSum = 0.;
      //System.out.println(j + "-th taper curve first point: " + taperCurve[0]);
      //System.out.println(j + "-th taper curve last point: " + taperCurve[l]);
      for (int i = 0; i < timeSeries.size(); ++i) {
        taperSum += Math.abs(taperCurve[i]);
        double point = timeSeries.get(i).doubleValue();
        toFFT[i] = point * taperCurve[i];
      }
      //System.out.println(j + "-th tapered-data first point: " + toFFT[0]);
      //System.out.println(j + "-th tapered-data last point: " + toFFT[l]);
      
      assertEquals(0., toFFT[0], 1E-10);
      assertEquals(0., toFFT[l], 1E-10);
    }
  }
  
  @Test
  public void testTrimAndDemean() {
    String dataFolderName = "data/random_cal_4/"; 
    String sensOutName = dataFolderName + "00_EHZ.512.seed";
    
    String metaName;
    try {
      metaName = TimeSeriesUtils.getMplexNameList(sensOutName).get(0);
      Pair<Long, Map<Long, double[]>> dataMap = 
          TimeSeriesUtils.getTimeSeriesMap(sensOutName, metaName);
      DataBlock sensor = TimeSeriesUtils.mapToTimeSeries(dataMap, metaName);
      // long interval = dataMap.getFirst();
      Map<Long, double[]> timeSeriesMap = dataMap.getSecond();
      List<Long> times = new ArrayList<Long>( timeSeriesMap.keySet() );
      Collections.sort(times);
      long initTime = times.get(0);

      SimpleDateFormat sdf = InputPanel.SDF;
      sdf.setTimeZone( TimeZone.getTimeZone("UTC") );
      Calendar cCal = Calendar.getInstance( sdf.getTimeZone() );
      cCal.setTimeInMillis( initTime );
      cCal.set(Calendar.HOUR_OF_DAY, 20);
      cCal.set(Calendar.MINUTE, 16);
      cCal.set(Calendar.SECOND, 0);
      cCal.set(Calendar.MILLISECOND, 0);
      long startTime = cCal.getTimeInMillis();
      cCal.set(Calendar.MINUTE, 16 + 8); // match ringler's test case's window
      long endTime = cCal.getTimeInMillis();

      sensor.trim(startTime, endTime);
      System.out.println("Trimmed!");
      // System.out.println(startTime + "," + sensor.getStartTime());

      List<Double> timeSeriesList = new ArrayList<Double>();

      for (int i = 0; i < times.size(); ++i) {
        long timeNow = times.get(i);

        if (timeNow < startTime) {
          continue;
        } else if (timeNow > endTime) {
          break;
        }

        // do a demean here so that we can add 0-values to empty points
        // without those values affecting the removal of a DC offset later
        for (Number sample : timeSeriesMap.get(timeNow)) {
          timeSeriesList.add( sample.doubleValue() );
        }

      }
      
      // System.out.println(sensor.size());
      // System.out.println(timeSeriesList.size());
      // System.out.println(endTime+","+sensor.getEndTime());

      double[] timeSeriesArrayMean = new double[timeSeriesList.size()];
      double[] timeSeriesArrayDemean = sensor.getData();
      for (int i = 0; i < timeSeriesArrayMean.length; ++i) {
        timeSeriesArrayMean[i] = timeSeriesList.get(i);
        // timeSeriesArrayDemean[i] = sensor.getData()[i];
      }
      
      Complex[] withMeanFFT = FFTResult.simpleFFT(timeSeriesArrayMean);
      Complex[] demeanedFFT = FFTResult.simpleFFT(timeSeriesArrayDemean);
      
      int len = withMeanFFT.length / 2 + 1;
      double nyquist = sensor.getSampleRate() / 2;
      double deltaFrq = nyquist / (len - 1);
      
      // System.out.println("PEAK FREQ: " + (len - 1) * deltaFrq);

      XYSeries meanSeries = new XYSeries(metaName + " w/ mean");
      XYSeries demeanSeries = new XYSeries(metaName + " demeaned");
      
      for (int i = 1; i < len; ++i) {
        double xValue = i * deltaFrq;
        
        meanSeries.add(xValue, 10 * Math.log10( withMeanFFT[i].abs() ) );
        demeanSeries.add(xValue, 10 * Math.log10( demeanedFFT[i].abs() ) );
      }

      XYSeriesCollection xysc = new XYSeriesCollection();
      xysc.addSeries(meanSeries);
      xysc.addSeries(demeanSeries);
      JFreeChart chart = ChartFactory.createXYLineChart(
          "Test demeaning operation",
          "frequency (hz)",
          "10 * log10 of FFT amplitude",
          xysc);
      chart.getXYPlot().setSeriesRenderingOrder(SeriesRenderingOrder.FORWARD);
      ValueAxis va = new LogarithmicAxis ("freq[hz]");
      chart.getXYPlot().setDomainAxis(va);
      
      
      String folderName = "testResultImages";
      File folder = new File(folderName);
      if ( !folder.exists() ) {
        System.out.println("Writing directory " + folderName);
        folder.mkdirs();
      }
      
      BufferedImage bi = ReportingUtils.chartsToImage(1280, 960, chart);
      File file = new File("testResultImages/demeaning-FFT-trim-test.png");
      ImageIO.write( bi, "png", file );
      
    } catch (FileNotFoundException e) {
      fail();
      e.printStackTrace();
    } catch (IOException e) {
      fail();
      e.printStackTrace();
    }
      
  }
  
}
