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
  public void testAutomateRingler() {
    String name = "data/random_cal_lowfrq/BHZ.512.seed";
    try {
      DataBlock db = TimeSeriesUtils.getFirstTimeSeries(name);

      Calendar cal = db.getStartCalendar();
      Calendar cal2 = db.getStartCalendar();
      cal2.set(Calendar.HOUR_OF_DAY, 7);
      cal2.set(Calendar.MINUTE, 30);
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
      for (int i = 0; i < data.length; ++i) {
        sb.append(data[i]);
        if (i + 1 < data.length) {
          sb.append("\n");
        }
      }
      PrintWriter out;
      /*
      String dataOut = "raw_data_MAJO_BHZ.txt";
      out = new PrintWriter(dir + dataOut);
      out.println(sb.toString());
      out.close();
      */
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
        sb = new StringBuilder();
        for (int i = 0; i < data.length; ++i) {
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
        
        String file = "taper_" + j + ".csv";
        out = new PrintWriter(dir + file);
        out.println(sb.toString());
        out.close();
        
        Complex[] freqSpace = fft.transform(convert, TransformType.FORWARD);
        file = "fft_" + j + ".csv";
        sb = new StringBuilder();
        
        for (int i = 0; i < firstHalf; ++i) {
          fftFinal[i] = fftFinal[i].add(freqSpace[i]);
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
      
      sb = new StringBuilder();
      String outputFinal = "FFT_result.csv";
      for (int i = 0; i < firstHalf; ++i) {
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
      }
      
      out = new PrintWriter(dir + outputFinal);
      out.write(sb.toString());
      out.close();
      
    } catch (FileNotFoundException e) {
      e.printStackTrace();
      fail();
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
  
  //@Test TODO: uncomment this line once procedure understood
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
  
}
