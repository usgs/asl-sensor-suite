package asl.sensor.test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
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
  public final void cosineTaperTest() throws Exception {
    Number[] x = { 5, 5, 5, 5, 5 };
    List<Number> toTaper = Arrays.asList(x);
    Double[] tapered = { 0d, 4.5d, 5d, 4.5d, 0d };

    double power = FFTResult.cosineTaper(toTaper, 0.25);
    
    assertEquals(new Double(Math.round(power)), new Double(4));
    
    for (int i = 0; i < x.length; i++) {
      // round to the first decimal (multiply by 10, round, divide by 10)
      Double result = 
          new Double(Math.round(toTaper.get(i).doubleValue()*10d)/10d);
      assertEquals(result, tapered[i]);
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
    List<Number> timeSeries = new ArrayList<Number>();
    List<Number> secondSeries = new ArrayList<Number>();
    
    for (int i = 0; i < 10000; ++i) {
      timeSeries.add( Math.sin(i) );
      secondSeries.add( Math.sin(i) + Math.sin(2 * i) );
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
      Pair<Long, Map<Long, Number>> dataMap = 
          TimeSeriesUtils.getTimeSeriesMap(sensOutName, metaName);
      DataBlock sensor = TimeSeriesUtils.mapToTimeSeries(dataMap, metaName);
      long interval = dataMap.getFirst();
      Map<Long, Number> timeSeriesMap = dataMap.getSecond();
      List<Long> times = new ArrayList<Long>( timeSeriesMap.keySet() );
      Collections.sort(times);
      long initTime = times.get(0);

      SimpleDateFormat sdf = InputPanel.SDF;
      sdf.setTimeZone( TimeZone.getTimeZone("UTC") );
      Calendar cCal = Calendar.getInstance( sdf.getTimeZone() );
      cCal.setTimeInMillis( initTime / 1000L );
      cCal.set(Calendar.HOUR_OF_DAY, 20);
      cCal.set(Calendar.MINUTE, 16);
      cCal.set(Calendar.SECOND, 0);
      cCal.set(Calendar.MILLISECOND, 0);
      long startTime = cCal.getTimeInMillis() * 1000L;
      cCal.set(Calendar.MINUTE, 16 + 8); // match ringler's test case's window
      long endTime = cCal.getTimeInMillis() * 1000L;

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
        timeSeriesList.add( timeSeriesMap.get(timeNow).doubleValue() );

      }
      
      // System.out.println(sensor.size());
      // System.out.println(timeSeriesList.size());
      // System.out.println(endTime+","+sensor.getEndTime());

      double[] timeSeriesArrayMean = new double[timeSeriesList.size()];
      double[] timeSeriesArrayDemean = new double[timeSeriesList.size()];
      for (int i = 0; i < timeSeriesArrayMean.length; ++i) {
        timeSeriesArrayMean[i] = timeSeriesList.get(i);
        timeSeriesArrayDemean[i] = sensor.getData().get(i).doubleValue();
      }
      
      Complex[] withMeanFFT = FFTResult.simpleFFT(timeSeriesArrayMean);
      Complex[] demeanedFFT = FFTResult.simpleFFT(timeSeriesArrayDemean);
      
      int len = withMeanFFT.length / 2 + 1;
      double nyquist = sensor.getSampleRate() / 2;
      double deltaFrq = nyquist / (len - 1);
      
      System.out.println("PEAK FREQ: " + (len - 1) * deltaFrq);

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
  
  // @Test
  public void testMultitaper() {
    final int TAPERS = 12;
    double[][] taper = FFTResult.getTaperSeries(2000, TAPERS);
    XYSeriesCollection xysc = new XYSeriesCollection();
    for (int j = 0; j < taper[0].length; ++j) {
      XYSeries xys = new XYSeries("Taper " + j);
      for (int i = 0; i < taper.length; ++i) {
        xys.add(i, taper[i][j]);
      }
      xysc.addSeries(xys);
    }
    
    JFreeChart chart = 
        ChartFactory.createXYLineChart("MULTITAPER", "taper series index", 
            "taper value", xysc);
    BufferedImage bi = ReportingUtils.chartsToImage(1280, 960, chart);
    File file = new File("testResultImages/multitaper plot.png");
    try {
      ImageIO.write( bi, "png", file );
    } catch (IOException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
      fail();
    }
  }
  
  @Test
  public void testDemeaning() {
    String dataFolderName = "data/random_cal/"; 
    String sensOutName = dataFolderName + "00_EHZ.512.seed";
    
    String metaName;
    try {
      metaName = TimeSeriesUtils.getMplexNameList(sensOutName).get(0);
      Pair<Long, Map<Long, Number>> dataMap = 
          TimeSeriesUtils.getTimeSeriesMap(sensOutName, metaName);
      DataBlock sensor = TimeSeriesUtils.mapToTimeSeries(dataMap, metaName);
      
      XYSeriesCollection xysc = new XYSeriesCollection();
      XYSeries meaned = new XYSeries(metaName + "FFT, mean kept");
      XYSeries demeaned = new XYSeries(metaName + "FFT, mean removed");
      
      Map<Long, Number> map = dataMap.getSecond();
      
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
        meanedTimeSeries[arrIdx] = map.get(timeNow).doubleValue();
        ++arrIdx;

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
      
      double[] demeanedTimeSeries = new double[padding];
      for (int i = 0; i < sensor.size(); ++i) {
        demeanedTimeSeries[i] = sensor.getData().get(i).doubleValue();
      }
      
      Complex[] meanedFFT = FFTResult.simpleFFT(meanedTimeSeries);
      Complex[] demeanedFFT = FFTResult.simpleFFT(demeanedTimeSeries);
      
      int numPoints = meanedFFT.length / 2 + 1;
      double nyquist = sensor.getSampleRate() / 2;
      double deltaFrq = nyquist / (numPoints - 1);
      
      double deltaFrqB = sensor.getSampleRate() / meanedFFT.length;
      
      System.out.println("PEAK FREQ: " + deltaFrq * (numPoints - 1) );
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
