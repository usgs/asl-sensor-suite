package asl.sensor.test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.imageio.ImageIO;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.transform.DftNormalization;
import org.apache.commons.math3.transform.FastFourierTransformer;
import org.apache.commons.math3.transform.TransformType;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.LogarithmicAxis;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.junit.Test;

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
  public void demeaningTest() {
    
    // tests that demean does what it says it does and that
    // the results are applied in-place
    
    Number[] numbers = {1,2,3,4,5};
    
    List<Number> numList = Arrays.asList(numbers);
    List<Number> demeaned = new ArrayList<Number>(numList);
    
    FFTResult.demeanInPlace(demeaned);
    
    for (int i = 0; i < numList.size(); ++i) {
      assertEquals(demeaned.get(i), numList.get(i).doubleValue()-3);
    }
    
  }
  
  @Test
  public void detrendingCycleTest() {
    
    Number[] x = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 
        18, 19, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 
        3, 2, 1 };
    
    List<Number> toDetrend = Arrays.asList(x);
    
    Number[] answer = { -9d, -8d, -7d, -6d, -5d, -4d, -3d, -2d, -1d, 0d, 1d, 2d,
        3d, 4d, 5d, 6d, 7d, 8d, 9d, 10d, 9d, 8d, 7d, 6d, 5d, 4d, 3d, 2d, 1d, 0d,
        -1d, -2d, -3d, -4d, -5d, -6d, -7d, -8d, -9d };

    
    FFTResult.detrend(toDetrend);
    
    for (int i = 0; i < x.length; i++) {
      assertEquals(
          new Double(Math.round(x[i].doubleValue())), 
          new Double(answer[i].doubleValue()));
    }
    
  }
  
  @Test
  public void detrendingLinearTest() {
    
    Number[] x = { 1, 2, 3, 4, 5, 6, 7, 8, 9};
    
    List<Number> toDetrend = Arrays.asList(x);
    FFTResult.detrend(toDetrend);
    
    for (Number num : toDetrend) {
      assertEquals(num.doubleValue(), 0.0, 0.001);
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
      
      assertEquals(timeSeries[i], result[i], 10.);
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
  
}
