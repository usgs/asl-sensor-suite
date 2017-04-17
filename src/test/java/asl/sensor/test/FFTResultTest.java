package asl.sensor.test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.transform.DftNormalization;
import org.apache.commons.math3.transform.FastFourierTransformer;
import org.apache.commons.math3.transform.TransformType;
import org.junit.Test;

import asl.sensor.utils.FFTResult;

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
