package asl.sensor.utils;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.transform.DftNormalization;
import org.apache.commons.math3.transform.FastFourierTransformer;
import org.apache.commons.math3.transform.TransformType;
import org.jfree.data.xy.XYSeries;

import asl.sensor.input.DataBlock;
import asl.sensor.input.InstrumentResponse;

/**
 * Holds the data returned from a power spectral density calculation
 * (The PSD data (without response correction) and frequencies of the FFT)
 * Most methods that either calculate or involve FFT calculations exist here,
 * such as raw PSD calculation, inverse and forward trimmed FFTs,
 * and band-pass filtering.
 * @author akearns
 *
 */
public class FFTResult {
  
  /**
   * Specifies the width of the cosine taper function used in windowing
   */
  private static final double TAPER_WIDTH = 0.10;
  
  /**
   * Filter out data outside of the range between the low and high frequencies;
   * can be used for a low-pass filter if low frequency is set to 0
   * and high-pass if higher frequency is set to sample rate
   * @param toFilt series of data to do a band-pass filter on
   * @param sps sample rate of the current data (samples / sec)
   * @param low low corner frequency of band-pass filter
   * @param high high corner frequency of band-pass filter
   * @return timeseries with band-pass filter applied
   */
  public static double[] 
  bandFilter(double[] toFilt, double sps, double low, double high) {
    
    return bandFilterWithCuts(toFilt, sps, low, high, 0., sps);
    /*
    Complex[] fft = simpleFFT(toFilt);
    
    int trim = fft.length / 2 + 1;
    
    Complex[] toInvert = new Complex[trim];
    
    double freqDelta = sps / trim;
    
    for (int i = 0; i < trim; ++i) {
      double x = i * freqDelta;
      double scale = 1;
      if (x < low) {
        scale = x / low;
      } else if (x > high) {
        scale = 1 - ( x / (sps - high) );
      }
      
      toInvert[i] = fft[i].multiply(scale);
    }
    
    return singleSidedInverseFFT(toInvert, toFilt.length);
    */
    
  }
  
  /**
   * Wrapper to do band filter on a list of data rather than an array.
   * For more details see other definition of bandFilter
   * @param toFilt timeseries data to be filtered
   * @param sps samples per second of input data
   * @param low low corner frequency for trim
   * @param high higher corner frequency for trim
   * @return timeseries data (list) that has gone through band-pass filter
   */
  public static List<Number> 
  bandFilter(List<Number> toFilt, double sps, double low, double high) {
    
    
    double[] toFFT = new double[toFilt.size()];
    
    for (int i = 0; i < toFFT.length; ++i) {
      toFFT[i] = toFilt.get(i).doubleValue();
    }
    
    toFFT = bandFilter(toFFT, sps, low, high);
    
    List<Number> out = new ArrayList<Number>();
    for (double value : toFFT) {
      out.add(value);
    }
    
    return out;
    
  }
  
  /**
   * Band-pass filter that creates hard-stop for values outside of range but
   * produces linear dropoff for points between the corner frequencies and
   * hard stop limits.
   * @param toFilt series of data to do a band-pass filter on
   * @param sps sample rate of the current data (samples / sec)
   * @param low low corner frequency of band-pass filter
   * @param high high corner frequency of band-pass filter
   * @param lowStop low frequency value beyond which to attenuate all signal
   * @param highStop high frequency value beyond which to attenuate all signal
   * @return timeseries with band-pass filter applied
   */
  public static double[] 
  bandFilterWithCuts(double[] toFilt, double sps, double low, double high, 
                     double lowStop, double highStop) {
    
    System.out.println("FILTERING OPERATION OCCURRING");
    
    Complex[] fft = simpleFFT(toFilt);
    
    int trim = fft.length/2 + 1;
    
    Complex[] toInvert = new Complex[trim];
    
    double freqDelta = sps / trim;
    
    for (int i = 0; i < trim; ++i) {
      double x = i * freqDelta;
      double scale = 1;
      if (x < lowStop || x > highStop) {
        scale = 0;
      } else if (x < low) {
        scale = (x - lowStop) / (low - lowStop);
      } else if (x > high) {
        scale = (x - highStop) / (high - highStop);
      }
      
      toInvert[i] = fft[i].multiply(scale);
    }
    
    return singleSidedInverseFFT(toInvert, toFilt.length);
  }
  
  /**
   * Calculates and performs an in-place cosine taper on an incoming data set.
   * Used for windowing for performing FFT.
   * @param dataSet The dataset to have the taper applied to.
   * @param taperW Width of taper to be used
   * @return Value corresponding to power loss from application of taper.
   */
  public static double cosineTaper(List<Number> dataSet, double taperW) {
    
    double ramp = taperW * dataSet.size();
    double taper;
    double wss = 0.0; // represents power loss
    
    for (int i = 0; i < ramp; i++) {
      taper = 0.5 * (1.0 - Math.cos( (double) i * Math.PI / ramp) );
      dataSet.set(i, dataSet.get(i).doubleValue() * taper);
      int idx = dataSet.size()-i-1;
      dataSet.set(idx, dataSet.get(idx).doubleValue() * taper );
      wss += 2.0 * taper * taper;
    }
    
    wss += ( dataSet.size() - (2 * ramp) );
    
    return wss;
  }
  
  public static double[][] getCosTaperVector(int len, double taperW) {
    
    double[] taperVec = new double[len];
    for (int i = 0; i < taperVec.length; ++i) {
      taperVec[i] = 1.;
    }
    
    double ramp = taperW * len;
    double taper;
    
    for (int i = 0; i < ramp; i++) {
      taper = 0.5 * (1.0 - Math.cos( (double) i * Math.PI / ramp) );
      taperVec[i] *= taper;
      int idx = taperVec.length-i-1;
      taperVec[idx] *= taper;
    }
    
    return new double[][]{taperVec};
    
  }
  
  public static double[][] getTaperSeries(int winLen, int numTapers) {
    double[][] taperMat = new double[winLen][numTapers];
    
    double denom = winLen + 1;
    double scale = Math.sqrt(2 / denom);
    
    for (int i = 0; i < winLen; ++i) {
      for (int j = 0; j < numTapers; ++j) {
        taperMat[i][j] = scale * Math.sin(Math.PI * (i + 1) * (j + 1) / denom);
      }
    }
    
    return taperMat;
  }
  
  /**
   * Root funtion for calculating crosspower. Gets spectral calculation of data
   * from inputted data series by calling the spectralCalc function, and then
   * applies the provided responses to that result. This is the Power Spectral
   * Density of the inputted data if both sets are the same.
   * @param data1 First data series
   * @param data2 Second data series
   * @param ir1 Response of instrument producing first series
   * @param ir2 Response of instrument producing second series
   * @return Data structure containing the crosspower result of the two data
   * sets as a complex array and the frequencies matched to them in a double
   * array. 
   */
  public static FFTResult crossPower(DataBlock data1, DataBlock data2,
      InstrumentResponse ir1, InstrumentResponse ir2) {
    
    FFTResult selfPSD = spectralCalc(data1, data2);
    Complex[] results = selfPSD.getFFT();
    double[] freqs = selfPSD.getFreqs();
    Complex[] out = new Complex[freqs.length];
    Complex[] freqRespd1 = ir1.applyResponseToInput(freqs);
    Complex[] freqRespd2 = ir2.applyResponseToInput(freqs);
    
    for (int j = 0; j < freqs.length; ++j) {
      Complex respMagnitude = 
          freqRespd1[j].multiply( freqRespd2[j].conjugate() );
      
      if (respMagnitude.abs() == 0) {
        respMagnitude = new Complex(Double.MIN_VALUE, 0);
      }
      
      out[j] = results[j].divide(respMagnitude);
    }
    
    return new FFTResult(out, freqs);
    
  }
  
  /**
   * Collects the data points in the Peterson new high noise model 
   * into a plottable format.
   * Assumes that there is a text file in the .resources folder that contains 
   * the NHNM data points for given input frequencies.
   * @param freqSpace True if the data's x-axis should be units of Hz
   * (otherwise it is units of seconds, the interval between samples)
   * @return Plottable data series representing the NHNM
   */
  public static XYSeries getHighNoiseModel(boolean freqSpace) {
    XYSeries xys = new XYSeries("NHNM");
    try {
      ClassLoader cl = FFTResult.class.getClassLoader();
      InputStream is = cl.getResourceAsStream("NHNM.txt");
      
      BufferedReader fr = new BufferedReader( new InputStreamReader(is) );
      String str = fr.readLine();
      while (str != null) {
        String[] values = str.split("\\s+");
        double x = Double.parseDouble(values[0]); // period, in seconds
        if (x > 1.0E3) {
          break;
        }
        double y = Double.parseDouble(values[1]);
        if (freqSpace) {
          xys.add(1/x, y);
        } else {
          xys.add(x, y);
        }
        
        str = fr.readLine();
      }
      fr.close();
    } catch (FileNotFoundException e) {
      e.printStackTrace();
    } catch (IOException e) {
      e.printStackTrace();
    }
    return xys;
  }
  
  
  /**
   * Collects the data points in the Peterson new low noise model 
   * into a plottable format.
   * Assumes that there is a text file in the .resources folder that contains 
   * the NLNM data points for given input frequencies. ("NLNM.txt")
   * @param freqSpace True if the data's x-axis should be units of Hz
   * (otherwise it is units of seconds, the interval between samples)
   * @return Plottable data series representing the NLNM
   */
  public static XYSeries getLowNoiseModel(boolean freqSpace) {
    XYSeries xys = new XYSeries("NLNM");
    try {
      
      ClassLoader cl = FFTResult.class.getClassLoader();
      
      InputStream is = cl.getResourceAsStream("NLNM.txt");
      
      BufferedReader fr = new BufferedReader( new InputStreamReader(is) );
      String str = fr.readLine();
      while (str != null) {
        String[] values = str.split("\\s+");
        double x = Double.parseDouble(values[0]); // period, in seconds
        if (x > 1.0E3) {
          break;
        }
        double y = Double.parseDouble(values[3]);
        if (freqSpace) {
          xys.add(1/x, y);
        } else {
          xys.add(x, y);
        }
        
        str = fr.readLine();
      }
      fr.close();
    } catch (FileNotFoundException e) {
      e.printStackTrace();
    } catch (IOException e) {
      e.printStackTrace();
    }
    return xys;
  }
  
  /**
   * Function for padding and returning the result of a forward FFT.
   * This does not trim the negative frequencies of the result; it returns
   * the full FFT result as an array of Complex numbers
   * @param dataIn Array of doubles representing timeseries data
   * @return Complex array representing forward FFT values, including
   * symmetric component (second half of the function)
   */
  public static Complex[] simpleFFT(double[] dataIn) {
    
    int padding = 2;
    while ( padding < dataIn.length ) {
      padding *= 2;
    }
    
    double[] toFFT = new double[padding];
    
    for (int i = 0; i < dataIn.length; ++i) {
      toFFT[i] = dataIn[i];
    }
    
    FastFourierTransformer fft = 
        new FastFourierTransformer(DftNormalization.STANDARD);
    
    Complex[] frqDomn = fft.transform(toFFT, TransformType.FORWARD);
    
    return frqDomn;
  }
  
  /**
   * Calculates the FFT of the timeseries data in a DataBlock
   * and returns the positive frequencies resulting from the FFT calculation
   * @param db DataBlock to get the timeseries data from
   * @param mustFlip True if signal from sensor is inverted (for step cal)
   * @return Complex array of FFT values and double array of corresponding 
   * frequencies 
   */
  public static FFTResult singleSidedFFT(DataBlock db, boolean mustFlip) {
    
    double[] data = new double[db.size()];
    
    for (int i = 0; i < db.size(); ++i) {
      data[i] = db.getData().get(i).doubleValue();
      
      if (mustFlip) {
        data[i] *= -1;
      }
    }
    
    data = TimeSeriesUtils.demean(data);
    
    // data = TimeSeriesUtils.normalize(data);
    
    Complex[] frqDomn = simpleFFT(data);
    
    int padding = frqDomn.length;
    int singleSide = padding/2 + 1;
    
    double nyquist = db.getSampleRate() / 2;
    double deltaFrq = nyquist / (singleSide - 1);
    
    Complex[] fftOut = new Complex[singleSide];
    double[] frequencies = new double[singleSide];
    
    for (int i = 0; i < singleSide; ++i) {
      fftOut[i] = frqDomn[i];
      frequencies[i] = i * deltaFrq;
    }
    
    // System.out.println(frequencies[singleSide - 1]);
    
    return new FFTResult(fftOut, frequencies);
    
  }
  
  /**
   * Calculates the FFT of the timeseries data in a DataBlock
   * and returns the positive frequencies resulting from the FFT calculation
   * @param db DataBlock to get the timeseries data from
   * @param mustFlip True if signal from sensor is inverted (for step cal)
   * @return Complex array of FFT values and double array of corresponding 
   * frequencies 
   */
  public static FFTResult 
  singleSidedFilteredFFT(DataBlock db, boolean mustFlip) {
    
    double[] data = new double[db.size()];
    
    for (int i = 0; i < db.size(); ++i) {
      data[i] = db.getData().get(i).doubleValue();
      
      if (mustFlip) {
        data[i] *= -1;
      }
    }
    
    long interval = db.getInterval();

    double sps = TimeSeriesUtils.ONE_HZ_INTERVAL / interval;
    
    data = bandFilter(data, sps, 0.0, 0.1);
    
    data = TimeSeriesUtils.demean(data);
    
    // data = TimeSeriesUtils.normalize(data);
    
    Complex[] frqDomn = simpleFFT(data);
    
    int padding = frqDomn.length;
    int singleSide = padding/2 + 1;
    
    double nyquist = db.getSampleRate() / 2;
    double deltaFrq = nyquist / (singleSide - 1);
    
    Complex[] fftOut = new Complex[singleSide];
    double[] frequencies = new double[singleSide];
    
    for (int i = 0; i < singleSide; ++i) {
      fftOut[i] = frqDomn[i];
      frequencies[i] = i * deltaFrq;
    }
    
    return new FFTResult(fftOut, frequencies);
    
  }
  
  /**
   * Do the inverse FFT on the result of a single-sided FFT operation.
   * The negative frequencies are reconstructed as the complex conjugates of
   * the positive corresponding frequencies
   * @param freqDomn Complex array (i.e., the result of a previous FFT calc)
   * @param trim How long the original input data was
   * @return A list of doubles representing the original timeseries of the FFT
   */
  public static double[] singleSidedInverseFFT(Complex[] freqDomn, int trim) {
    FastFourierTransformer fft = 
        new FastFourierTransformer(DftNormalization.STANDARD);
     
    int padding = (freqDomn.length - 1) * 2;
    
    Complex[] padded = new Complex[padding];
    for (int i = 0; i < freqDomn.length; ++i) {
      padded[i] = freqDomn[i];
    }
    for (int i = 1; i < padding/2; ++i) {
      // System.out.println(freqDomn.length+","+i);
      padded[padded.length - i] = padded[i].conjugate();
    }
    
    Complex[] timeSeriesCpx = 
        fft.transform(padded, TransformType.INVERSE);
    
    double[] timeSeries = new double[trim];
    for (int i = 0; i < trim; ++i) {
      timeSeries[i] = timeSeriesCpx[i].getReal();
    }
    
    return timeSeries;
  }
  
  /**
   * Helper function to calculate power spectral density / crosspower.
   * Takes in two time series data and produces the windowed FFT over each.
   * The first is multiplied by the complex conjugate of the second.
   * If the two series are the same, this is the PSD of that series. If they
   * are different, this result is the crosspower.
   * The result is smoothed but does not have the frequency response applied,
   * and so does not give a full result -- this is merely a helper function
   * for the crossPower function.
   * @param data1 DataBlock with relevant time series data
   * @param data2 DataBlock with relevant time series data
   * @return A structure with two arrays: an array of Complex numbers 
   * representing the PSD result, and an array of doubles representing the
   * frequencies of the PSD.
   */
  public static FFTResult spectralCalc(DataBlock data1, DataBlock data2) {

    List<Number> list1 = data1.getData();
    
    // divide into windows of 1/4, moving up 1/16 of the data at a time
    int range = list1.size()/4;
    int slider = range/4;
    
    FFTResult data = spectralCalc(data1, data2, range, slider, TaperType.COS);
    Complex[] powSpectDens = data.getFFT();
    double[] frequencies = data.getFreqs();
    
    // do smoothing over neighboring frequencies; values taken from 
    // asl.timeseries' PSD function
    int singleSide = frequencies.length;
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
    
    return new FFTResult(psdCFSmooth, frequencies);
    
  }
  
  public static FFTResult spectralCalc(DataBlock data1, DataBlock data2, 
      int range, int slider, TaperType taper) {
    
    // this is ugly logic here, but this saves us issues with looping
    // and calculating the same data twice
    boolean sameData = data1.getName().equals( data2.getName() );
    
    List<Number> list1 = data1.getData();
    List<Number> list2 = list1;
    if (!sameData) {
      list2 = data2.getData();
    }
    
    // period is 1/sample rate in seconds
    // since the interval data is just that multiplied by a large number
    // let's divide it by that large number to get our period
    
    // shouldn't need to worry about a cast here
    
    int padding = 2;
    while (padding < range) {
      padding *= 2;
    }
    
    final int TAPER_COUNT = 12;
    double period = 1.0 / TimeSeriesUtils.ONE_HZ_INTERVAL;
    period *= data1.getInterval();
    
    int singleSide = padding / 2 + 1;
    double deltaFreq = 1. / (padding * period);
    
    Complex[] powSpectDens = new Complex[singleSide];
    double wss = 1.;
    
    int segsProcessed = 0;
    int rangeStart = 0;
    int rangeEnd = range;
    
    for (int i = 0; i < powSpectDens.length; ++i) {
      powSpectDens[i] = Complex.ZERO;
    }
    
    do {
      
      Complex[] fftResult1 = new Complex[singleSide]; // first half of FFT reslt
      Complex[] fftResult2 = fftResult1;
      
      int upperBound = Math.min( rangeEnd, list1.size() );
      
      if (!sameData) {
        fftResult2 = new Complex[singleSide];
      }
      
      // give us a new list we can modify to get the data of
      List<Number> data1Range = 
          new ArrayList<Number>(
              list1.subList(rangeStart, upperBound) );
      
      List<Number> data2Range = data1Range;
      if (!sameData) {
        data2Range = 
            new ArrayList<Number>(
                list2.subList(rangeStart, upperBound) );
      }
       
      // double arrays initialized with zeros, set as a power of two for FFT
      // (i.e., effectively pre-padded on initialization)
      double[] toFFT1 = new double[padding];
      double[] toFFT2 = toFFT1;
      
      // demean and detrend work in-place on the list
      TimeSeriesUtils.detrend(data1Range);
      TimeSeriesUtils.demeanInPlace(data1Range);
      switch (taper) {
      case COS:
        wss = cosineTaper(data1Range, TAPER_WIDTH);
        break;
      case MULT:
      default:
        wss = 0;
        double[][] taperMat = getTaperSeries(data1Range.size(), TAPER_COUNT);
        for (int i = 0; i < data1Range.size(); ++i) {
          double point = data1Range.get(i).doubleValue();
          double sum = 0.;
          for (int j = 0; j < taperMat[0].length; ++j) {
            sum += point * taperMat[i][j];
            wss += Math.pow(taperMat[i][j], 2);
          }
          // need to re-assign point here; primitives not assigned by reference
          // point = data1Range.get(i).doubleValue();
          wss /= TAPER_COUNT;
          data1Range.set(i, sum / TAPER_COUNT);
        }
        break;
      }
      // presumably we only need the last value of wss
      
      if (!sameData) {
        TimeSeriesUtils.demeanInPlace(data2Range);
        TimeSeriesUtils.detrend(data2Range);
        switch (taper) {
        case COS:
          wss = cosineTaper(data2Range, TAPER_WIDTH);
          break;
        case MULT:
        default:
          // wss is already calculated above
          double[][] taperMat = getTaperSeries(data2Range.size(), TAPER_COUNT);
          for (int i = 0; i < data2Range.size(); ++i) {
            double point = data2Range.get(i).doubleValue();
            double sum = 0.;
            for (int j = 0; j < taperMat[0].length; ++j) {
              sum += point * taperMat[i][j];
            }
            // need to re-assign point here; not evaluated by reference
            // point = data2Range.get(i).doubleValue();
            data2Range.set(i, sum / TAPER_COUNT);
          }
          break;
        }
        toFFT2 = new double[padding];
      }

      for (int i = 0; i < data1Range.size(); ++i) {
        // no point in using arraycopy -- must make sure each Number's a double
        toFFT1[i] = data1Range.get(i).doubleValue();
        if (!sameData) {
          toFFT2[i] = data2Range.get(i).doubleValue();
        }
      }
      
      FastFourierTransformer fft = 
          new FastFourierTransformer(DftNormalization.STANDARD);

      Complex[] frqDomn1 = fft.transform(toFFT1, TransformType.FORWARD);
      // use arraycopy now (as it's fast) to get the first half of the fft
      System.arraycopy(frqDomn1, 0, fftResult1, 0, fftResult1.length);
      
      Complex[] frqDomn2 = frqDomn1;
      if (!sameData) {
        frqDomn2 = fft.transform(toFFT2, TransformType.FORWARD);
        System.arraycopy(frqDomn2, 0, fftResult2, 0, fftResult2.length);
      }
      
      for (int i = 0; i < singleSide; ++i) {
        
        Complex val1 = fftResult1[i];
        Complex val2 = val1;
        if (!sameData) {
          val2 = fftResult2[i];
        }
        
        powSpectDens[i] = 
            powSpectDens[i].add( 
                val1.multiply( 
                    val2.conjugate() ) );
      }
      
      ++segsProcessed;
      rangeStart  += slider;
      rangeEnd    += slider;
      
    } while ( rangeEnd <= data1.size() );
    
    // normalization time!
    double psdNormalization = 2.0 * period / padding;
    double windowCorrection = wss / (double) range;
    // it only uses the last value of wss, but that was how the original
    // code was
    
    psdNormalization /= windowCorrection;
    // System.out.println(segsProcessed);
    psdNormalization /= segsProcessed; // NOTE: divisor here should be 13
    
    double[] frequencies = new double[singleSide];
    
    for (int i = 0; i < singleSide; ++i) {
      powSpectDens[i] = powSpectDens[i].multiply(psdNormalization);
      frequencies[i] = i * deltaFreq;
    }
    
    return new FFTResult(powSpectDens, frequencies);
    
  }
  
  final private Complex[] transform; // the FFT data
  
  final private double[] freqs; // array of frequencies matching the fft data
  
  /**
   * Instantiate the structure holding an FFT and its frequency range
   * (Used to return data from the spectral density calculations)
   * Holds results of an FFT calculation already performed, usable in return
   * statements
   * @param inPSD Precalculated FFT result for some timeseries
   * @param inFreq Frequencies matched up to each FFT value
   */
  public FFTResult(Complex[] inPSD, double[] inFreq) {
    transform = inPSD;
    freqs = inFreq;
  }

  /**
   * Get the FFT for some sort of previously calculated data
   * @return Array of FFT results, as complex numbers
   */
  public Complex[] getFFT() {
    return transform;
  }
  
  /**
   * Return the value of the FFT at the given index
   * @param idx Index to get the FFT value at
   * @return FFT value at index
   */
  public Complex getFFT(int idx) {
    return transform[idx];
  }
  
  /**
   * Get the frequency value at the given index
   * @param idx Index to get the frequency value at
   * @return Frequency value at index
   */
  public double getFreq(int idx) {
    return freqs[idx];
  }
  
  /**
   * Get the frequency range for the (previously calculated) FFT
   * @return Array of frequencies (doubles), matching index to each FFT point
   */
  public double[] getFreqs() {
    return freqs;
  }
  
  /**
   * Get the size of the complex array of FFT values, also the size of the
   * double array of frequencies for the FFT at each index
   * @return int representing size of thi's object's arrays
   */
  public int size() {
    return transform.length;
  }

  public enum TaperType {
    COS, // cosine taper
    MULT; // multi-taper function
  }
  
}
