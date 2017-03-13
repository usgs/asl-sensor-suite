package asl.sensor.experiment;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.complex.Complex;
import org.jfree.data.xy.XYSeriesCollection;

import asl.sensor.input.DataBlock;
import asl.sensor.input.DataStore;
import asl.sensor.input.InstrumentResponse;
import asl.sensor.utils.FFTResult;

/**
 * Gain experiment does tests to determine a relative gain value of a sensor's
 * PSD as compared to a reference sensor over the same range of frequencies.
 * The ratio of the means of the data over that range is taken, as is the
 * standard deviation of the ratio; assuming fixed reference center gain, the
 * result of the calculated gain is given as gain2/ratio where gain2 is the
 * gain of the sensor we want to calculate (that is, not the reference sensor).
 * @author akearns
 *
 */
public class GainExperiment extends Experiment {

  /**
   * Gets the indices denoting the inclusive range of a frequency bound on
   * a list of input frequencies
   * @param freqs Frequency series resulting from FFT calculation
   * (must be pre-sorted)
   * @param freqBoundaries Array of size 2 denoting frequency 
   * upper and lower bounds (in Hz), lower bound first
   * @return Array of size 2 denoting the lower and upper indices of the 
   * sub-range of frequencies including the given boundaries
   */
  private static int[] getRange(double[] freqs, double[] freqBoundaries) {
    
    double lowFrq = freqBoundaries[0];
    double highFrq = freqBoundaries[1];
    
    int[] indices = new int[2];
    
    for (int i = 1; i < freqs.length; ++i) {
      if (freqs[i] == lowFrq || 
          (freqs[i-1] < lowFrq && freqs[i] > lowFrq) ) {
        indices[0] = i;
      }
    }
    
    indices[1] = freqs.length-1;
    
    for (int i = indices[0]; i < freqs.length-1; ++i) {
      if (freqs[i] == highFrq || 
          (freqs[i] < highFrq && freqs[i+1] > highFrq) ) {
        indices[1] = i+1;
      }
    }
    
    return indices;
  }
  
  /**
   * Get the mean of the PSD calculation within the specified range
   * @param psd PSD calculation (frequency space / FFT)
   * @param lower Starting index of window
   * @param higher Ending index of window
   * @return Mean of datapoints in the specific range
   */
  private static double mean(FFTResult psd, int lower, int higher) {
    double result = 0;
    int range = higher-lower;
    
    Complex[] data = psd.getFFT();
    
    for (int i = lower; i <= higher; ++i) {
      if ( data[i].abs() >= Double.POSITIVE_INFINITY ) {
        continue;
      }
      result += data[i].abs();
    }
    
    return result/range; // since result is a double, no cast needed?
    
  }
  
  /**
   * PSD Standard deviation calculation given mean ratio and specified range
   * (Get the result for mean calculation using the mean function before
   * calling this calculation; does not call the function directly)
   * @param fft1 first PSD (numerator) of relative gain standard deviation
   * @param fft2 second PSD (denominator)
   * @param meanRatio mean(fft1)/mean(fft2)
   * @param lower Starting index of window (should be same used in mean calc)
   * @param higher Ending index of window (should be same used in mean calc)
   * @return double corresponding to standard deviation ratio over the window
   */
  private static double sdev(
      FFTResult fft1, FFTResult fft2, 
      double meanRatio, 
      int lower, int higher) {
    Complex[] density1 = fft1.getFFT();
    Complex[] density2 = fft2.getFFT();
    double sigma = 0.;
    for (int i = lower; i <= higher; ++i) {
      double value1 = density1[i].abs();
      double value2 = density2[i].abs();
      
      if (value1 >= Double.POSITIVE_INFINITY || 
          value2 >= Double.POSITIVE_INFINITY) {
        continue;
      }
      
      sigma += Math.pow( (value1 / value2) - meanRatio, 2 );
    }
    
    return Math.sqrt(sigma);
  }
  
  private List<Double> gainStage1;
  private List<Double> otherGainStages; // product of gain stages 2 and up
  
  private List<FFTResult> fftResults;
  
  private double ratio, sigma;
  
  /**
   * Constructor for the gain experiment; effectively the same as that of the
   * superclass, as all components unique to this class are populated at time
   * of calculation
   */
  public GainExperiment() {
    super();
  }

  /**
   * Populate XYSeriesCollection with all input data (will be trimmed on plot)
   * This gets the (possibly cached) PSDs from available data in the
   * passed-in DataStore object, builds plottable series from them, and
   * generates statistics about gain calculation, defaulting to the
   * octave around the peak-value frequency.
   * The data to be plotted consists of the data from a reference series held
   * as constant and a second series for which a gain estimation is done.
   * The specific statistics being calculated are the mean and standard dev.
   * of the ratios of the series within the given frequency range, used to
   * calculate an estimation of gain from the second series.
   */
  @Override
  protected void backend(final DataStore ds, final boolean freqSpace) {
    
    int numberToLoad = ds.numberFullySet();
    
    int[] indices = new int[numberToLoad];
    
    for (int i = 1; i <= numberToLoad; ++i) {
      indices[i-1] = ds.getXthFullyLoadedIndex(i);
    }
    
    gainStage1 = new ArrayList<Double>();
    otherGainStages = new ArrayList<Double>();
    
    InstrumentResponse[] resps = ds.getResponses();
    for (int i = 0; i < 2; ++i) {
      if (resps[i] == null) {
        continue;
      }
      List<Double> gains = resps[i].getGain();
      gainStage1.add( gains.get(1) );
      double accumulator = 1.0;
      for (int j = 2; j < gains.size(); ++j) {
        accumulator *= gains.get(j);
      }
      otherGainStages.add(accumulator);
    }
    
    fftResults = new ArrayList<FFTResult>();
    ArrayList<DataBlock> blocksPlotting = new ArrayList<DataBlock>();
    
    for (int i = 0; i < indices.length; ++i) {
      int idx = indices[i];
      fftResults.add( ds.getPSD(idx) );
      blocksPlotting.add( ds.getBlock(idx) );
    }
    
    XYSeriesCollection xysc = (XYSeriesCollection) xySeriesData;
    
    xysc.setAutoWidth(true);
    
    addToPlot(ds, false, indices, xysc); 
    // false, because we don't want to plot in frequency space
    
    xysc.addSeries( FFTResult.getLowNoiseModel(false) );
    
  }
  
  /**
   * Gets the octave centered around the frequency at the plotted PSD peak
   * @param idx Index of inputted data to get the peak of
   * @return Array containing 2 elements, the values of the low and high 
   * frequencies bounding the octave
   */
  public double[] getOctaveCenteredAtPeak(int idx) {
    
    int center = getPeakIndex(idx);
    double[] freqs = fftResults.get(idx).getFreqs();
    int max = freqs.length - 1;
    double peakFreq = freqs[center];
    
    double lowFreq = Math.max( peakFreq / Math.sqrt(2), freqs[0] );
    double highFreq = Math.min( peakFreq * Math.sqrt(2), freqs[max] );
    
    return new double[]{lowFreq, highFreq};
  }
  
  /**
   * Finds the maximum value of PSD plot curve, by its index in the array
   * @param fft PSD calculation, including both FFT function and matching
   * frequencies
   * @return The index of the peak location
   */
  private int getPeakIndex(int idx) {
    
    FFTResult fft = fftResults.get(idx);
    
    Complex[] timeSeries = fft.getFFT();
    double[] freqs = fft.getFreqs();
    
    double max = Double.NEGATIVE_INFINITY;
    int index = 0;
    for (int i = 0; i < timeSeries.length; ++i) {
      if (freqs[i] < 0.001) {
        continue;
      }
      Complex temp = timeSeries[i].multiply(Math.pow(2*Math.PI*freqs[i],4));
      double result = 10*Math.log10( temp.abs() );
      if ( result < Double.POSITIVE_INFINITY && result > max ) {
        max = result;
        index = i;
      }
    }
    return index;
  }

  /**
   * Given indices to specific PSD data sets and frequency boundaries, gets
   * the mean and standard deviation ratios 
   * @param idx0 Index of numerator PSD
   * @param idx1 Index of denominator PSD
   * @param lowFreq Lower-bound of frequency window of PSD
   * @param highFreq Upper-bound of frequency window of PSD
   * @return 2-entry array of form {mean, standard deviation}
   */
  public double[] getStatsFromFreqs(int idx0, int idx1, 
      double lowFreq, double highFreq) {
    
    // TODO: check indices are valid
    
    FFTResult plot0 = fftResults.get(idx0);
    
    double[] freqBoundaries = new double[2];
    freqBoundaries[0] = Math.min(lowFreq, highFreq);
    freqBoundaries[1] = Math.max(lowFreq, highFreq);
    
    int[] indices = getRange( plot0.getFreqs(), freqBoundaries );
    
    return getStatsFromIndices(idx0, idx1, indices[0], indices[1]);
  }
  
  /**
   * Given indices to specific PSD data sets and indices to the corresponding
   * frequency boundaries, gets the mean and standard deviation ratios
   * @param idx0 Index of numerator PSD (reference sensor)
   * @param idx1 Index of denominator PSD
   * @param lowBnd Lower-bound index of PSDs' frequency array
   * @param higBnd Upper-bound index of PSDs' frequency array
   * @return 2-entry array of form {mean, standard deviation}
   */
  private double[] getStatsFromIndices(int idx0, int idx1,
      int lowBnd, int higBnd) {
    
    // make sure lowInd really is the lower index
    int temp = Math.min(lowBnd, higBnd);
    higBnd = Math.max(lowBnd, higBnd);
    lowBnd = temp;
    
    FFTResult plot0 = fftResults.get(idx0);
    FFTResult plot1 = fftResults.get(idx1);
    
    double mean0 = mean(plot0, lowBnd, higBnd);
    // since both datasets must have matching interval, PSDs have same freqs
    double mean1 = mean(plot1, lowBnd, higBnd);
    
    // double MIN_VALUE field is effectively java's machine epsilon
    // calculate ratio and sigma over the range
    ratio = (mean0+Double.MIN_VALUE) / (mean1+Double.MIN_VALUE); 
      // added terms exist to prevent division by 0
    
    sigma = sdev(plot0, plot1, ratio, lowBnd, higBnd);
    
    double gain1 = gainStage1.get(idx0);
    double gain2 = gainStage1.get(idx1)/Math.sqrt(ratio);
    
    return new double[]{Math.sqrt(ratio), sigma, gain1, gain2};
  }
  
  /**
   * Find the peak frequency of the reference series and use it to get the
   * gain statistics
   * @param idx0 Index of the reference sensor's FFT data
   * @param idx1 Index of the other sensor's FFT data
   * @return
   */
  public double[] getStatsFromPeak(int idx0, int idx1) {
    double[] freqBounds = getOctaveCenteredAtPeak(idx0);
    double freq1 = freqBounds[0];
    double freq2 = freqBounds[1];
    return getStatsFromFreqs(idx0, idx1, freq1, freq2);
  }

  @Override
  public boolean hasEnoughData(DataStore ds) {
    return ( ds.bothComponentsSet(0) && ds.bothComponentsSet(1) );
  }

}
