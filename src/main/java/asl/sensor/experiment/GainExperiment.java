package asl.sensor.experiment;

import asl.sensor.input.DataStore;
import asl.sensor.input.InstrumentResponse;
import asl.sensor.utils.FFTResult;
import asl.sensor.utils.NumericUtils;
import org.apache.commons.math3.complex.Complex;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

/**
 * Gain experiment does tests to determine a relative gain value of a sensor's
 * PSD as compared to a reference sensor over the same range of frequencies.
 * The ratio of the means of the data over that range is taken, as is the
 * standard deviation of the ratio; assuming fixed reference center gain, the
 * result of the calculated gain is given as gain2/ratio where gain2 is the
 * gain of the sensor we want to calculate (that is, not the reference sensor).
 *
 * @author akearns - KBRWyle
 */
public class GainExperiment extends Experiment {

  private static final int NUMBER_TO_LOAD = 2;

  private double[] gainStage1;
  private FFTResult[] fftResults;
  private int[] indices; // indices of valid data sources (i.e., 0 and 1)

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
  protected void backend(final DataStore dataStore) {

    indices = new int[NUMBER_TO_LOAD];
    // indices here is a linear array pointing to where
    // the data is in the passed-in data store, mainly to find relevant resps

    for (int i = 0; i < NUMBER_TO_LOAD; ++i) {
      // XthFullyLoaded starts at 1 (i.e., get first full-loaded), not 0
      int idx = dataStore.getXthFullyLoadedIndex(i + 1);
      indices[i] = idx;
      dataNames.add(dataStore.getBlock(idx).getName());
      dataNames.add(dataStore.getResponse(idx).getName());
    }

    gainStage1 = new double[NUMBER_TO_LOAD];

    fireStateChange("Accumulating gain values...");
    // InstrumentResponse[] resps = ds.getResponses();
    for (int i = 0; i < indices.length; ++i) {
      InstrumentResponse ir = dataStore.getResponse(indices[i]);
      double[] gains = ir.getGain();
      gainStage1[i] = gains[1];
    }

    fftResults = new FFTResult[NUMBER_TO_LOAD];
    XYSeriesCollection xysc = new XYSeriesCollection();
    xysc.setAutoWidth(true);

    for (int i = 0; i < indices.length; ++i) {
      fireStateChange("Getting PSD " + i + "...");
      int idx = indices[i];
      String name = "PSD " + dataStore.getBlock(idx).getName() + " [" + idx + "]";
      XYSeries xys = new XYSeries(name);
      fftResults[i] = dataStore.getPSD(idx);
      Complex[] fft = fftResults[i].getFFT();
      double[] freqs = fftResults[i].getFreqs();
      // false, because we don't want to plot in frequency space
      addToPlot(xys, fft, freqs, false, xysc);
    }

    fireStateChange("Getting NLNM data...");
    xysc.addSeries(FFTResult.getLowNoiseModel(false));

    xySeriesData.add(xysc);
  }

  @Override
  public int blocksNeeded() {
    return 2;
  }

  /**
   * Gets the octave centered around the frequency at the plotted PSD peak
   *
   * @param index Index of inputted data to get the peak of
   * @return Array containing 2 elements, the values of the low and high
   * frequencies bounding the octave
   */
  private double[] getOctaveCenteredAtPeak(int index) {

    int center = getPeakIndex(index);
    double[] freqs = fftResults[index].getFreqs();
    int max = freqs.length - 1;
    double peakFreq = freqs[center];

    double lowFreq = Math.max(peakFreq / Math.sqrt(2), 0.001);
    double highFreq = Math.min(peakFreq * Math.sqrt(2), freqs[max]);

    return new double[]{lowFreq, highFreq};
  }

  /**
   * Finds the maximum value of PSD plot curve, by its index in the array
   *
   * @param fftIndex Index of array to be loaded in from result set to be plotted
   * @return The index of the peak location
   */
  private int getPeakIndex(int fftIndex) {

    FFTResult fft = fftResults[fftIndex];

    Complex[] timeSeries = fft.getFFT();
    double[] freqs = fft.getFreqs();

    double max = Double.NEGATIVE_INFINITY;
    int peakIndex = 0;
    for (int i = 0; i < timeSeries.length; ++i) {
      if (freqs[i] < 0.001) {
        continue;
      }

      double result = 10 * Math.log10(timeSeries[i].abs());
      if (result < Double.POSITIVE_INFINITY && result > max) {
        max = result;
        peakIndex = i;
      }
    }
    return peakIndex;
  }

  /**
   * Given indices to specific PSD data sets and frequency boundaries, gets
   * the mean and standard deviation ratios
   *
   * @param refIndex Index of first curve to be plotted (numerator PSD)
   * @param lowerBound Lower-bound of frequency window of PSD
   * @param upperBound Upper-bound of frequency window of PSD
   * @return Array of form {mean, standard deviation, ref. gain, calc. gain}
   */
  public double[] getStatsFromFreqs(int refIndex, double lowerBound, double upperBound) {
    FFTResult plot0 = fftResults[refIndex];

    int lowIndex =
        FFTResult.getIndexOfFrequency(plot0.getFreqs(), Math.min(lowerBound, upperBound));
    int highIndex =
        FFTResult.getIndexOfFrequency(plot0.getFreqs(), Math.max(lowerBound, upperBound));
    return getStatsFromIndices(refIndex, lowIndex, highIndex);
  }

  /**
   * Given indices to specific PSD data sets and indices to the corresponding
   * frequency boundaries, gets the mean and standard deviation ratios
   *
   * @param refIndex Index of first curve to be plotted (numerator PSD)
   * @param lowerBound Lower-bound index of PSDs' frequency array
   * @param upperBound Upper-bound index of PSDs' frequency array
   * @return Array of form {mean, standard deviation, ref. gain, calc. gain}
   */
  private double[] getStatsFromIndices(int refIndex, int lowerBound, int upperBound) {

    int refIndexPlusOne = (refIndex + 1) % NUMBER_TO_LOAD;

    // make sure lowInd really is the lower index
    int temp = Math.min(lowerBound, upperBound);
    upperBound = Math.max(lowerBound, upperBound);
    lowerBound = temp;

    FFTResult plot0 = fftResults[refIndex];
    FFTResult plot1 = fftResults[refIndexPlusOne];

    double mean0 = NumericUtils.getFFTMean(plot0, lowerBound, upperBound);
    // since both datasets must have matching interval, PSDs have same frequencies
    double mean1 = NumericUtils.getFFTMean(plot1, lowerBound, upperBound);

    // double MIN_VALUE field is effectively java's machine epsilon
    // calculate ratio and sigma over the range
    double ratio = (mean0 + Double.MIN_VALUE) / (mean1 + Double.MIN_VALUE);
    // added terms exist to prevent division by 0

    double sigma = NumericUtils.getFFTSDev(plot0, plot1, ratio, lowerBound, upperBound);

    double refGain = gainStage1[refIndex];
    double calcGain = gainStage1[refIndexPlusOne] / Math.sqrt(ratio);

    return new double[]{Math.sqrt(ratio), sigma, refGain, calcGain};
  }

  /**
   * Find the peak frequency of the reference series and use it to get the
   * gain statistics
   *
   * @param refIndex Index of the reference sensor's FFT data
   * @return Array of form {mean, standard deviation, ref. gain, calc. gain}
   */
  public double[] getStatsFromPeak(int refIndex) {
    double[] freqBounds = getOctaveCenteredAtPeak(refIndex);
    return getStatsFromFreqs(refIndex, freqBounds[0], freqBounds[1]);
  }

  @Override
  public boolean hasEnoughData(DataStore dataStore) {
    return (dataStore.bothComponentsSet(0) && dataStore.bothComponentsSet(1));
  }

  @Override
  public int[] listActiveResponseIndices() {
    return indices;
  }

}
