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
 * @author akearns
 */
public class GainExperiment extends Experiment {

  private static final int NUMBER_TO_LOAD = 2;

  /**
   * Gets the indices denoting the inclusive range of a frequency bound on
   * a list of input frequencies
   *
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
          (freqs[i - 1] < lowFrq && freqs[i] > lowFrq)) {
        indices[0] = i;
      }
    }

    indices[1] = freqs.length - 1;

    for (int i = indices[0]; i < freqs.length - 1; ++i) {
      if (freqs[i] == highFrq ||
          (freqs[i] < highFrq && freqs[i + 1] > highFrq)) {
        indices[1] = i + 1;
      }
    }

    return indices;
  }

  private double[] gainStage1;
  private FFTResult[] fftResults;
  private int[] indices; // indices of valid data sources (i.e., 0 and 1)
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
  protected void backend(final DataStore ds) {

    indices = new int[NUMBER_TO_LOAD];
    // indices here is a linear array pointing to where
    // the data is in the passed-in data store, mainly to find relevant resps

    for (int i = 0; i < NUMBER_TO_LOAD; ++i) {
      // XthFullyLoaded starts at 1 (i.e., get first full-loaded), not 0
      int idx = ds.getXthFullyLoadedIndex(i + 1);
      indices[i] = idx;
      dataNames.add(ds.getBlock(idx).getName());
      dataNames.add(ds.getResponse(idx).getName());
    }

    gainStage1 = new double[NUMBER_TO_LOAD];

    fireStateChange("Accumulating gain values...");
    // InstrumentResponse[] resps = ds.getResponses();
    for (int i = 0; i < indices.length; ++i) {
      InstrumentResponse ir = ds.getResponse(indices[i]);
      double[] gains = ir.getGain();
      gainStage1[i] = gains[1];
      double accumulator = 1.0;
      for (int j = 2; j < gains.length; ++j) {
        accumulator *= gains[j];
      }
    }

    fftResults = new FFTResult[NUMBER_TO_LOAD];
    XYSeriesCollection xysc = new XYSeriesCollection();
    xysc.setAutoWidth(true);

    for (int i = 0; i < indices.length; ++i) {
      fireStateChange("Getting PSD " + i + "...");
      int idx = indices[i];
      String name = "PSD " + ds.getBlock(idx).getName() + " [" + idx + "]";
      XYSeries xys = new XYSeries(name);
      fftResults[i] = ds.getPSD(idx);
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
   * @param idx Index of inputted data to get the peak of
   * @return Array containing 2 elements, the values of the low and high
   * frequencies bounding the octave
   */
  public double[] getOctaveCenteredAtPeak(int idx) {

    int center = getPeakIndex(idx);
    double[] freqs = fftResults[idx].getFreqs();
    int max = freqs.length - 1;
    double peakFreq = freqs[center];

    double lowFreq = Math.max(peakFreq / Math.sqrt(2), 0.001);
    double highFreq = Math.min(peakFreq * Math.sqrt(2), freqs[max]);

    return new double[]{lowFreq, highFreq};
  }

  /**
   * Finds the maximum value of PSD plot curve, by its index in the array
   *
   * @param idx Index of array to be loaded in from result set to be plotted
   * @return The index of the peak location
   */
  private int getPeakIndex(int idx) {

    FFTResult fft = fftResults[idx];

    Complex[] timeSeries = fft.getFFT();
    double[] freqs = fft.getFreqs();

    double max = Double.NEGATIVE_INFINITY;
    int index = 0;
    for (int i = 0; i < timeSeries.length; ++i) {
      if (freqs[i] < 0.001) {
        continue;
      }

      double result = 10 * Math.log10(timeSeries[i].abs());
      if (result < Double.POSITIVE_INFINITY && result > max) {
        max = result;
        index = i;
      }
    }
    return index;
  }

  /**
   * Given indices to specific PSD data sets and frequency boundaries, gets
   * the mean and standard deviation ratios
   *
   * @param refIdx Index of first curve to be plotted (numerator PSD)
   * @param lowFq Lower-bound of frequency window of PSD
   * @param highFq Upper-bound of frequency window of PSD
   * @return Array of form {mean, standard deviation, ref. gain, calc. gain}
   */
  public double[] getStatsFromFreqs(int refIdx, double lowFq, double highFq) {

    FFTResult plot0 = fftResults[refIdx];

    double[] freqBoundaries = new double[2];
    freqBoundaries[0] = Math.min(lowFq, highFq);
    freqBoundaries[1] = Math.max(lowFq, highFq);

    int[] indices = getRange(plot0.getFreqs(), freqBoundaries);

    return getStatsFromIndices(refIdx, indices[0], indices[1]);
  }

  /**
   * Given indices to specific PSD data sets and indices to the corresponding
   * frequency boundaries, gets the mean and standard deviation ratios
   *
   * @param refIdx Index of first curve to be plotted (numerator PSD)
   * @param lowBnd Lower-bound index of PSDs' frequency array
   * @param higBnd Upper-bound index of PSDs' frequency array
   * @return Array of form {mean, standard deviation, ref. gain, calc. gain}
   */
  private double[] getStatsFromIndices(int refIdx, int lowBnd, int higBnd) {

    int idx0 = refIdx;
    int idx1 = (refIdx + 1) % NUMBER_TO_LOAD;

    // make sure lowInd really is the lower index
    int temp = Math.min(lowBnd, higBnd);
    higBnd = Math.max(lowBnd, higBnd);
    lowBnd = temp;

    FFTResult plot0 = fftResults[idx0];
    FFTResult plot1 = fftResults[idx1];

    double mean0 = NumericUtils.getFFTMean(plot0, lowBnd, higBnd);
    // since both datasets must have matching interval, PSDs have same freqs
    double mean1 = NumericUtils.getFFTMean(plot1, lowBnd, higBnd);

    // double MIN_VALUE field is effectively java's machine epsilon
    // calculate ratio and sigma over the range
    ratio = (mean0 + Double.MIN_VALUE) / (mean1 + Double.MIN_VALUE);
    // added terms exist to prevent division by 0

    sigma = NumericUtils.getFFTSDev(plot0, plot1, ratio, lowBnd, higBnd);

    double refGain = gainStage1[idx0];
    double calcGain = gainStage1[idx1] / Math.sqrt(ratio);

    return new double[]{Math.sqrt(ratio), sigma, refGain, calcGain};
  }

  /**
   * Find the peak frequency of the reference series and use it to get the
   * gain statistics
   *
   * @param refIdx Index of the reference sensor's FFT data
   * @return Array of form {mean, standard deviation, ref. gain, calc. gain}
   */
  public double[] getStatsFromPeak(int refIdx) {
    double[] freqBounds = getOctaveCenteredAtPeak(refIdx);
    double freq1 = freqBounds[0];
    double freq2 = freqBounds[1];
    return getStatsFromFreqs(refIdx, freq1, freq2);
  }

  @Override
  public boolean hasEnoughData(DataStore ds) {
    return (ds.bothComponentsSet(0) && ds.bothComponentsSet(1));
  }

  @Override
  public int[] listActiveResponseIndices() {
    return indices;
  }

}
