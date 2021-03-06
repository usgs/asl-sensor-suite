package asl.sensor.experiment;

import static asl.utils.FFTResult.spectralCalc;
import static asl.utils.NumericUtils.TAU;
import static asl.utils.NumericUtils.detrend;
import static asl.utils.NumericUtils.getFFTMean;
import static asl.utils.NumericUtils.getFFTSDev;
import static asl.utils.response.ChannelMetadata.EMBED_STRING;

import asl.sensor.input.DataStore;
import asl.utils.FFTResult;
import asl.utils.response.ChannelMetadata;
import asl.utils.response.ChannelMetadata.ResponseStageException;
import asl.utils.timeseries.DataBlock;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.Arrays;
import org.apache.commons.math3.complex.Complex;
import org.apache.log4j.Logger;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

/**
 * Gain experiment does tests to determine a relative gain value of a sensor's PSD as compared to a
 * reference sensor over the same range of frequencies. The ratio of the means of the data over that
 * range is taken, as is the standard deviation of the ratio; assuming fixed reference center gain,
 * the result of the calculated gain is given as gain2/ratio where gain2 is the gain of the sensor
 * we want to calculate (that is, not the reference sensor).
 *
 * All response stages are compared over in general use, but embedded responses don't have
 * stages specified other than pole-zero stages (or, at least, any such stages are trivially
 * defined, such as a coefficient stage with a single value of 1). In cases where an embedded
 * response is used, we skip over non PZ stages.
 *
 * Some additional details in the backend relate to how these calculations are done if the A0
 * value in the RESP file is deemed inaccurate; see {@link #ERROR_LOW_BOUND} and
 * {@link GainExperiment#backend(DataStore)}.
 *
 * @author akearns - KBRWyle
 */
public class GainExperiment extends SpectralAnalysisExperiment {

  private static final Logger logger = Logger.getLogger(GainExperiment.class);

  /**
   * Percentage error lower-bound for A0 values; if the error is above this the resp-specified
   * A0 value is not used and the experimentally-determined A0 is used instead.
   */
  private static final double ERROR_LOW_BOUND = 1.;

  /**
   * Number of sensor data inputs to load (reference sensor +
   */
  private static final int NUMBER_TO_LOAD = 2;

  /**
   * Upper bound of initial region to calculate gain estimation over (9 seconds period)
   */
  public static final double DEFAULT_UP_BOUND = 9.;

  /**
   * Lower bound of initial region to calculate gain estimation over (3 seconds period)
   */
  public static final double DEFAULT_LOW_BOUND = 3.;

  private double[] stage1Gains, A0Freqs, calcA0s, respA0s;
  private FFTResult[] fftResults;
  private int[] indices; // indices of valid data sources (i.e., 0 and 1)
  private int referenceIndex;
  private double lowPeriod, highPeriod; // initialized as DEFAULT_UP/LOW_BOUND as above
  private boolean usesOnlyOneEmbeddedResp;


  /**
   * Constructor for the gain experiment; effectively the same as that of the superclass, as all
   * components unique to this class are populated at time of calculation
   */
  public GainExperiment() {
    super();
    referenceIndex = 0;
    lowPeriod = DEFAULT_LOW_BOUND;
    highPeriod = DEFAULT_UP_BOUND;
    usesOnlyOneEmbeddedResp = false;
  }


  String getResultString() {
    double[] varResultArray = getStatsFromFreqs();

    double mean = varResultArray[0];
    double standardDeviation = varResultArray[1];
    double referenceGain = varResultArray[2];
    double calculatedGain = varResultArray[3];
    double referenceFrequency = varResultArray[4];
    double calculatedFrequency = varResultArray[5];
    double refErrorA0 = varResultArray[6];
    double testErrorA0 = varResultArray[7];

    // Scientific notation format for printing experimentally-derived A0
    // will allow for printing both very large and very small values of A0
    NumberFormat sciFormat = new DecimalFormat("0.00000E00");
    NumberFormat trimFormat = new DecimalFormat("00.00E0");

    StringBuilder sb = new StringBuilder("Ratio: " + DECIMAL_FORMAT.get().format(mean)
        + "\nSigma: " + DECIMAL_FORMAT.get().format(standardDeviation)
        + "\nRef. gain: " + DECIMAL_FORMAT.get().format(referenceGain)
        + " [resp. A0 " + sciFormat.format(referenceFrequency) + "Hz]"
        + "\n** CALCULATED GAIN: " + DECIMAL_FORMAT.get().format(calculatedGain)
        + " [resp. A0 " + sciFormat.format(calculatedFrequency) + "Hz]");

    if (refErrorA0 > ERROR_LOW_BOUND) {
      // these values are not part of the result if the percent error is 1% or lower
      double calcA0 = varResultArray[8];
      sb.append("\n**REF SENSOR A0 VALUE ERROR ")
          .append(trimFormat.format(refErrorA0)).append("%**")
          .append("\nUsed estimated A0: ").append(sciFormat.format(calcA0));
    }
    if (testErrorA0 > ERROR_LOW_BOUND) {
      double calcA0 = varResultArray[9];
      sb.append("\n**TEST SENSOR A0 VALUE ERROR ")
          .append(trimFormat.format(testErrorA0)).append("%**")
          .append("\nUsed estimated A0: ").append(sciFormat.format(calcA0));
    }

    return sb.toString();
  }

  @Override
  public String[] getDataStrings() {
    return new String[]{getResultString()};
  }

  @Override
  public String[] getInsetStrings() {
    return new String[]{getResultString() + '\n' + getFormattedDateRange()};
  }

  /**
   * Populate XYSeriesCollection with all input data (will be trimmed on plot) This gets the
   * (possibly cached) PSDs from available data in the passed-in DataStore object, builds plottable
   * series from them, and generates statistics about gain calculation, defaulting to the octave
   * around the peak-value frequency. The data to be plotted consists of the data from a reference
   * series held as constant and a second series for which a gain estimation is done. The specific
   * statistics being calculated are the mean and standard dev. of the ratios of the series within
   * the given frequency range, used to calculate an estimation of gain from the second series.
   *
   * Note that because the gain values are dependent both on the actual gain stages specified
   * and the pole-zero stage's A0 value, in the event that a RESP file has a bad A0 parameter,
   * the calculations are performed on an experimentally derived A0 value (done by getting
   * 1/r(A0f) where r() is the unscaled response curve and A0f is the value defined in the resp
   * for the A0 frequency. This value is applied both to the results published in the index and
   * also to the plotted results; the inset will include text warning the user this has been done.
   */
  @Override
  protected void backend(final DataStore dataStore) {

    usesOnlyOneEmbeddedResp = false;

    lowPeriod = DEFAULT_LOW_BOUND;
    highPeriod = DEFAULT_UP_BOUND;

    indices = new int[NUMBER_TO_LOAD];
    // indices here is a linear array pointing to where
    // the data is in the passed-in data store, mainly to find relevant resps

    int maxLength = Integer.MAX_VALUE;

    for (int i = 0; i < NUMBER_TO_LOAD; ++i) {
      // XthFullyLoaded starts at 1 (i.e., get first full-loaded), not 0
      int idx = dataStore.getXthFullyLoadedIndex(i + 1);
      indices[i] = idx;
      dataNames.add(dataStore.getBlock(idx).getName());
      dataNames.add(dataStore.getResponse(idx).getName());
      maxLength = Math.min(maxLength, dataStore.getBlock(idx).size());
      // XOR equivalent -- if embedded response is found and bool is false, make true
      // else if embedded response is found and bool is *already* true, make false
      if (dataStore.getResponse(idx).getFullFilePath().equals(EMBED_STRING)) {
        usesOnlyOneEmbeddedResp = !usesOnlyOneEmbeddedResp;
      }
    }

    stage1Gains = new double[NUMBER_TO_LOAD];
    A0Freqs = new double[NUMBER_TO_LOAD];
    calcA0s = new double[NUMBER_TO_LOAD];
    respA0s = new double[NUMBER_TO_LOAD];

    fireStateChange("Accumulating gain values...");
    // InstrumentResponse[] resps = ds.getResponses();
    for (int i = 0; i < indices.length; ++i) {
      ChannelMetadata ir = dataStore.getResponse(indices[i]);
      {
        double[] gains = ir.getGain();
        stage1Gains[i] = gains[1];
      }
      try {
        A0Freqs[i] = ir.getPoleZeroStage().getNormalizationFreq();
        respA0s[i] = ir.getPoleZeroStage().getNormalizationFactor();
        calcA0s[i] = ir.getEstimatedA0();
      } catch (ResponseStageException e) {
        // this probably shouldn't happen
        String errMsg = "There is no Pole-Zero response stage for " + ir.getName();
        throw new RuntimeException(errMsg);
      }
    }

    fftResults = new FFTResult[NUMBER_TO_LOAD];
    XYSeriesCollection xysc = new XYSeriesCollection();
    xysc.setAutoWidth(true);

    if (usesOnlyOneEmbeddedResp) {
      logger.warn("Calculating PSDs using only response stage 1 due to comparison with embed");
    }

    for (int i = 0; i < indices.length; ++i) {
      fireStateChange("Getting PSD " + i + "...");
      int idx = indices[i];
      String name = "PSD " + dataStore.getBlock(idx).getName() + " [" + idx + "]";
      XYSeries xys = new XYSeries(name);
      if (!usesOnlyOneEmbeddedResp) {
        fftResults[i] = dataStore.getPSD(idx, maxLength);
      } else {
        DataBlock block = dataStore.getBlock(idx);
        double[] timeSeries = Arrays.copyOfRange(dataStore.getBlock(idx).getData(), 0, maxLength);
        timeSeries = detrend(timeSeries);
        FFTResult uncorrectedPSD = spectralCalc(timeSeries, timeSeries, block.getInterval());
        double[] freqs = uncorrectedPSD.getFreqs();
        // note that maxStage here is set to 2 as it is an exclusive upper bound
        Complex[] responseCurve =
            dataStore.getResponse(idx).applyResponseToInput(freqs, 1, 2);
        Complex[] out = new Complex[freqs.length];
        for (int j = 0; j < freqs.length; ++j) {
          // response curves in velocity, put them into acceleration
          Complex scaleFactor =
              new Complex(0., -1.0 / (TAU * freqs[j]));
          Complex resp1 = responseCurve[j].multiply(scaleFactor);
          Complex respMagnitude =
              resp1.multiply(resp1.conjugate());

          if (respMagnitude.abs() == 0) {
            respMagnitude = new Complex(Double.MIN_VALUE, 0);
          }

          out[j] = uncorrectedPSD.getFFT(j).divide(respMagnitude);
        }
        fftResults[i] = new FFTResult(out, freqs);
      }
      Complex[] fft = fftResults[i].getFFT();
      double[] freqs = fftResults[i].getFreqs();
      double errorOnA0 = Math.abs(calcA0s[i] - respA0s[i]) / Math.abs(respA0s[i]) * 100.;
      // an important point to consider here -- if the A0 error is high enough
      // we are recalculating the A0 for the PSD for the given sensor
      // and re-assigning that value back to the FFT array we started from
      // so if we do any future statistics on this data we don't need to re-correct for A0
      // because this value is what's pointed to from fftResults array for this index
      if (errorOnA0 > ERROR_LOW_BOUND) {
        for (int j = 0; j < fft.length; ++j) {
          // since A0 is a property of the response curve gain, we divide it out of the FFT
          double respA0Sqd = Math.pow(respA0s[i], 2);
          double calcA0Sqd = Math.pow(calcA0s[i], 2);
          fft[j] = fft[j].multiply(respA0Sqd).divide(calcA0Sqd);
        }
      }
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
   * @return Array containing 2 elements, the values of the low and high frequencies bounding the
   * octave
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
   * Given indices to specific PSD data sets and frequency boundaries, gets the mean and standard
   * deviation ratios. For more details on the specifics of this calculation see the linked
   * method where the values are actually calculated.
   *
   * @see #getStatsFromIndices(int, int, int)
   * @return Array of form {mean, standard deviation, ref. gain, calc. gain, ref. A0 freq., calc. A0
   * freq., ref. A0 error, calc. A0 error, ref. sensor's exp. A0, calc. sensor's exp. A0}
   */
  public double[] getStatsFromFreqs() {
    FFTResult plot0 = fftResults[referenceIndex];
    double lowerBound = 1. / highPeriod; // high period = low frequency
    double upperBound = 1. / lowPeriod; // low period = high frequency

    int lowIndex =
        FFTResult.getIndexOfFrequency(plot0.getFreqs(), lowerBound);
    int highIndex =
        FFTResult.getIndexOfFrequency(plot0.getFreqs(), upperBound);
    return getStatsFromIndices(referenceIndex, lowIndex, highIndex);
  }

  /**
   * Find the peak frequency of the reference series and use it to get the gain statistics. Only
   * used as a quick reference call from an automated test.
   *
   * @see #getStatsFromIndices(int, int, int)
   * @param refIndex Index of the reference sensor's FFT data
   * @return Array of form {mean, standard deviation, ref. gain, calc. gain, ref. A0 freq., calc. A0
   * freq., ref. A0 error, calc. A0 error, ref. sensor's exp. A0, calc. sensor's exp. A0}
   */
  public double[] getStatsFromPeak(int refIndex) {
    double[] freqBounds = getOctaveCenteredAtPeak(refIndex);
    setReferenceIndex(refIndex);
    setRangeForStatistics(1. / freqBounds[0], 1. / freqBounds[1]);
    return getStatsFromFreqs();
  }

  /**
   * Given indices to specific PSD data sets and indices to the corresponding frequency boundaries,
   * gets the mean and standard deviation ratios as well as the frequencies the normalizations are
   * taken from. The percent error on the RESP vs. experimental A0 is also provided, along with what
   * the experimental A0 values are for each sensor.
   *
   * Note that if an A0 error value is above {@link #ERROR_LOW_BOUND}, the calculations
   * for gain, etc. (along with the plotted data) will use the experimental A0 value.
   *
   * @param refIndex Index of first curve to be plotted (numerator PSD)
   * @param lowerBound Lower-bound index of PSDs' frequency array
   * @param upperBound Upper-bound index of PSDs' frequency array
   * @return Array of form {mean, standard deviation, ref. gain, calc. gain, ref. A0 freq., calc. A0
   * freq., ref. A0 error, calc. A0 error, ref. sensor's exp. A0, calc. sensor's exp. A0}
   */
  private double[] getStatsFromIndices(int refIndex, int lowerBound, int upperBound) {

    int refIndexPlusOne = (refIndex + 1) % NUMBER_TO_LOAD;

    // make sure lower bound really is the lower bound
    int temp = Math.min(lowerBound, upperBound);
    upperBound = Math.max(lowerBound, upperBound);
    lowerBound = temp;

    FFTResult plot0 = fftResults[refIndex];
    FFTResult plot1 = fftResults[refIndexPlusOne];

    double mean0 = getFFTMean(plot0, lowerBound, upperBound);
    // since both datasets must have matching interval, PSDs have same frequencies
    double mean1 = getFFTMean(plot1, lowerBound, upperBound);

    // double MIN_VALUE field is effectively java's machine epsilon
    // calculate ratio and sigma over the range
    double ratio = (mean0 + Double.MIN_VALUE) / (mean1 + Double.MIN_VALUE);
    // added terms exist to prevent division by 0

    double sigma = getFFTSDev(plot0, plot1, ratio, lowerBound, upperBound);

    double refGain = stage1Gains[refIndex];
    double calcGain = stage1Gains[refIndexPlusOne] / Math.sqrt(ratio);

    double normalFreqRef = A0Freqs[refIndex];
    double normalFreqTest = A0Freqs[refIndexPlusOne];

    double calcA0Ref = calcA0s[refIndex];
    double respA0Ref = respA0s[refIndex];
    double errorOnA0Ref = Math.abs(calcA0Ref - respA0Ref) / Math.abs(respA0Ref) * 100.;

    double calcA0Test = calcA0s[refIndexPlusOne];
    double respA0Test = respA0s[refIndexPlusOne];
    double errorOnA0Test = Math.abs(calcA0Test - respA0Test) / Math.abs(respA0Test) * 100.;

    return new double[]{
        Math.sqrt(ratio), sigma, refGain, calcGain, normalFreqRef, normalFreqTest,
        errorOnA0Ref, errorOnA0Test, calcA0Ref, calcA0Test
    };
  }

  @Override
  public boolean hasEnoughData(DataStore dataStore) {
    return (dataStore.bothComponentsSet(0) && dataStore.bothComponentsSet(1));
  }

  @Override
  public int[] listActiveResponseIndices() {
    return indices;
  }

  /**
   * Choose what index of data to hold as reference for gain estimation, either the first input (0)
   * or second (1). If a number other than this is chosen, the correct value mod {@link
   * #NUMBER_TO_LOAD NUMBER_TO_LOAD} will be chosen.
   *
   * @param newIndex Index of data to be used as reference
   */
  public void setReferenceIndex(int newIndex) {
    // reference index must be positive
    referenceIndex = ((newIndex % NUMBER_TO_LOAD) + NUMBER_TO_LOAD) % NUMBER_TO_LOAD;
  }

  /**
   * Select the range of data to be used as the range over which gain statistics are calculated.
   * This is by default the range from 3 to 9 seconds period.
   *
   * @param lowPeriod New low period to get range over
   * @param highPeriod New high period to get range over
   */
  public void setRangeForStatistics(double lowPeriod, double highPeriod) {
    this.highPeriod = Math.max(lowPeriod, highPeriod);
    this.lowPeriod = Math.min(lowPeriod, highPeriod);
  }

  /**
   * This is true if a comparison is being made between an embedded response (which has no FIR
   * stages) and a field response (which does have them). In order to produce a reasonable
   * comparison, the field response's FIR stages are ignored during processing.
   * @return true if such condition is met
   */
  public boolean skipsFIRStages() {
    return usesOnlyOneEmbeddedResp;
  }
}
