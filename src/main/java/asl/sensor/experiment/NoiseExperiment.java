package asl.sensor.experiment;

import asl.sensor.input.DataStore;
import asl.utils.FFTResult;
import asl.utils.input.DataBlock;
import asl.utils.input.InstrumentResponse;
import java.util.Arrays;
import java.util.stream.IntStream;
import org.apache.commons.math3.complex.Complex;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

/**
 * Produces the data for a self-noise test. Calculates PSD to get cross-power.
 * These calculations are based around the 3-component self noise calculation
 * rather than the 9-component.
 * Based on code in the seedscan timeseries package, see
 * https://github.com/usgs/seedscan/tree/master/src/main/java/asl/timeseries
 * See also Ringler, Hutt: 'Self-Noise Models of Seismic Instruments', Seismological Research
 * Letters (2010).
 * The self-noise formula is derived from Sleeman, Wettum, Trampert: 'Three-Channel Correlation
 * Anlysis: A New Technique to Measure Instrumental Noise of Digitizers and Seismic Sensors' (2006).
 * @author akearns - KBRWyle
 * @author jholland - USGS
 */
public class NoiseExperiment extends Experiment {


  /**
   * Controls plotting in Hz vs. time interval between samples
   */
  boolean freqSpace;

  protected static final int DATA_NEEDED = 3;

  /**
   * To keep track of the response data used in this experiment
   */
  int[] respIndices;

  /**
   * Instantiates a noise experiment -- axis titles and scales
   */
  public NoiseExperiment() {
    super();
    respIndices = new int[DATA_NEEDED];
    freqSpace = false;
  }

  @Override
  public String[] getInsetStrings() {
    return new String[]{getFormattedDateRange()};
  }

  /**
   * Generates power spectral density of each inputted file, and calculates
   * self-noise based on that result.
   * The overhead view is as follows:
   * Take a window of size 1/4 incrementing through 1/16 of the data and
   * calculate the FFT of that region. Average these results together.
   * Apply the magnitude of the frequency response (relative to the FFT indices)
   * to that result and then take the complex conjugate.
   * This produces the PSD plots.
   * Then, take the cross-powers of each of the terms (same calculation, but
   * multiply one result by the complex conjugate of the other), producing the
   * remaining terms for the formula for the self-noise results.
   */
  @Override
  protected void backend(final DataStore dataStore) {

    XYSeriesCollection xysc = new XYSeriesCollection();
    xysc.setAutoWidth(true);

    // there are 3 inputs required in order to do this calculation correctly
    respIndices = new int[DATA_NEEDED]; // first 3 fully-loaded data sets

    // get the first (index.length) seed/resp pairs. while we expect to
    // have the first three plots be the ones with loaded data, in general
    // it is probably better to keep the program flexible against valid input
    int maxLength = Integer.MAX_VALUE;
    for (int i = 0; i < respIndices.length; ++i) {
      // xth fully loaded function begins at 1
      int idx = dataStore.getXthFullyLoadedIndex(i + 1);
      respIndices[i] = idx;
      dataNames.add(dataStore.getBlock(idx).getName());
      maxLength = Math.min(maxLength, dataStore.getBlock(idx).size());
      dataNames.add(dataStore.getResponse(idx).getName());
    }

    DataBlock[] dataIn = new DataBlock[respIndices.length];
    InstrumentResponse[] responses = new InstrumentResponse[respIndices.length];

    for (int i = 0; i < respIndices.length; ++i) {
      dataIn[i] = dataStore.getBlock(respIndices[i]);
      responses[i] = dataStore.getResponse(respIndices[i]);
    }

    Complex[][] spectra = new Complex[DATA_NEEDED][];
    double[] freqs = new double[]{}; // initialize to prevent later errors

    // gets the PSDs of each given index for given freqSpace
    for (int i = 0; i < respIndices.length; ++i) {
      int idx = respIndices[i];
      fireStateChange("Getting PSDs of data " + (idx + 1) + "...");
      String name = "PSD " + dataStore.getBlock(idx).getName() + " [" + idx + "]";
      XYSeries powerSeries = new XYSeries(name);
      FFTResult psdCalc = dataStore.getPSD(idx, maxLength);
      Complex[] fft = psdCalc.getFFT();
      spectra[i] = fft;
      freqs = psdCalc.getFreqs();
      addToPlot(powerSeries, fft, freqs, freqSpace, xysc);
    }

    String getting = "Getting crosspower of series ";

    // spectra[i] is crosspower pii, now to get pij terms for i!=j
    fireStateChange(getting + "1 & 3");
    FFTResult fft =
        FFTResult.crossPower(dataIn[0], dataIn[2], responses[0], responses[2], maxLength);
    Complex[] c13 = fft.getFFT();

    fireStateChange(getting + "2 & 1");
    fft =
        FFTResult.crossPower(dataIn[1], dataIn[0], responses[1], responses[0], maxLength);
    Complex[] c21 = fft.getFFT();

    fireStateChange(getting + "2 & 3");
    fft =
        FFTResult.crossPower(dataIn[1], dataIn[2], responses[1], responses[2], maxLength);
    Complex[] c23 = fft.getFFT();

    // WIP: use PSD results to get noise at each point see spectra
    XYSeries[] noiseSeriesArr = new XYSeries[DATA_NEEDED];
    for (int j = 0; j < DATA_NEEDED; ++j) {
      // initialize each xySeries with proper name for the data
      noiseSeriesArr[j] =
          new XYSeries("Noise " + dataIn[j].getName() + " [" + j + "]");
    }

    fireStateChange("Doing noise estimation calculations...");

    Complex[][] noise = new Complex[DATA_NEEDED][];
    for (int i = 0; i < noise.length; ++i) {
      noise[i] = new Complex[spectra[0].length];
      noise[i][0] = Complex.ZERO;
    }

    final double[] finalFreqs = freqs;
    for (int i = 1; i < freqs.length; ++i) {
      if ((1. / finalFreqs[i]) > MAX_PLOT_PERIOD) {
        return;
      }

      Complex p11 = spectra[0][i];
      Complex p22 = spectra[1][i];
      Complex p33 = spectra[2][i];

      Complex p13 = c13[i];
      Complex p21 = c21[i];
      Complex p23 = c23[i];

      // n_ii = p_ii - p_ji * h_ij
      // where h_ij = p_ik / p_jk
      // and p_ik is the complex conjugate of p_ki
      noise[0][i] =
          p11.subtract(p21.multiply(p13).divide(p23)); // n_11

      noise[1][i] =
          p22.subtract(
              (p23.conjugate()).multiply(p21).divide(p13.conjugate())); // n_22

      noise[2][i] =
          p33.subtract(
              p23.multiply(p13.conjugate()).divide(p21)); //n_33

    }

    // now that we've calculated the noise, add it to the plot w/ smoothing applied
    double[] noiseRange = Arrays.copyOfRange(freqs, 1, freqs.length);
    for (int i = 0; i < noiseSeriesArr.length; ++i) {
      Complex[] noiseData = noise[i];
      XYSeries noiseSeries = noiseSeriesArr[i];
      addToPlot(noiseSeries, noiseData, noiseRange, freqSpace, xysc);
    }

    xysc.addSeries(FFTResult.getLowNoiseModel(freqSpace));
    xysc.addSeries(FFTResult.getHighNoiseModel(freqSpace));

    xySeriesData.add(xysc);

  }

  @Override
  public int blocksNeeded() {
    return DATA_NEEDED;
  }

  @Override
  public boolean hasEnoughData(DataStore dataStore) {
    for (int i = 0; i < blocksNeeded(); ++i) {
      if (!dataStore.bothComponentsSet(i)) {
        return false;
      }
    }
    return true;
  }

  /**
   * NOTE: not used by corresponding panel, overrides with active indices of components in the combo-box
   *
   * @return response indices
   */
  @Override
  public int[] listActiveResponseIndices() {
    return respIndices;
  }

  /**
   * Used to set the x-axis over which the PSDs / cross-powers are plotted,
   * either frequency (Hz) units or sample-interval (s) units
   *
   * @param freqSpace True if the plot should use units of Hz
   */
  public void setFreqSpace(boolean freqSpace) {
    this.freqSpace = freqSpace;
  }

}
