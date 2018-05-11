package asl.sensor.experiment;

import asl.sensor.input.DataStore;
import asl.sensor.utils.FFTResult;
import org.jfree.data.xy.XYSeriesCollection;

/**
 * Calculates PSD to get cross-power. Can be done on 1-3 different series of data; does not estimate
 * noise parameters taken from these PSDs due to having such limited data.
 * Based on code in the seedscan timeseries package, see
 * https://github.com/usgs/seedscan/tree/master/src/main/java/asl/timeseries
 *
 * @author akearns - KBRWyle
 * @author jholland - USGS
 */
public class SpectrumExperiment extends Experiment {

  /**
   * True if plotting using Hz, False if sampleRate
   */
  private boolean freqSpace;

  private int[] respIndices;

  /**
   * Instantiates a noise experiment -- axis titles and scales
   */
  SpectrumExperiment() {
    super();
    respIndices = new int[3];
    freqSpace = false;
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

    XYSeriesCollection plotTimeseries = new XYSeriesCollection();
    plotTimeseries.setAutoWidth(true);

    int loadedDataCount = 0;
    for (int i = 0; i < 3; ++i) {
      if (dataStore.bothComponentsSet(i)) {
        ++loadedDataCount;
      }
    }
    respIndices = new int[loadedDataCount];

    // get the first (index.length) seed/resp pairs. while we expect to
    // have the first three plots be the ones with loaded data, in general
    // it is probably better to keep the program flexible against valid input
    for (int i = 0; i < respIndices.length; ++i) {
      // xth fully loaded function begins at 1
      int idx = dataStore.getXthFullyLoadedIndex(i + 1);
      respIndices[i] = idx;
      dataNames.add(dataStore.getBlock(idx).getName());
      dataNames.add(dataStore.getResponse(idx).getName());
    }

    fireStateChange("Getting PSDs of each series...");

    // gets the PSDs of each given index for given freqSpace
    for (int index : respIndices) {
      fireStateChange("Getting PSDs of data " + index + "...");
      addToPlot(dataStore, freqSpace, index, plotTimeseries);
    }

    plotTimeseries.addSeries(FFTResult.getLowNoiseModel(freqSpace));
    plotTimeseries.addSeries(FFTResult.getHighNoiseModel(freqSpace));

    xySeriesData.add(plotTimeseries);

  }

  @Override
  public int blocksNeeded() {
    return 3;
  }

  @Override
  public boolean hasEnoughData(DataStore dataStore) {
    for (int i = 0; i < 3; ++i) {
      if (dataStore.bothComponentsSet(i)) {
        return true;
      }
    }
    return false;
  }

  /**
   * NOTE: not used by corresponding panel, overrides with active indices
   * of components in the combo-box
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

  /**
   * Used for testing.
   * @return freqSpace
   */
  boolean getFreqSpace(){
    return freqSpace;
  }

}
