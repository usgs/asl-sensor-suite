package asl.sensor.experiment;

import org.jfree.data.xy.XYSeriesCollection;

import asl.sensor.input.DataStore;
import asl.sensor.input.InstrumentResponse;
import asl.sensor.utils.FFTResult;

/**
 * Calculates PSD to get cross-power.
 * Based on code in the seedscan timeseries package, see
 * https://github.com/usgs/seedscan/tree/master/src/main/java/asl/timeseries
 * @author akearns, jholland 
 *
 */
public class SpectrumExperiment extends Experiment {
  

  protected boolean freqSpace;
  
  protected int[] respIndices;
  
  /**
   * Instantiates a noise experiment -- axis titles and scales
   */
  public SpectrumExperiment() {
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
  protected void backend(final DataStore ds) {
    
    XYSeriesCollection xysc = new XYSeriesCollection();
    xysc.setAutoWidth(true);
    
    int loadedDataCount = 0;
    for (int i = 0; i < 3; ++i) {
      if ( ds.bothComponentsSet(i) ) {
        ++loadedDataCount;
      }
    }
    respIndices = new int[loadedDataCount];
    
    // get the first (index.length) seed/resp pairs. while we expect to
    // have the first three plots be the ones with loaded data, in general
    // it is probably better to keep the program flexible against valid input
    for (int i = 0; i < respIndices.length; ++i) {
      // xth fully loaded function begins at 1
      int idx = ds.getXthFullyLoadedIndex(i+1);
      respIndices[i] = idx;
      dataNames.add( ds.getBlock(idx).getName() );
      dataNames.add( ds.getResponse(idx).getName() );
    }
    
    InstrumentResponse[] responses = new InstrumentResponse[respIndices.length];
    
    for (int i = 0; i < respIndices.length; ++i) {
      responses[i] = ds.getResponse(respIndices[i]);
    }
    
    fireStateChange("Getting PSDs of each series...");
    
    // gets the PSDs of each given index for given freqSpace
    for (int i = 0; i < respIndices.length; ++i) {
      int idx = respIndices[i];
      fireStateChange("Getting PSDs of data " + idx + "...");
      addToPlot(ds, freqSpace, idx, xysc);
    }
    
    xysc.addSeries( FFTResult.getLowNoiseModel(freqSpace) );
    xysc.addSeries( FFTResult.getHighNoiseModel(freqSpace) );
    
    xySeriesData.add(xysc);

  }

  @Override
  public int blocksNeeded() {
    return 3;
  }
  
  @Override
  public boolean hasEnoughData(DataStore ds) {
    for (int i = 0; i < 3; ++i) {
      if ( ds.bothComponentsSet(i) ) {
        return true;
      }
    }
    return false;
  }

  @Override
  public int[] listActiveResponseIndices() {
    // NOTE: not used by corresponding panel, overrides with active indices
    // of components in the combo-box
    return respIndices;
  }

  /**
   * Used to set the x-axis over which the PSDs / cross-powers are plotted,
   * either frequency (Hz) units or sample-interval (s) units
   * @param freqSpace True if the plot should use units of Hz
   */
  public void setFreqSpace(boolean freqSpace) {
    this.freqSpace = freqSpace;
  }
  
}
