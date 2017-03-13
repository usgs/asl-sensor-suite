package asl.sensor.experiment;

import org.apache.commons.math3.complex.Complex;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import asl.sensor.input.DataBlock;
import asl.sensor.input.DataStore;

/**
 * This function defines template patterns for each type of sensor experiment
 * (we use this term in the code to prevent confusion with the overloaded term
 * "test"). Concrete extensions of this class are used to define a backend
 * for the calculations of the experiment to be passed into a class that does
 * plotting for the data that results from this object.
 * 
 * Some experiment implementations may not only produce XY series data to be
 * plotted by the corresponding GUI implemntation, but also produce additional
 * statistical data relevant to the plots, such as residual calculations
 * or the values of best-fit parameters given a series of inputs.
 * 
 * @author akearns
 *
 */
public abstract class Experiment {
  
  // defines template pattern for each type of test, given by backend
  // each test returns new (set of) timeseries data from the input data
  
  // TODO: include axis scale definitions constructed along with axis titles
  // (i.e., logarithmic, etc.)
  
  /**
   * Helper function to add data from a datastore object (the PSD calculation)
   * into an XYSeriesCollection to eventually be plotted
   * Used in both self-noise and relative gain calculations
   * @param ds DataStore to collect data from
   * @param freqSpace True if using units of Hz, False if units of s
   * (sample rate vs. interval between points)
   * @param indices Specifies which of the data can be loaded in
   * @param xysc
   */
  public static void addToPlot(
      final DataStore ds,
      final boolean freqSpace,
      final int[] indices,
      XYSeriesCollection xysc) {
    
    // DataBlock[] dataIn = ds.getData();
    
    for (int i = 0; i < indices.length; ++i) {
      
      int idx = indices[i];
      
      XYSeries powerSeries = 
          new XYSeries( "PSD " + ds.getBlock(idx).getName() + " [" + idx +"]" );

      Complex[] resultPSD = ds.getPSD(idx).getFFT();
      double[] freqs = ds.getPSD(idx).getFreqs();

      for (int j = 0; j < freqs.length; ++j) {
        if (1/freqs[j] > 1.0E3) {
          continue;
        }

        // TODO: is this right (seems to be)
        Complex temp = resultPSD[j].multiply(Math.pow(2*Math.PI*freqs[j],4));

        if (freqSpace) {
          powerSeries.add( freqs[j], 10*Math.log10( temp.abs() ) );
        } else {
          powerSeries.add( 1/freqs[j], 10*Math.log10( temp.abs() ) );
        }
      }

      xysc.addSeries(powerSeries);

    }
    
  }
  
  protected XYDataset xySeriesData;
  
  /**
   * Abstract function that
   * (Overwritten by concrete experiments with specific operations)
   * @param ds Object containing the raw timeseries data to process
   * @param freqSpace True if the data should be plotted by frequency (False
   * if it should be plotted by period) for domain axis
   */
  protected abstract void backend(
      final DataStore ds, final boolean freqSpace);
  
  /**
   * Return the plottable data for this experiment, populated in the backend
   * function of an implementing class; calling this class before running the
   * setData function / backend will produce initialization errors (NPE).
   * @return Plottable data 
   */
  public XYDataset getData() {
    return xySeriesData;
  }
  
  /**
   * Driver to do data processing on inputted data (calls a concrete backend
   * method which is different for each type of experiment)
   * This function specifically (rather than the backend implementation) is
   * where interval consistency is checked before doing calculations.
   * @param ds Timeseries data to be processed
   * @param freqSpace True if the x-axis should be frequency units (False if it
   * should be units of time, for the period)
   */
  public void setData(final DataStore ds, final boolean freqSpace) {
    
    if ( hasEnoughData(ds) && ds.numberOfBlocksSet() == 0) {
      backend(ds, freqSpace);
      return;
    }
    
    final DataBlock db = ds.getXthLoadedBlock(1);

    long interval = db.getInterval();
    
    final DataBlock[] dataIn = ds.getData();
    
    // most backend functions will use an xyseries collection
    // (this can, of course, be overridden for the few that won't)
    xySeriesData = new XYSeriesCollection();
    
    // int length = dataIn[0].size();
    for (final DataBlock data : dataIn) {
      
      if ( data == null) {
        // we can have null blocks, but can't get interval from a null block
        continue;
      }
      
      if ( data.getInterval() != interval ) {
        System.out.println( interval+","+data.getInterval() );
        ds.matchIntervals();
      }
    }
    
    backend(ds, freqSpace);
  }
  
  public abstract boolean hasEnoughData(final DataStore ds);
  
}
