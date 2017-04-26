package asl.sensor.experiment;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.complex.Complex;
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

        // unit conversion with response
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
  
  long start;
  long end;
  protected List<XYSeriesCollection> xySeriesData;
  
  public Experiment() {
    start = 0L; end = 0L;
  }
  
  /**
   * Abstract function that
   * (Overwritten by concrete experiments with specific operations)
   * @param ds Object containing the raw timeseries data to process
   */
  protected abstract void backend(final DataStore ds);
  
  /**
   * Return the number of data blocks needed by the experiment
   * (Used in determining the number of input plots needed to be shown)
   * @return Number of blocks needed as integer
   */
  public abstract int blocksNeeded();
  
  /**
   * Return the plottable data for this experiment, populated in the backend
   * function of an implementing class; calling this class before running the
   * setData function / backend will produce initialization errors (NPE).
   * The results are returned as a list, where each list is the data to be
   * placed into a separate chart.
   * @return Plottable data 
   */
  public List<XYSeriesCollection> getData() {
    return xySeriesData;
  }
  
  /**
   * Get the end time of the data sent into this experiment
   * @return End time, in microseconds
   */
  public long getEnd() {
    return end;
  }
  
  /**
   * Get the start time of the data sent into this experiment
   * @return Start time, in microseconds
   */
  public long getStart() {
    return start;
  }
  
  /**
   * Used to check if the current input has enough data to do the calculation
   * @param ds DataStore to be fed into experiment calculation
   * @return True if there is enough data to be run
   */
  public abstract boolean hasEnoughData(final DataStore ds);
  
  /**
   * Driver to do data processing on inputted data (calls a concrete backend
   * method which is different for each type of experiment)
   * This function specifically (rather than the backend implementation) is
   * where interval consistency is checked before doing calculations.
   * @param ds Timeseries data to be processed
   */
  public void setData(final DataStore ds) {
    
    if ( hasEnoughData(ds) && ( blocksNeeded() == 0 ) ) {
      // prevent null issue 
      xySeriesData = new ArrayList<XYSeriesCollection>();
      start = 0L;
      end = 0L;
      backend(ds);
      return;
    }
    
    final DataBlock db = ds.getXthLoadedBlock(1);
    
    start = db.getStartTime();
    end = db.getEndTime();

    long interval = db.getInterval();
    
    final DataBlock[] dataIn = ds.getData();
    
    xySeriesData = new ArrayList<XYSeriesCollection>();
    
    // int length = dataIn[0].size();
    for (final DataBlock data : dataIn) {
      
      if ( data == null) {
        // we can have null blocks, but can't get interval from a null block
        continue;
      }
      
      if ( data.getInterval() != interval ) {
        // System.out.println( interval+","+data.getInterval() );
        ds.matchIntervals();
      }
    }
    
    backend(ds);
  }
  
  /**
   * Return an array of indices of responses used by an index, to include
   * data in report generation
   * @return Indices in which responses are required
   */
  public int[] listActiveResponseIndices() {
    // override this in functions that use a backend including responses
    return new int[]{};
  }
   
}
