package asl.sensor.experiment;

import java.util.ArrayList;
import java.util.List;

import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.event.EventListenerList;

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
 * Experiments work in a manner similar to builder patterns: experiments that
 * rely on variables to determine how their calculations are run, such as
 * the randomized experiment using a boolean to determine what frequency poles
 * should be set or the noise experiment using a boolean to determine if it
 * should return results in frequency or period space for the x-axis. Set these
 * values first, and then call "runExperimentOnData" with a given DataStore
 * containing the relevant values.
 * 
 * Some experiment implementations may not only produce XY series data to be
 * plotted by the corresponding GUI implemntation, but also produce additional
 * statistical data relevant to the plots, such as residual calculations
 * or the values of best-fit parameters given a series of inputs. These should
 * not be called unless the experiment has already been run, as they will
 * otherwise not be populated with valid results. Because the primary use case
 * of an experiment is to be run by the GUI panel containing it, which will
 * read in all the additional results as soon as the calculation completes,
 * there is no particular safeguard against doing so at this time. The GUI
 * panels should read in data during their updateData routine, which is where
 * the experiment should be run.
 * 
 * @author akearns
 *
 */
public abstract class Experiment {
  
  // defines template pattern for each type of test, given by backend
  // each test returns new (set of) timeseries data from the input data
  
  public static final String STATUS = "status";
  
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
  private String status;
  protected List<String> dataNames; // list of filenames of seed, resp files
  // NOTE: if implementing new experiment, best to use consistent ordering with
  // current set of experiments for this list: 
  // SEED, RESP (if used), SEED, RESP (if used), etc.
  // That is, place response files after their associated timeseries
  
  private EventListenerList eventHelper;
  
  public Experiment() {
    start = 0L; end = 0L;
    dataNames = new ArrayList<String>();
    status = "";
    eventHelper = new EventListenerList();
  }
  
  /**
   * Add an object to the list of objects to be notified when the experiment's
   * status changes
   * @param listener ChangeListener to be notified (i.e., parent panel)
   */
  public void addChangeListener(ChangeListener listener) {
     eventHelper.add(ChangeListener.class, listener);
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
   * Update processing status and notify listeners of change
   * @param newStatus Status change message to notify listeners of
   */
  protected void fireStateChange(String newStatus) {
    status = newStatus;
    ChangeListener[] lsners = eventHelper.getListeners(ChangeListener.class);
    if (lsners != null && lsners.length > 0) {
      ChangeEvent evt = new ChangeEvent(this);
      for (ChangeListener lsnr : lsners) {
        lsnr.stateChanged(evt);
      }
    }
  }
  
  /**
   * Return newest status message produced by this program
   * @return String representing status of program
   */
  public String getStatus() {
    return status;
  }
  
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
   * Get the names of data sent into program (set during backend calculations)
   * @return
   */
  public List<String> getInputNames() {
    return dataNames;
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
   * Return an array of indices of responses used by an index, to include
   * data in report generation
   * @return Indices in which responses are required
   */
  public int[] listActiveResponseIndices() {
    // override this in functions that use a backend including responses
    return new int[]{};
  }
   
  /**
   * Remove changelistener from list of listeners notified on status change
   * @param listener Listener to remove from list
   */
  public void removeChangeListener(ChangeListener listener) {
      eventHelper.remove(ChangeListener.class, listener);
  }
  
  /**
   * Driver to do data processing on inputted data (calls a concrete backend
   * method which is different for each type of experiment)
   * This function specifically (rather than the backend implementation) is
   * where interval consistency is checked before doing calculations.
   * @param ds Timeseries data to be processed
   */
  public void runExperimentOnData(final DataStore ds) {
    
    status = "";
    
    fireStateChange("Beginning loading data...");
    
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

    dataNames = new ArrayList<String>();
    
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
        fireStateChange("Downsampling data...");
        // System.out.println( interval+","+data.getInterval() );
        ds.matchIntervals();
      }
    }
    
    fireStateChange("Beginning calculations...");
    
    backend(ds);
    
    fireStateChange("Calculations done!");
  }
   
}
