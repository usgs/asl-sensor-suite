package asl.sensor.experiment;

import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.event.EventListenerList;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.util.Pair;
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
   * @param idx Specifies which of the datablocks in the datastore to get the data from
   * @param xysc
   */
  public static void addToPlot(
      final DataStore ds,
      final boolean freqSpace,
      final int idx,
      XYSeriesCollection xysc) {

    XYSeries powerSeries =
        new XYSeries( "PSD " + ds.getBlock(idx).getName() + " [" + idx +"]" );

    Complex[] resultPSD = ds.getPSD(idx).getFFT();
    double[] freqs = ds.getPSD(idx).getFreqs();

    addToPlot(powerSeries, resultPSD, freqs, freqSpace, xysc);
  }

  /**
   *
   * Helper function to add data from a PSD calculation (Complex freq. space series and
   * corresponding frequencies taken from the FFTResult object produced by the calculation)
   * into an XYSeriesCollection to eventually be plotted.
   * Used in both self-noise and relative gain calculations
   * @param powerSeries XYSeries data to load given PSD calculations into
   * @param resultPSD FFT data produced by the PSD
   * @param freqs Frequency corresponding to PSD value at given index
   * @param freqSpace True if using units of Hz, False if units of s
   * @param xysc XYSeriesCollection the given XYSeries will be loaded into
   */
  public static void addToPlot(
      final XYSeries powerSeries,
      final Complex[] resultPSD,
      final double[] freqs,
      final boolean freqSpace,
      XYSeriesCollection xysc) {

    for (int j = 0; j < freqs.length; ++j) {
      if (1/freqs[j] > 1.0E3) {
        continue;
      }

      double temp = 10 * Math.log10( resultPSD[j].abs() );
      if (freqSpace) {
        powerSeries.add(freqs[j], temp);
      } else {
        powerSeries.add(1/freqs[j], temp);
      }
    }

    xysc.addSeries(powerSeries);

  }

  protected long start;
  protected long end;
  protected boolean statusToTerminal;
  protected List<XYSeriesCollection> xySeriesData;
  private String status;
  protected List<String> dataNames; // list of filenames of seed, resp files
  // NOTE: if implementing new experiment, best to use consistent ordering with
  // current set of experiments for this list:
  // SEED, RESP (if used), SEED, RESP (if used), etc.
  // That is, place response files after their associated timeseries
  protected Map<String, List<Pair<Date, Date>>> gapRegions;

  private EventListenerList eventHelper;

  /**
   * Initialize all fields common to experiment objects
   */
  public Experiment() {
    start = 0L; end = 0L;
    dataNames = new ArrayList<String>();
    status = "";
    eventHelper = new EventListenerList();
    statusToTerminal = false;
  }

  public void setTerminalPrintStatus(boolean print) {
    statusToTerminal = print;
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
   * Abstract function that runs the calculations specific to a given procedure
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
   * (Used to show messages displaying the progress of the function on the GUI)
   * @param newStatus Status change message to notify listeners of
   */
  protected void fireStateChange(String newStatus) {

    if (statusToTerminal) {
      System.out.println(newStatus);
    }

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
   * Get the end time of the data sent into this experiment, when timeseries data is used
   * @return End time, in microseconds
   */
  public long getEnd() {
    return end;
  }

  /**
   * Returns a map from datablock name identifiers to pairs of dates
   * representing regions within the given window where data does not exist
   * Gaps in data can show up in experiments as noise or other glitches and
   * may be evidence of a failing sensor in some cases
   * @return Map from data names to list of paired values representing gap durations (start, end)
   */
  public Map<String, List<Pair<Date, Date>>> getGapRegions() {
    return gapRegions;
  }

  /**
   * Get the names of data sent into program (set during backend calculations),
   * mainly used in report metadata generation
   * @return Names of data sent into the program (SNCLs, response filenames)
   */
  public List<String> getInputNames() {
    return dataNames;
  }

  /**
   * Get the start time of the data sent into this experiment, when timeseries data is used
   * @return Start time, in microseconds
   */
  public long getStart() {
    return start;
  }

  /**
   * Return newest status message produced by this program
   * (Used to get the actual status messages that should be displayed in the GUI while processing)
   * @return String representing status of program
   */
  public String getStatus() {
    return status;
  }

  /**
   * Used to check if the current input has enough data to do the calculation.
   * (i.e., if data requires 2 timeseries, check the first two datablocks in the datastore have
   * data and return true if so)
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
   * (Generally not used but may be useful for implementing forwarding status messages from data
   * subcomponents)
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

    xySeriesData = new ArrayList<XYSeriesCollection>();

    ds.matchIntervals( blocksNeeded() );

    gapRegions = new HashMap<String, List<Pair<Date, Date>>>();
    for (int i = 0; i < blocksNeeded(); ++i) {
      DataBlock block = ds.getBlock(i);
      String name = block.getName();
      // gaps already is calculated based on trimmed start and end times
      List<Pair<Long, Long>> gaps = block.getGapBoundaries();
      List<Pair<Date, Date>> gapsAsDates = new ArrayList<Pair<Date, Date>>();
      for (Pair<Long, Long> gap : gaps) {
        Date start = new Date( gap.getFirst() );
        Date end = new Date( gap.getSecond() );
        gapsAsDates.add( new Pair<Date, Date>(start, end) );
      }
      gapRegions.put(name, gapsAsDates);
    }

    fireStateChange("Beginning calculations...");

    backend(ds);

    fireStateChange("Calculations done!");
  }

}
