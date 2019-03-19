package asl.sensor.experiment;

import asl.sensor.utils.TimeSeriesUtils;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.time.Instant;
import java.util.ArrayList;
import java.util.Date;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.TimeZone;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.event.EventListenerList;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.util.Pair;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import asl.sensor.input.DataBlock;
import asl.sensor.input.DataStore;
import asl.sensor.utils.NumericUtils;

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
 * plotted by the corresponding GUI implementation, but also produce additional
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
 * @author akearns - KBRWyle
 */
public abstract class Experiment {

  // defines template pattern for each type of test, given by backend
  // each test returns new (set of) timeseries data from the input data

  public static final ThreadLocal<DecimalFormat> DECIMAL_FORMAT =
      ThreadLocal.withInitial(() -> {
        DecimalFormat format = new DecimalFormat("#.###");
        NumericUtils.setInfinityPrintable(format);
        return format;
      });
  static final ThreadLocal<SimpleDateFormat> DATE_FORMAT =
      ThreadLocal.withInitial(() -> {
        SimpleDateFormat format = new SimpleDateFormat("YYYY.DDD");
        format.setTimeZone(TimeZone.getTimeZone("UTC"));
        return format;
      });

  /**
   * Frequency plots should be limited by this value as max resolution of period (1E6 seconds)
   */
  static final double MAX_PLOT_PERIOD = 1.0E6;
  private final EventListenerList eventHelper;
  long start;
  long end;
  List<XYSeriesCollection> xySeriesData;
  /**
   * list of filenames of seed, resp files
   * NOTE: if implementing new experiment, best to use consistent ordering with
   * current set of experiments for this list:
   * SEED, RESP (if used), SEED, RESP (if used), etc.
   * That is, place response files after their associated timeseries
   */
  List<String> dataNames;
  private String status;
  private Map<String, List<Pair<Date, Date>>> gapRegions;
  /**
   * Initialize all fields common to experiment objects
   */
  Experiment() {
    start = 0L;
    end = 0L;
    dataNames = new ArrayList<>();
    status = "";
    eventHelper = new EventListenerList();
  }

  /**
   * Helper function to add data from a DataStore object (the PSD calculation)
   * into an XYSeriesCollection to eventually be plotted
   * Used in both self-noise and relative gain calculations
   *
   * @param dataStore DataStore to collect data from
   * @param freqSpace True if using units of Hz, False if units of s
   * (sample rate vs. interval between points)
   * @param index Specifies which of the DataBlocks in the DataStore to get the data from
   */
  static void addToPlot(
      final DataStore dataStore,
      final boolean freqSpace,
      final int index,
      XYSeriesCollection xysc) {

    XYSeries powerSeries =
        new XYSeries("PSD " + dataStore.getBlock(index).getName() + " [" + index + "]");

    Complex[] resultPSD = dataStore.getPSD(index).getFFT();
    double[] freqs = dataStore.getPSD(index).getFreqs();

    addToPlot(powerSeries, resultPSD, freqs, freqSpace, xysc);
  }

  /**
   * Stub method to be overridden for other methods to produce String data for experiment result.
   * Includes formatting of numeric data. This may not be used for all experiments.
   * @return String containing human-readable data
   */
  String[] getDataStrings() {
    return new String[]{""};
  }

  /**
   * Stub method to be overridden for other methods to produce String data for plot data.
   * Includes formatting of numeric data and is usually designed to include dates of data that
   * are plotted in frequency space (i.e., Hz or period seconds).
   * This may not be used for all plots.
   * @return String containing human-readable data
   */
  public String[] getInsetStrings() {
    return getDataStrings();
  }

  /**
   * Stub method to be overridden for other methods to produce String data for reports.
   * Includes formatting of numeric data. This may not be used for all plots.
   * @return String containing human-readable data
   */
  public String getReportString() {
    StringBuilder sb = new StringBuilder();
    String[] strings = getDataStrings();
    for (int i = 0; i < strings.length; ++i) {
      String insetString = strings[i];
      sb.append(insetString);
      // add space between inset strings
      if (i + 1 < strings.length) {
        sb.append('\n');
      }
    }
    return sb.toString();
  }

  /**
   * Method to get a formatted string with start and end dates of data, to be used in
   * producing reports of the given data. This is empty if start and end are set to the same value
   * (i.e., the experiment does not produce a result on timeseries data, such as RESP plots)
   * @return String of formatted start and end, if they have been set
   */
  public String getFormattedDateRange() {
    StringBuilder sb = new StringBuilder();
    if (start != end) {
      sb.append(getFormattedStartDate());
      sb.append('\n');
      sb.append(getFormattedEndDate());
    }
    return sb.toString();
  }

  /**
   * Return a formatted string with the given data's start time
   * @return formatted start time (Julian day)
   */
  String getFormattedStartDate() {
    return "Data start time:\n" + TimeSeriesUtils.formatEpochMillis(start);
  }

  String getFormattedEndDate() {
    return "Data end time:\n" + TimeSeriesUtils.formatEpochMillis(end) + '\n';
  }

  /**
   * Helper function to add data from a PSD calculation (Complex freq. space series and
   * corresponding frequencies taken from the FFTResult object produced by the calculation)
   * into an XYSeriesCollection to eventually be plotted.
   * Used in both self-noise and relative gain calculations
   *
   * @param powerSeries XYSeries data to load given PSD calculations into
   * @param resultPSD FFT data produced by the PSD
   * @param freqs Frequency corresponding to PSD value at given index
   * @param freqSpace True if using units of Hz, False if units of s
   * @param xysc XYSeriesCollection the given XYSeries will be loaded into
   */
  static void addToPlot(
      final XYSeries powerSeries,
      final Complex[] resultPSD,
      final double[] freqs,
      final boolean freqSpace,
      XYSeriesCollection xysc) {

    // Smooth the PSD data before it goes out to the plots
    Complex[] smoothedPSD = NumericUtils.multipointMovingAverage(resultPSD, 9, false);
    // for the last 3 points, do 7, 5, 3 last points
    Complex last3 = Complex.ZERO;
    Complex last5 = Complex.ZERO;
    Complex last7 = Complex.ZERO;
    for (int i = 0; i < 7; ++i) {
      if (i < 3) {
        last3 = last3.add(resultPSD[i]);
      }
      if (i < 5) {
        last5 = last5.add(resultPSD[i]);
      }
      last7 = last7.add(resultPSD[i]);
    }

    int idx = smoothedPSD.length - 1;
    smoothedPSD[idx] = last3.divide(3);
    --idx;
    smoothedPSD[idx] = last5.divide(5);
    --idx;
    smoothedPSD[idx] = last7.divide(7);

    for (int j = 0; j < freqs.length; ++j) {
      if (1 / freqs[j] > MAX_PLOT_PERIOD) {
        continue;
      }
      double temp = 10 * Math.log10(smoothedPSD[j].abs());
      if (freqSpace) {
        powerSeries.add(freqs[j], temp);
      } else {
        powerSeries.add(1 / freqs[j], temp);
      }
    }

    xysc.addSeries(powerSeries);
  }

  /**
   * Add an object to the list of objects to be notified when the experiment's
   * status changes
   *
   * @param listener ChangeListener to be notified (i.e., parent panel)
   */
  public void addChangeListener(ChangeListener listener) {
    eventHelper.add(ChangeListener.class, listener);
  }

  /**
   * Abstract function that runs the calculations specific to a given procedure,
   * overwritten by concrete experiments with specific operations.
   * Information on what an experiment's implementation does is in its documentation intro section.
   *
   * @param dataStore Object containing the raw timeseries data to process
   */
  protected abstract void backend(final DataStore dataStore);

  /**
   * Return the number of data blocks needed by the experiment
   * (Used in determining the number of input plots needed to be shown)
   *
   * @return Number of blocks needed as integer
   */
  public abstract int blocksNeeded();

  /**
   * Update processing status and notify listeners of change
   * (Used to show messages displaying the progress of the function on the GUI)
   *
   * @param newStatus Status change message to notify listeners of
   */
  void fireStateChange(String newStatus) {
    status = newStatus;
    ChangeListener[] listeners = eventHelper.getListeners(ChangeListener.class);
    if (listeners != null && listeners.length > 0) {
      ChangeEvent event = new ChangeEvent(this);
      for (ChangeListener listener : listeners) {
        listener.stateChanged(event);
      }
    }
  }

  /**
   * Return the plottable data for this experiment, populated in the backend
   * function of an implementing class; calling this class before running the
   * setData function / backend will produce initialization errors (NPE).
   * The results are returned as a list, where each list is the data to be
   * placed into a separate chart.
   *
   * @return Plottable data
   */
  public List<XYSeriesCollection> getData() {
    return xySeriesData;
  }

  /**
   * Get the end time of the data sent into this experiment, when timeseries data is used
   *
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
   *
   * @return Map from data names to list of paired values representing gap durations (start, end)
   */
  public Map<String, List<Pair<Date, Date>>> getGapRegions() {
    return gapRegions;
  }

  /**
   * Get the names of data sent into program (set during backend calculations),
   * mainly used in report metadata generation
   *
   * @return Names of data sent into the program (SNCLs, response filenames)
   */
  public List<String> getInputNames() {
    return dataNames;
  }

  /**
   * Get the start time of the data sent into this experiment, when timeseries data is used
   *
   * @return Start time, in microseconds
   */
  public long getStart() {
    return start;
  }

  /**
   * Return newest status message produced by this program
   * (Used to get the actual status messages that should be displayed in the GUI while processing)
   *
   * @return String representing status of program
   */
  public String getStatus() {
    return status;
  }

  /**
   * Used to check if the current input has enough data to do the calculation.
   * (i.e., if data requires 2 timeseries, check the first two datablocks in the datastore have
   * data and return true if so)
   *
   * @param dataStore DataStore to be fed into experiment calculation
   * @return True if there is enough data to be run
   */
  public abstract boolean hasEnoughData(final DataStore dataStore);

  /**
   * Return an array of indices of responses used by an index, to include
   * data in report generation
   *
   * @return Indices in which responses are required
   */
  public int[] listActiveResponseIndices() {
    // override this in functions that use a backend including responses
    return new int[]{};
  }

  /**
   * Driver to do data processing on inputted data (calls a concrete backend
   * method which is different for each type of experiment)
   * This function specifically (rather than the backend implementation) is
   * where interval consistency is checked before doing calculations.
   *
   * @param dataStore Timeseries data to be processed
   */
  public void runExperimentOnData(final DataStore dataStore) {

    fireStateChange("Beginning loading data...");

    dataNames = new ArrayList<>();
    xySeriesData = new ArrayList<>();
    gapRegions = new LinkedHashMap<>();

    if (hasEnoughData(dataStore) && (blocksNeeded() == 0)) {
      // prevent null issue when doing response data, which does not really have times
      start = 0L;
      end = 0L;
      backend(dataStore);
      return;
    }

    final DataBlock db = dataStore.getXthLoadedBlock(1);

    start = db.getStartTime();
    end = db.getEndTime();

    dataStore.matchIntervals(blocksNeeded());

    // populate gapregions data
    for (int i = 0; i < blocksNeeded(); ++i) {
      // in the case of spectrum panel, not all data inputs may be set
      // lockout check to make sure enough needed data is set already done, so this is OK
      if (!dataStore.blockIsSet(i)) {
        continue;
      }
      DataBlock block = dataStore.getBlock(i);
      String name = block.getName();
      // gaps already is calculated based on trimmed start and end times
      List<Pair<Long, Long>> gaps = block.getGapBoundaries();
      List<Pair<Date, Date>> gapsAsDates = new ArrayList<>();
      for (Pair<Long, Long> gap : gaps) {
        Date start = new Date(gap.getFirst());
        Date end = new Date(gap.getSecond());
        gapsAsDates.add(new Pair<>(start, end));
      }
      gapRegions.put(name, gapsAsDates);
    }

    fireStateChange("Beginning calculations...");

    backend(dataStore);

    fireStateChange("Calculations done!");
  }
}
