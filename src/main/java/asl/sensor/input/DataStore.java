package asl.sensor.input;

import static asl.utils.TimeSeriesUtils.ONE_HZ_INTERVAL;
import static asl.utils.TimeSeriesUtils.formatEpochMillis;
import static asl.utils.TimeSeriesUtils.getMplexNameList;
import static asl.utils.TimeSeriesUtils.getTimeSeries;

import asl.utils.FFTResult;
import asl.utils.input.DataBlock;
import asl.utils.input.InstrumentResponse;
import edu.iris.dmc.seedcodec.CodecException;
import edu.sc.seis.seisFile.mseed.SeedFormatException;
import java.io.IOException;
import java.time.Instant;
import org.apache.commons.math3.util.Pair;

/**
 * Holds the inputted data from miniSEED files both as a simple struct
 * (see DataBlock) and as a plottable format for the DataPanel.
 * Also serves as container for loaded-in instrument response files.
 *
 * This structure is build around a series of arrays that hold the data loaded
 * in from a miniseed file, from a response file, from PSD calculations
 * (power spectral density) using the miniseed and response data, and for
 * plotting in a separate class. Each of these arrays is of the same length
 * (see the FILE_COUNT parameter).
 *
 * Each subcomponent of this class is designed to match up by index. That is,
 * for index i, the DataBlock at index i will be associated with a response
 * at i. Both of these will be used to calculate the power-spectral
 * density at i.
 *
 * This structure also includes functions to get the Xth (i.e., first, second)
 * set of valid data. While other code in this program expects data to be
 * loaded in sequentially from index 0 to whatever the maximum needed input is,
 * this allows the datastore to remain somewhat flexible about the means in
 * which the data it lays out is stored.
 *
 * This class also has means with which to check that data is properly trimmed
 * to the same range and has the same sample rate, which is necessary for most
 * experiments.
 *
 * @author akearns
 */
public class DataStore {

  /**
   * Defines the maximum number of plots to be shown
   */
  public final static int FILE_COUNT = 9;
  private final DataBlock[] dataBlockArray;
  private final InstrumentResponse[] responses;

  // these are used to check to make sure data has been loaded
  private final boolean[] thisBlockIsSet;
  private final boolean[] thisResponseIsSet;

  /**
   * Instantiate the collections, including empty datasets to be sent to
   * charts for plotting (see DataPanel)
   */
  public DataStore() {
    dataBlockArray = new DataBlock[FILE_COUNT];
    responses = new InstrumentResponse[FILE_COUNT];
    thisBlockIsSet = new boolean[FILE_COUNT];
    thisResponseIsSet = new boolean[FILE_COUNT];
    for (int i = 0; i < FILE_COUNT; ++i) {
      thisBlockIsSet[i] = false;
      thisResponseIsSet[i] = false;
    }
  }


  /**
   * Create a copy of the current datastore
   *
   * @param ds datastore to copy
   */
  public DataStore(DataStore ds) {
    dataBlockArray = new DataBlock[FILE_COUNT];
    responses = new InstrumentResponse[FILE_COUNT];
    thisBlockIsSet = new boolean[FILE_COUNT];
    thisResponseIsSet = new boolean[FILE_COUNT];
    boolean[] setBlocks = ds.dataIsSet();
    boolean[] setResps = ds.responsesAreSet();
    for (int i = 0; i < FILE_COUNT; ++i) {
      if (setBlocks[i]) {
        dataBlockArray[i] = new DataBlock(ds.getBlock(i));
        thisBlockIsSet[i] = true;
      }

      if (setResps[i]) {
        responses[i] = ds.getResponse(i);
        thisResponseIsSet[i] = true;
      }
    }
  }


  /**
   * Determine if any data blocks in the data store has been initialized
   *
   * @return true if at least one datablock has been added to the store object
   */
  public boolean areAnyBlocksSet() {

    return areAnyBlocksSet(FILE_COUNT);
  }

  /**
   * Determine if any data blocks in the current active plotting set have been initialized
   * @param activePlots limit of plots shown by GUI currently
   * @return true if at least one datablock visible in GUI has had data added
   */
  public boolean areAnyBlocksSet(int activePlots) {

    activePlots = Math.min(activePlots, FILE_COUNT);
    for (int i = 0; i < activePlots; ++i) {
      if (blockIsSet(i)) {
        return true;
      }
    }
    return false;
  }

  /**
   * Checks whether the data block (timeseries) at a given index
   * contains data or is null/empty
   *
   * @param idx Index of the data to look at (between 0 and DATA_STORE)
   * @return True if a seed file has been loaded in at the index
   */
  public boolean blockIsSet(int idx) {
    return thisBlockIsSet[idx];
  }

  /**
   * Checks if both components at a specific index are set
   *
   * @param idx Index of data to check if set or not
   * @return True if a seed and response have been both loaded in
   */
  public boolean bothComponentsSet(int idx) {
    return (thisBlockIsSet[idx] && thisResponseIsSet[idx]);
  }

  /**
   * Returns the boolean array where each index shows if there is a DataBlock
   * loaded into this datastore object at that index
   * has had data loaded into it (see blockIsSet method)
   *
   * @return Array of booleans where an entry is true if data is loaded at
   * the corresponding index
   */
  private boolean[] dataIsSet() {
    return thisBlockIsSet;
  }

  /**
   * Return a single data block according to the passed index
   *
   * @param idx Index of datablock, corresponding to data panel plot index
   * @return Timeseries data for corresponing plot
   */
  public DataBlock getBlock(int idx) {
    return dataBlockArray[idx];
  }

  public Pair<Long, Long> getCommonTime() {
    return getCommonTime(FILE_COUNT);
  }

  public Pair<Long, Long> getCommonTime(int limit) {
    if (numberOfBlocksSet() < 1) {
      return new Pair<>(Long.MIN_VALUE, Long.MAX_VALUE);
    } else {
      long lastStartTime = Long.MIN_VALUE;
      long firstEndTime = Long.MAX_VALUE;

      // first pass to get the limits of the time data
      for (int i = 0; i < limit; ++i) {
        DataBlock data = dataBlockArray[i];
        if (!thisBlockIsSet[i]) {
          continue;
        }
        long start = data.getInitialStartTime();
        if (start > lastStartTime) {
          lastStartTime = start;
        }
        long end = data.getInitialEndTime();
        if (end < firstEndTime) {
          firstEndTime = end;
        }
      }
      return new Pair<>(lastStartTime, firstEndTime);
    }
  }

  /**
   * Gets the power-spectral density of an index in this object.
   * If a PSD has already been calculated, this will return that. If not,
   * it will calculate the result, store it, and then return that data.
   *
   * @param idx Index of data to get the PSD of
   * @return Complex array of frequency values and a
   * double array of the frequencies
   */
  public FFTResult getPSD(int idx) {
    double[] data = dataBlockArray[idx].getData();
    long interval = dataBlockArray[idx].getInterval();
    InstrumentResponse ir = responses[idx];
    return FFTResult.crossPower(data, data, ir, ir, data.length, interval);
  }

  /**
   * Gets the power-spectral density of an index in this object.
   * If a PSD has already been calculated, this will return that. If not,
   * it will calculate the result, store it, and then return that data.
   *
   * This version of the function is meant to be used in cases where exact length
   * of the PSD must be specified in advance, such as when working with data that
   * may have timing differences due to quantization that mean trimming makes
   * some data be of different lengths from what is expected.
   *
   * @param idx Index of data to get the PSD of
   * @param maxLength Maximum number of points to calculate PSD over -- range 0 to maxLength
   * @return Complex array of frequency values and a
   * double array of the frequencies
   */
  public FFTResult getPSD(int idx, int maxLength) {
    double[] data = dataBlockArray[idx].getData();
    long interval = dataBlockArray[idx].getInterval();
    InstrumentResponse ir = responses[idx];
    return FFTResult.crossPower(data, data, ir, ir, maxLength, interval);
  }

  /**
   * Get the instrument response object at a given index
   *
   * @param idx Index to get the response for
   * @return The instrument response data (gain, poles, zeros, etc.)
   */
  public InstrumentResponse getResponse(int idx) {
    return responses[idx];
  }

  /**
   * Used to get the first, second, etc. data set loaded. Used when operations
   * reading in data don't require all the inputs to be loaded.
   * Requires both SEED and RESP to be loaded for this to be valid.
   *
   * @param x x-th set of loaded data to get, starting at 1 (NOT 0)
   * @return index of the loaded data
   */
  public int getXthFullyLoadedIndex(int x) {
    if (x < 1) {
      throw new IndexOutOfBoundsException("Parameter must be >= 1");
    }

    int loaded = 0;
    for (int i = 0; i < FILE_COUNT; ++i) {
      if (bothComponentsSet(i)) {
        ++loaded;
        if (loaded == x) {
          return i;
        }
      }
    }

    String errMsg = "Not enough data loaded in (found " + loaded + ")";
    throw new IndexOutOfBoundsException(errMsg);
  }

  /**
   * Used to get the first, second, etc. loaded block, whether or not it has
   * a loaded response file as well.
   * Used to find the panel where a step calibration is loaded
   *
   * @param x x-th set of data to get, starting at 1 (NOT 0)
   * @return The Xth DataBlock in this object that is not null
   */
  public DataBlock getXthLoadedBlock(int x) {
    if (x < 1) {
      throw new IndexOutOfBoundsException("Parameter must be >= 1");
    }

    int count = 0;
    for (int i = 0; i < FILE_COUNT; ++i) {
      if (thisBlockIsSet[i]) {
        ++count;
        if (count == x) {
          return dataBlockArray[i];
        }
      }
    }

    String errMsg = "Not enough data loaded in (found " + count + ")";
    throw new IndexOutOfBoundsException(errMsg);
  }

  /**
   * Checks if there is any data at all loaded into this object so far,
   * either data or response
   *
   * @return True if there is anything loaded into this datastore object.
   */
  public boolean isAnythingSet() {
    for (int i = 0; i < FILE_COUNT; ++i) {
      if (thisBlockIsSet[i] || thisResponseIsSet[i]) {
        return true;
      }
    }
    return false;
  }

  /**
   * Get lowest-frequency data and downsample all data to it
   */
  public void matchIntervals() {
    matchIntervals(FILE_COUNT);
  }

  /**
   * Math the first [limit] inputs' intervals to the lowest-frequency used by
   * any of the blocks within that range
   *
   * @param limit Index of last block to match intervals to
   */
  public void matchIntervals(int limit) {
    long interval = 0;
    // first loop to get lowest-frequency data
    for (int i = 0; i < limit; ++i) {
      if (thisBlockIsSet[i]) {
        interval = Math.max(interval, getBlock(i).getInterval());
      }
    }
    // second loop to downsample
    for (int i = 0; i < limit; ++i) {
      if (thisBlockIsSet[i] && getBlock(i).getInterval() != interval) {
        getBlock(i).resample(interval);
      }
    }

    trimToCommonTime();
  }

  /**
   * Gives the number of DataBlocks (miniseed files) read in to this object
   *
   * @return number of files read in
   */
  public int numberOfBlocksSet() {
    int loaded = 0;
    for (int i = 0; i < FILE_COUNT; ++i) {
      if (thisBlockIsSet[i]) {
        ++loaded;
      }
    }
    return loaded;
  }


  /**
   * Removes all data at a specific index -- miniseed, response, and any
   * data generated from them
   */
  public void removeData(int idx) {
    removeBlock(idx);
    responses[idx] = null;
    thisResponseIsSet[idx] = false;
  }

  /**
   * Used to unset the time series data at a given index, mainly to be used in case of an error on
   * loading in data, so that it is cleared correctly
   *
   * @param idx Index of data to be removed
   */
  public void removeBlock(int idx) {
    dataBlockArray[idx] = null;
    thisBlockIsSet[idx] = false;
  }

  /**
   * Resample data to a given sample rate
   *
   * @param newSampleRate new sample rate (Hz)
   */
  public void resample(double newSampleRate) {
    long newInterval = (long) (ONE_HZ_INTERVAL / newSampleRate);
    // make sure all data over range gets set to the same interval (and don't upsample)
    for (int i = 0; i < FILE_COUNT; ++i) {
      if (thisBlockIsSet[i]) {
        newInterval = Math.max(newInterval, getBlock(i).getInitialInterval());
      }
    }
    for (int i = 0; i < FILE_COUNT; ++i) {
      if (thisBlockIsSet[i] && getBlock(i).getInitialInterval() != newInterval) {
        getBlock(i).resample(newInterval);
      }
    }
  }

  /**
   * Tests if a given response is loaded or not
   *
   * @param i Index of response to load in
   * @return True if that index has a response loaded in
   */
  public boolean responseIsSet(int i) {
    return thisResponseIsSet[i];
  }

  /**
   * Returns the boolean array where each index shows if there is a response
   * loaded into this datastore object at that index
   * has had data loaded into it (see blockIsSet, dataIsSet methods)
   *
   * @return Array of booleans where an entry is true if a response is loaded
   * at the corresponding index
   */
  private boolean[] responsesAreSet() {
    return thisResponseIsSet;
  }

  /**
   * Adds a pre-constructed datablock to this data store object at the
   * specified index
   *
   * @param idx Index to place the data into
   * @param db Datablock to place into idx
   */
  public void setBlock(int idx, DataBlock db) {
    thisBlockIsSet[idx] = true;
    dataBlockArray[idx] = db;
  }

  public void setBlock(int idx, DataBlock db, int activePlots) throws TimeRangeException {

    thisBlockIsSet[idx] = true;
    dataBlockArray[idx] = db;

    synchronized (this) {
      if (numberOfBlocksSet() > 1) {
        // don't trim data here, that way we don't lose data
        long start = dataBlockArray[idx].getStartTime();
        long end = dataBlockArray[idx].getEndTime();

        // there's clearly already another block loaded, let's make sure they
        // actually have an intersecting time range
        for (int i = 0; i < FILE_COUNT; ++i) {
          if (i != idx && thisBlockIsSet[i]) {
            // whole block either comes before or after the data set
            if (end <= dataBlockArray[i].getInitialStartTime() ||
                start >= dataBlockArray[i].getInitialEndTime()) {

              if (i < activePlots) {
                thisBlockIsSet[idx] = false;
                dataBlockArray[idx] = null;
                throw new TimeRangeException(i + 1);
              } else {
                // unload data that we aren't currently using
                thisBlockIsSet[i] = false;
              }
            }
          }
        }
      }
    }
  }

  /**
   * Loads in a miniseed file with assumption that it is not multiplexed and only includes a single
   * channel's range of data. Mainly used as a shortcut for non-GUI calls (for test cases, etc.)
   *
   * @param idx The plot (range 0 to FILE_COUNT) to be given new data
   * @param filepath Full address of file to be loaded in
   */
  public void setBlock(int idx, String filepath)
      throws SeedFormatException, CodecException,
      IOException {
    String nameFilter = getMplexNameList(filepath).get(0);
    setBlock(idx, filepath, nameFilter, FILE_COUNT);
  }

  /**
   * Takes a loaded miniSEED data series and loads it in as a datablock into
   * this datastore object
   *
   * @param idx The plot (range 0 to FILE_COUNT) to be given new data
   * @param filepath Full address of file to be loaded in
   * @param nameFilter Station ID (SNCL) to load in from multiplexed file
   */
  public void setBlock(int idx, String filepath, String nameFilter)
      throws SeedFormatException, CodecException,
      IOException {
    setBlock(idx, filepath, nameFilter, FILE_COUNT);
  }

  /**
   * Takes a loaded miniSEED data series and loads it in as a datablock into
   * this datastore object. Attempts to fit data in if it overlaps data that is currently unused
   *
   * @param idx The plot (range 0 to FILE_COUNT) to be given new data
   * @param filepath Full address of file to be loaded in
   * @param nameFilter Station ID (SNCL) to load in from multiplexed file
   * @param activePlots Max index of active panel to check as active
   */
  public void setBlock(int idx, String filepath, String nameFilter, int activePlots)
      throws SeedFormatException, CodecException,
      IOException {

    DataBlock xy = getTimeSeries(filepath, nameFilter);
    thisBlockIsSet[idx] = true;
    dataBlockArray[idx] = xy;

    synchronized (this) {
      if (numberOfBlocksSet() > 1) {
        // don't trim data here, that way we don't lose data
        long start = dataBlockArray[idx].getStartTime();
        long end = dataBlockArray[idx].getEndTime();

        // there's clearly already another block loaded, let's make sure they
        // actually have an intersecting time range
        for (int i = 0; i < FILE_COUNT; ++i) {
          if (i != idx && thisBlockIsSet[i]) {
            // whole block either comes before or after the data set
            if (end <= dataBlockArray[i].getInitialStartTime() ||
                start >= dataBlockArray[i].getInitialEndTime()) {

              if (i < activePlots) {
                thisBlockIsSet[idx] = false;
                dataBlockArray[idx] = null;
                throw new TimeRangeException(i+1);
              } else {
                // unload data that we aren't currently using
                thisBlockIsSet[i] = false;
              }
            }
          }
        }
      }
    }

  }

  /**
   * Place an already-constructed instrument response at the index idx
   *
   * @param idx index in this object to place the response at
   * @param ir InstrumentResponse to have placed into this object
   */
  public void setResponse(int idx, InstrumentResponse ir) {
    responses[idx] = ir;
    thisResponseIsSet[idx] = true;
  }

  /**
   * Sets the response of a sensor's dataseries matched by index
   *
   * @param idx Index of plot for which response file matches
   * @param filepath Full address of file to be loaded in
   */
  public void setResponse(int idx, String filepath) throws IOException {
    responses[idx] = new InstrumentResponse(filepath);
    thisResponseIsSet[idx] = true;
  }

  /**
   * Trim all data according to DateTime objects. Converts into epoch milliseconds, which are then
   * used to get the trim range for the underlying datablocks. This trims all data.
   *
   * @param start Start time to trim data to
   * @param end End time to trim data to
   */
  public void trim(Instant start, Instant end) {
    long startTime = start.toEpochMilli();
    long endTime = end.toEpochMilli();
    trim(startTime, endTime, FILE_COUNT);
  }

  /**
   * Trim data according to epoch millisecond longs
   *
   * @param start Start time to trim data to
   * @param end End time to trim data to
   */
  public void trim(long start, long end) {
    trim(start, end, FILE_COUNT);
  }

  /**
   * Trims all data blocks to be within a certain time range.
   * Used for getting a sub-range specified by sliding-bar window.
   *
   * @param start Start time, relative to epoch (nanoseconds)
   * @param end End time, relative to epoch (nanoseconds)
   */
  public void trim(long start, long end, int limit)
      throws IndexOutOfBoundsException {

    // check that the time range is valid to trim all set data
    for (int i = 0; i < limit; ++i) {
      if (!thisBlockIsSet[i]) {
        continue;
      }
      DataBlock db = getBlock(i);

      if (end < db.getStartTime() || start > db.getEndTime()) {

        String trimStartFormatted = formatEpochMillis(start);
        String trimEndFormatted = formatEpochMillis(end);
        String blockStartFormatted = formatEpochMillis(db.getStartTime());
        String blockEndFormatted = formatEpochMillis(db.getEndTime());

        String errMessage = "Trim range outside of valid data window for " + db.getName() + '\n'
            + "Attempted to trim to range (" + trimStartFormatted
            + ", " + trimEndFormatted + ")\n"
            + "Data range is only from (" + blockStartFormatted
            + ", " + blockEndFormatted + ")\n";
        throw new IndexOutOfBoundsException(errMessage);
      }

      if (start < db.getStartTime()) {
        start = db.getStartTime();
      }
      if (end > db.getEndTime()) {
        end = db.getEndTime();
      }
    }

    for (int i = 0; i < FILE_COUNT; ++i) {
      if (thisBlockIsSet[i]) {
        getBlock(i).trim(start, end);
      }
    }
  }

  /**
   * Trims this object's data blocks to hold only points in their common range
   * WARNING: assumes each plot has its data taken at the same point in time
   * (that is, that a common time range exists to be trimmed to)
   */
  public void trimToCommonTime() {
    trimToCommonTime(FILE_COUNT);
  }

  /**
   * Trim the first [limit] blocks of data to a common time range
   *
   * @param limit upper bound of blocks to do trimming on
   */
  public void trimToCommonTime(int limit) {
    // trims the data to only plot the overlapping time of each input set

    if (numberOfBlocksSet() <= 1) {
      return;
    }

    long lastStartTime = Long.MIN_VALUE;
    long firstEndTime = Long.MAX_VALUE;

    // first pass to get the limits of the time data
    for (int i = 0; i < limit; ++i) {
      DataBlock data = dataBlockArray[i];
      if (!thisBlockIsSet[i]) {
        continue;
      }
      long start = data.getStartTime();
      if (start > lastStartTime) {
        lastStartTime = start;
      }
      long end = data.getEndTime();
      if (end < firstEndTime) {
        firstEndTime = end;
      }
    }

    // second pass to trim the data to the limits given
    for (int i = 0; i < limit; ++i) {
      if (!thisBlockIsSet[i]) {
        continue;
      }
      DataBlock data = dataBlockArray[i];
      data.trim(lastStartTime, firstEndTime);
    }
  }

  /**
   * Check if the current trim level (specified by GUI wrapper) encloses all available active data
   * @param start Start time value of zoom/trim set in the GUI
   * @param end End time value of zoom/trim set in the GUI
   * @param limit Number of currently active plots
   * @return true if there is no data that exists beyond current specified trim length
   */
  public boolean currentTrimIsMaximum(long start, long end, int limit) {
    boolean isZoomedOut = false;
    for (int i = 0; i < limit; ++i) {
      if (thisBlockIsSet[i]) {
        // is there at least one block that has a start and end matching the current trim level?
        isZoomedOut |= (start == dataBlockArray[i].getInitialStartTime()
            && end == dataBlockArray[i].getInitialEndTime());
      }
    }
    return isZoomedOut;
  }

  /**
   * Zoom out on current data (include all available)
   *
   * @param limit Number of blocks to return to original time bounds
   */
  public void untrim(int limit) {
    for (int i = 0; i < limit; ++i) {
      if (!thisBlockIsSet[i]) {
        continue;
      }
      DataBlock data = dataBlockArray[i];
      data.untrim();
    }
    trimToCommonTime(limit);
  }

  public void appendBlock(int idx, DataBlock dataBlock, int activePlots) {
    if (!thisBlockIsSet[idx]) {
      setBlock(idx, dataBlock);
      return;
    }

    dataBlockArray[idx].appendTimeSeries(dataBlock);

    synchronized (this) {
      if (numberOfBlocksSet() > 1) {
        // don't trim data here, that way we don't lose data
        long start = dataBlockArray[idx].getStartTime();
        long end = dataBlockArray[idx].getEndTime();

        // there's clearly already another block loaded, let's make sure they
        // actually have an intersecting time range
        for (int i = 0; i < FILE_COUNT; ++i) {
          if (i != idx && thisBlockIsSet[i]) {
            // whole block either comes before or after the data set
            // note that if data ends when another starts, then the data has no overlap --
            // the end time is effectively when the next sample should start
            if (end <= dataBlockArray[i].getInitialStartTime() ||
                start >= dataBlockArray[i].getInitialEndTime()) {

              if (i < activePlots) {
                thisBlockIsSet[idx] = false;
                dataBlockArray[idx] = null;
                throw new TimeRangeException(i+1);
              } else {
                // unload data that we aren't currently using
                thisBlockIsSet[i] = false;
              }
            }
          }
        }
      }
    }
  }

  public void appendBlock(int idx, String filepath, String nameFilter, int activePlots)
      throws SeedFormatException, CodecException,
      IOException {

    if (!thisBlockIsSet[idx]) {
      setBlock(idx, filepath, nameFilter, activePlots);
      return;
    }

    dataBlockArray[idx].appendTimeSeries(filepath);

    synchronized (this) {
      if (numberOfBlocksSet() > 1) {
        // don't trim data here, that way we don't lose data
        long start = dataBlockArray[idx].getStartTime();
        long end = dataBlockArray[idx].getEndTime();

        // there's clearly already another block loaded, let's make sure they
        // actually have an intersecting time range
        for (int i = 0; i < FILE_COUNT; ++i) {
          if (i != idx && thisBlockIsSet[i]) {
            // whole block either comes before or after the data set
            // note that if data ends when another starts, then the data has no overlap --
            // the end time is effectively when the next sample should start
            if (end <= dataBlockArray[i].getInitialStartTime() ||
                start >= dataBlockArray[i].getInitialEndTime()) {

              if (i < activePlots) {
                thisBlockIsSet[idx] = false;
                dataBlockArray[idx] = null;
                throw new TimeRangeException(i+1);
              } else {
                // unload data that we aren't currently using
                thisBlockIsSet[i] = false;
              }
            }
          }
        }
      }
    }
  }

  public class TimeRangeException extends RuntimeException {

    TimeRangeException(int input) {
      super("This data's time range has no overlap with input " + input + ".");
    }

  }

}
