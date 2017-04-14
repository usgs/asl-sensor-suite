package asl.sensor.input;

import java.io.FileNotFoundException;
import java.io.IOException;

import org.jfree.data.xy.XYSeries;

import asl.sensor.utils.FFTResult;
import asl.sensor.utils.TimeSeriesUtils;

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
 *
 */
public class DataStore {

  /**
   * Defines the maximum number of plots to be shown
   */
  public final static int FILE_COUNT = 9;
  private DataBlock[] dataBlockArray;
  private InstrumentResponse[] responses;
  private XYSeries[] outToPlots; // used to cache graph data
  
  // these are used to check to make sure data has been loaded
  private boolean[] thisBlockIsSet;
  private boolean[] thisResponseIsSet;
  
  /**
   * Instantiate the collections, including empty datasets to be sent to
   * charts for plotting (see DataPanel)
   */
  public DataStore() {
   dataBlockArray = new DataBlock[FILE_COUNT];
   responses = new InstrumentResponse[FILE_COUNT];
   outToPlots = new XYSeries[FILE_COUNT];
   thisBlockIsSet = new boolean[FILE_COUNT];
   thisResponseIsSet = new boolean[FILE_COUNT];
   for (int i = 0; i < FILE_COUNT; ++i) {
     outToPlots[i] = new XYSeries("(EMPTY) " + i);
     thisBlockIsSet[i] = false;
     thisResponseIsSet[i] = false;
   }
  }
  
  /**
   * Create a full copy of the current datastore
   * @param ds datastore to copy
   */
  public DataStore(DataStore ds) {
    this(ds, FILE_COUNT);
  }
  
  /**
   * Create a copy of only some of the current datastore, up to upperBound
   * @param ds datastore to copy
   * @param upperBound max index to copy data from
   */
  public DataStore(DataStore ds, int upperBound) {
    dataBlockArray = new DataBlock[FILE_COUNT];
    responses = new InstrumentResponse[FILE_COUNT];
    thisBlockIsSet = new boolean[FILE_COUNT];
    thisResponseIsSet = new boolean[FILE_COUNT];
    outToPlots = new XYSeries[FILE_COUNT];
    boolean[] setBlocks = ds.dataIsSet();
    boolean[] setResps = ds.responsesAreSet();
    for (int i = 0; i < upperBound; ++i) {
      if (setBlocks[i]) {
        dataBlockArray[i] = new DataBlock( ds.getBlock(i) );
        thisBlockIsSet[i] = true;
        outToPlots[i] = ds.getPlotSeries(i);
      }
      if (setResps[i]) {
        responses[i] = ds.getResponse(i);
        thisResponseIsSet[i] = true;
      }
    }
  }
  
  /**
   * Create a copy of a datablock object that has been trimmed to a
   * smaller range for the use of showing plot zooms, etc.
   * @param ds DataStore object to be copied and trimmed from
   * @param start Start time trim to
   * @param end End time to trim to
   */
  public DataStore(DataStore ds, long start, long end, int upperBound) {
    this(ds);
    
    this.trimAll(start, end);
  }
  
  /**
   * Checks whether the data block (timeseries) at a given index
   * contains data or is null/empty
   * @param idx Index of the data to look at (between 0 and DATA_STORE)
   * @return True if a seed file has been loaded in at the index
   */
  public boolean blockIsSet(int idx) {
    return thisBlockIsSet[idx];
  }
  
  /**
   * Checks if both components at a specific index are set
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
   * @return Array of booleans where an entry is true if data is loaded at
   * the corresponding index
   */
  public boolean[] dataIsSet() {
    return thisBlockIsSet;
  }
  
  /**
   * Return a single data block according to the passed index 
   * @param idx Index of datablock, corresponding to data panel plot index
   * @return Timeseries data for corresponing plot
   */
  public DataBlock getBlock(int idx) {
    return dataBlockArray[idx];
  }
  
  /**
   * Returns the set of structures used to hold the loaded miniSeed data sets
   * @return An array of DataBlocks (time series and metadata)
   */
  public DataBlock[] getData() {
    return dataBlockArray;
  }
  
  /**
   * Returns the plottable format of the data held in the arrays at 
   * the specified index
   * @param idx Which of this structure's plottable series should be loaded
   * @return The time series data at given index, to be sent to a chart
   */
  public XYSeries getPlotSeries(int idx) {
    return outToPlots[idx];
  }
  
  /**
   * Gets the power-spectral density of an index in this object.
   * If a PSD has already been calculated, this will return that. If not,
   * it will calculate the result, store it, and then return that data.
   * @param idx Index of data to get the PSD of
   * @return Complex array of frequency values and a 
   * double array of the frequencies
   */
  public synchronized FFTResult getPSD(int idx) {
    DataBlock db = dataBlockArray[idx];
    InstrumentResponse ir = responses[idx];
    return FFTResult.crossPower(db, db, ir, ir);
  }
  
  /**
   * Get the instrument response object at a given index
   * @param idx Index to get the response for
   * @return The instrument response data (gain, poles, zeros, etc.)
   */
  public InstrumentResponse getResponse(int idx) {
    return responses[idx];
  }
  
  /**
   * Returns the instrument responses as an array, ordered such that the
   * response at index i is the response associated with the DataBlock at i
   * in the array of DataBlocks
   * @return Array of instrument responses
   */
  public InstrumentResponse[] getResponses() {
    return responses;
  }
  
  /**
   * Used to get the first, second, etc. data set loaded. Used when operations
   * reading in data don't require all the inputs to be loaded.
   * Requires both SEED and RESP to be loaded for this to be valid.
   * @param x x-th set of loaded data to get, starting at 1 (NOT 0) 
   * @return index of the loaded data
   */
  public int getXthFullyLoadedIndex(int x) {
    if (x < 1) {
      throw new IndexOutOfBoundsException("Parameter must be >= 1");
    }
    
    int loaded = 0;
    for (int i = 0; i < FILE_COUNT; ++i) {
      if ( bothComponentsSet(i) ) {
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
   * @return True if there is anything loaded into this datastore object.
   */
  public boolean isAnythingSet() {
    for (int i = 0; i< FILE_COUNT; ++i) {
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
   * @param limit Index of last block to match intervals to
   */
  public void matchIntervals(int limit){
    long interval = 0;
    // first loop to get lowest-frequency data
    for (int i = 0; i < limit; ++i) {
      if ( thisBlockIsSet[i] && getBlock(i).getInterval() > interval ) {
        interval = getBlock(i).getInterval();
      }
    }
    // second loop to downsample
    for (int i = 0; i < limit; ++i) {
      if ( thisBlockIsSet[i] && getBlock(i).getInterval() != interval ) {
        getBlock(i).resample(interval);
      }
    }
    
    trimToCommonTime();
    
  }

  /**
   * Gives the count of indices where both a miniseed and response are loaded
   * @return the number of entries of miniseeds with a matching response
   */
  public int numberFullySet() {
    int loaded = 0;
    for (int i = 0; i < FILE_COUNT; ++i) {
      if ( bothComponentsSet(i) ) {
        ++loaded;
      }
    }
    return loaded;
  }
  
  /**
   * Gives the number of DataBlocks (miniseed files) read in to this object
   * @return number of files read in
   */
  public int numberOfBlocksSet() {
    int loaded = 0;
    for (int i = 0; i < FILE_COUNT; ++i) {
      if (thisBlockIsSet[i]){
        ++loaded;
      }
    }
    return loaded;
  }
  
  /**
   * Removes all data at a specific index -- miniseed, response, and any
   * data generated from them
   * @param idx
   */
  public void removeData(int idx) {
    dataBlockArray[idx] = null;
    responses[idx] = null;
    outToPlots[idx] = null;
    thisBlockIsSet[idx] = false;
    thisResponseIsSet[idx] = false;
  }
  
  /**
   * Tests if a given response is loaded or not
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
   * @return Array of booleans where an entry is true if a response is loaded 
   * at the corresponding index
   */
  public boolean[] responsesAreSet() {
    return thisResponseIsSet;
  }
  
  /**
   * Adds a pre-constructed datablock to this data store object at the 
   * specified index
   * @param idx Index to place the data into
   * @param db Datablock to place into idx
   */
  public void setData(int idx, DataBlock db) {
    thisBlockIsSet[idx] = true;
    dataBlockArray[idx] = db;
    outToPlots[idx] = db.toXYSeries();
  }
  
  /**
   * Takes a loaded miniSEED data series and converts it to a plottable format
   * @param idx The plot (range 0 to FILE_COUNT) to be given new data
   * @param filepath Full address of file to be loaded in
   * @return The miniSEED data, in a plottable format
   */
  public void setData(int idx, String filepath, String nameFilter) {
    
    try {
      DataBlock xy = TimeSeriesUtils.getTimeSeries(filepath, nameFilter);
      thisBlockIsSet[idx] = true;
      dataBlockArray[idx] = xy;
      outToPlots[idx] = xy.toXYSeries();
    } catch (FileNotFoundException e) {
      e.printStackTrace();
    }

    synchronized(this) {
      if (numberOfBlocksSet() > 1) {
        // don't trim data here, that way we don't lose data
        long start = dataBlockArray[idx].getStartTime();
        long end = dataBlockArray[idx].getEndTime();

        // there's clearly already another block loaded, let's make sure they
        // actually have an intersecting time range
        for (int i = 0; i < FILE_COUNT; ++i) {
          if (i != idx && thisBlockIsSet[i]) {
            // whole block either comes before or after the data set
            if (end < dataBlockArray[i].getStartTime() || 
                start > dataBlockArray[i].getEndTime() ) {
              thisBlockIsSet[idx] = false;
              outToPlots[idx] = null;
              dataBlockArray[idx] = null;
              throw new RuntimeException("Time range does not intersect");
            }
          }
        }
      }
    }
    
  }
  
  /**
   * Place an already-constructed instrument response at the index idx 
   * @param idx index in this object to place the response at
   * @param ir InstrumentResponse to have placed into this object
   */
  public void setResponse(int idx, InstrumentResponse ir) {
    responses[idx] = ir;
    thisResponseIsSet[idx] = true;
  }

  /**
   * Sets the response of a sensor's dataseries matched by index
   * @param idx Index of plot for which response file matches
   * @param filepath Full address of file to be loaded in
   */
  public void setResponse(int idx, String filepath) {
    try {
      responses[idx] = new InstrumentResponse(filepath);
      thisResponseIsSet[idx] = true;
    } catch (IOException e) {
      e.printStackTrace();
    }
  }
  
  /**
   * Alias to blockIsSet function
   * @param idx Index of a datablock to check
   * @return True if a seed file has been loaded in there
   */
  public boolean timeSeriesSet(int idx) {
    return thisBlockIsSet[idx];
  }
  
  /**
   * Trims all data blocks to be within a certain time range.
   * Used for getting a sub-range specified by sliding-bar window.
   * @param start Start time, relative to epoch (nanoseconds)
   * @param end End time, relative to epoch (nanoseconds)
   */
  public void trimAll(long start, long end) {
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
   * @param limit upper bound of blocks to do trimming on
   */
  public void trimToCommonTime(int limit) {
    // trims the data to only plot the overlapping time of each input set
    
    if ( numberOfBlocksSet() <= 1 ) {
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
      long end = start + data.getInterval() * data.size();
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
      outToPlots[i] = data.toXYSeries();
    }
    
  }
}
