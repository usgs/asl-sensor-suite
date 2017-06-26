package asl.sensor.input;

import java.util.ArrayList;
import java.util.List;

import org.jfree.data.xy.XYSeries;

import asl.sensor.utils.TimeSeriesUtils;

/**
 * Holds the time series and metadata for a miniSEED file loaded in by the user.
 * Includes methods for resampling data (backed by the TimeSeriesUtils decimate
 * function) and for trimming to specific time regions. 
 * @author akearns
 *
 */
public class DataBlock {
  
  private static final int MAX_POINTS = 100000;
  
  private List<Number> data;
  private long interval;
  private String name;
  private long startTime;
  
  /**
   * Creates a copy of a given DataBlock, which has the same parameters
   * @param in The datablock to be copied
   */
  public DataBlock(DataBlock in) {
    setInterval( in.getInterval() );
    setData( new ArrayList<Number>( in.getData() ) );
    name = in.getName();
    setStartTime(in.getStartTime());
  }
  
  /**
   * Creates a trimmed copy of a given DataBlock, with same metadata
   * @param in The datablock to be copied
   * @param start Start time to trim data to
   * @param end End time to trim data to
   */
  public DataBlock(DataBlock in, long start, long end) {
    setInterval( in.getInterval() );
    int startIdx = in.getTrimStartIndex(start);
    int endIdx = in.getTrimEndIndex(end);
    
    data = in.getData();
    setData( new ArrayList<Number>( data.subList(startIdx, endIdx) ) );
    
    name = in.getName();
    if (startIdx == 0) {
      setStartTime( in.getStartTime() );
    } else {
      setStartTime(start);
    }
    
  }
  
  /**
   * Creates a new DataBlock based on the given parameters
   * @param dataIn Time series data, as a list of numeric objects
   * @param intervalIn The interval between samples (1/sample rate in Hz)
   * @param nameIn The name of the data series, taken from file's metadata
   * @param timeIn The start time of the data series
   */
  public DataBlock
        (List<Number> dataIn, long intervalIn, String nameIn, long timeIn) {
    this.data = dataIn;
    setInterval(intervalIn);
    name = nameIn;
    setStartTime(timeIn);
  }
  

  /**
   * Returns the time series object as a list of Java numeric types.
   * The underlying data can be any numeric type, likely Doubles or Integers
   * @return A list representing the miniSEED time series
   */
  public List<Number> getData() {
    return data;
  }
  
  /**
   * Returns the end time of the data, used mainly in getting range bounds
   * for things like setting end markers, etc.
   * @return The time between two samples of data in microseconds
   */
  public long getEndTime() {
    return startTime + ( interval * data.size() );
  }
  
  /**
   * Get the interval of the data. The timestamp for a given data point in the
   * block can be calculated by startTime + (index * interval).
   * @return The time between two samples of data in microseconds
   */
  public long getInterval() {
    return interval;
  }

  
  /**
   * Get the name of the dataset (used as a key in dataseries, etc.)
   * The name is defined by station-location-channel.
   * @return Name of dataset
   */
  public String getName() {
    return name;
  }
  
  /**
   * Return the sample rate of the data, in Hz. This is the inverse of the
   * interval, scaled to be in seconds rather than microseconds.
   * @return Sample rate in Hz.
   */
  public double getSampleRate() {
    return (double) TimeSeriesUtils.ONE_HZ_INTERVAL / (double) interval;
  }

  /**
   * Gives the start timestamp of the miniSEED data. This is a long compatible
   * with the Java System Library's Date and Calendar objects and expressed
   * as microseconds from the UTC epoch
   * @return When the miniSEED data logging started in microseconds 
   */
  public long getStartTime() {
    return startTime;
  }

  /**
   * Converts an end time to the last index of data to include in a trimmed
   * data series  
   * @param end Terminal cutoff point for data (time, microseconds)
   * @return Index of last data point to include in trimmed subset + 1
   */
  public int getTrimEndIndex(long end) {
    
    int endIdx = 0;
    for (int i = 0; i <= data.size(); ++i) {
      endIdx = i; // exclusive upper bound, so iterate, then check above bound
      if ( startTime + (i * interval) >= end ) {
        break;
      }
    }
    return endIdx;
    
  }

  /**
   * Converts a start time to the first index of data to include in a trimmed
   * data series
   * @param start Initial cutoff point for data (time, microseconds)
   * @return Index of first data point to include in trimmed subset
   */
  public int getTrimStartIndex(long start) {
    int startIdx = 0;
    long endTime = getEndTime();
    if (startTime < start && endTime > start) {
      long diff = start - startTime;
      startIdx = (int) (diff / interval);
    }
    return startIdx;
  }

  /**
   * Checks to see if the sensor's calibration is wired positively or not
   * (i.e., if the result of a step-calibration is upside-down)
   * @return True if sign appears to be incorrect compared to expected step cal
   */
  public boolean needsSignFlip() {
    
    double max = Math.abs( data.get(0).doubleValue() );
    int idx = 0;
    
    for (int i = 1; i < data.size() / 2; ++i) {
      if ( Math.abs( data.get(i).doubleValue() ) > max ) {
        max = Math.abs( data.get(i).doubleValue() );
        idx = i;
      }
    }
    
    return Math.signum( data.get(idx).doubleValue() ) < 0;
    
  }

  /**
   * Replace the time series data with a new list of Java numeric types.
   * This is the preferred call when data is resampled as it requires a new
   * interval to be specified, so that the old sample rate does not persist.
   * @param newInterval The new interval (time between samples in microseconds)
   */
  public void resample(long newInterval) {
    
    if (interval == newInterval) {
      return;
    }
    
    data = TimeSeriesUtils.decimate(data, interval, newInterval);
    interval = newInterval;
  }
  
  /**
   * Replace the time series data with a new list of Java numeric types.
   * If the data has been resampled somehow, the resample method is preferred,
   * as it requires a new interval to be resampled
   * @param dataIn New time series data, as from a miniSEED file
   */
  public void setData(List<Number> dataIn) {
    data = dataIn;
  }

  /**
   * Used to set the interval of the data (to be used, for example, when the
   * time series has had decimation applied)
   * @param interval The time between two samples of data in microseconds
   */
  public void setInterval(long interval) {
    this.interval = interval;
  }
  
  /**
   * Set the start timestamp of the data with a given long, expressed as
   * microseconds from UTC epoch (compatible with Java System Library Date and
   * Calendar objects)
   * @param startTime The start time in microseconds
   */
  public void setStartTime(long startTime) {
    this.startTime = startTime;
  }
  
  /**
   * Alias to function to get size of list in this data block.
   * Equivalent to calling getData().size()
   * @return The size of the list in this structure
   */
  public int size() {
    return data.size();
  }
  
  /**
   * Converts this object's time series data into a form plottable by a chart.
   * The format is a pair of data: the time of a sample and that sample's value.
   * @return JFreeChart XYSeries representation of the data
   */
  public XYSeries toXYSeries() {
    
    // doing a quick decimation here on the displays for datapanel
    // so that we can do the sliding/zooming operations relatively expediently
    // trying to draw the charts with too much data slows it down terribly
    
    int skipFactor = data.size() / MAX_POINTS + 1; // must be >= 1
    
    // 1000 milliseconds in a microsecond
    long divisor = TimeSeriesUtils.ONE_HZ_INTERVAL / 1000;
    
    XYSeries out = new XYSeries(name);
    long thisTime = startTime;
    for (int i = 0; i < data.size(); i+=skipFactor) {
      Number point = data.get(i);
      double xTime = (double) thisTime / divisor;
      out.add(xTime, point);
      thisTime += skipFactor*interval;
    }
    
    return out;
  }

  /**
   * Trim data to a given range (start time, end time)
   * @param start Start time to trim to in microseconds from epoch
   * @param end End time to trim to in microseconds from epoch
   */
  public void trim(long start, long end) {
    
    int startIdx = getTrimStartIndex(start);
    int endIdx = getTrimEndIndex(end);
    
    if ( startIdx == 0 && endIdx >= data.size() ){
      return;
    }
    
    data = new ArrayList<Number>( data.subList(startIdx, endIdx) );
    startTime = start;
    
  }
  
}
