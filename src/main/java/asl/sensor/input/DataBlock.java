package asl.sensor.input;

import java.util.ArrayList;
import java.util.Calendar;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TimeZone;

import org.apache.commons.math3.util.Pair;
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
  
  private long interval, targetInterval;
  private String name;
  private long startTime, endTime;
  private Map<Long, List<Number>> dataMap;
  private long trimmedStart, trimmedEnd;
  private List<Number> cachedTimeSeries;
  boolean rebuildList;
  
  /**
   * Creates a copy of a given DataBlock, which has the same parameters
   * @param in The datablock to be copied
   */
  public DataBlock(DataBlock in) {
    
    interval = in.getInitialInterval();
    targetInterval = in.getInterval();
    dataMap = in.getDataMap();
    name = in.getName();
    startTime = in.getInitialStartTime();
    trimmedStart = in.getStartTime();
    endTime = in.getInitialEndTime();
    trimmedEnd = in.getEndTime();
    
    cachedTimeSeries = in.getData();
    rebuildList = false;
    
  }
  
  /**
   * Creates a trimmed copy of a given DataBlock, with same metadata
   * @param in The datablock to be copied
   * @param start Start time to trim data to
   * @param end End time to trim data to
   */
  public DataBlock(DataBlock in, long start, long end) {
    interval = in.getInitialInterval();
    targetInterval = in.getInterval();
    
    startTime = in.getInitialStartTime();
    trimmedStart = Math.max(startTime, start);
    endTime = in.getInitialEndTime();
    trimmedEnd = Math.min(endTime, end);
    
    dataMap = in.getDataMap();
    name = in.getName();
    
    cachedTimeSeries = in.getData();
    rebuildList = false;
  }
  
  public 
  DataBlock(Map<Long, List<Number>> dataIn, long intervalIn, String nameIn) {
    interval = intervalIn;
    targetInterval = intervalIn;
    
    List<Long> times = new ArrayList<Long>( dataIn.keySet() );
    Collections.sort(times);
    startTime = times.get(0);
    trimmedStart = startTime;
    // TODO: this is wrong
    long lastListStart = times.get( times.size() - 1 );
    int pointsToEnd = dataIn.get(lastListStart).size();
    endTime = lastListStart + (pointsToEnd * intervalIn);
    trimmedEnd = endTime;
    
    name = nameIn;
    dataMap = dataIn;
    
    mergeContiguousTimes();
    rebuildList = true;
  }
  
  public 
  DataBlock(List<Number> dataIn, long intervalIn, String nameIn, long start) {
    interval = intervalIn;
    targetInterval = intervalIn;
    startTime = start;
    dataMap = new HashMap<Long, List<Number>>();
    dataMap.put( startTime, dataIn );
    
    trimmedStart = startTime;
    endTime = interval * ( dataIn.size() - 1 );
    trimmedEnd = endTime;
    
    name = nameIn;
    cachedTimeSeries = new ArrayList<Number>(dataIn);
    rebuildList = false;
  }
  
  public Map<Long, List<Number>> getDataMap() {
    return dataMap;
  }
  
  public void setDataMap(Map<Long, List<Number>> dataMap) {
    this.dataMap = dataMap;
    mergeContiguousTimes();
    rebuildList = true;
  }
  
  public void setData(List<Number> data, long interval, long start) {
    this.interval = interval;
    targetInterval = interval;
    startTime = start;
    trimmedStart = start;
    dataMap = new HashMap<Long, List<Number>>();
    for (int i = 0; i < data.size(); ++i) {
      dataMap.put( startTime, data );
    }
    endTime = interval * ( data.size() - 1 );
    this.interval = interval;
    rebuildList = true;
  }
  
  public void setData(List<Number> data) {
    setData(data, targetInterval, trimmedStart);
    interval = targetInterval;
    rebuildList = true;
  }
  
  public List<Pair<Long, Long>> getGapBoundaries() {
    
    List<Pair<Long, Long>> gapList = new ArrayList<Pair<Long, Long>>();
    
    List<Long> times = new ArrayList<Long>( dataMap.keySet() );
    Collections.sort(times);

    for (int i = 0; i < times.size(); ++i) {
      long timeNow = times.get(i);
      if (timeNow < trimmedStart) {
        continue;
      } else if (timeNow > trimmedEnd) {
        break;
      }
      
      if ( (i + 1) < times.size() ) {
        long blockEnd = ( dataMap.get(timeNow).size() * interval ) + timeNow;
        long timeNext = times.get(i + 1);
        // is there a discrepancy, and is it big enough for us to care?
        if (timeNext - blockEnd > interval * 2) {
           gapList.add( new Pair<Long, Long>(blockEnd, timeNext) );
        }
      }
    }
    
    return gapList;
  }
  
  /**
   * Identifies whether or not input of signal starts positive. Used
   * in step calibration solver to figure out if the data's signs are inverted.
   * @return True if initial step response is negative
   */
  public boolean needsSignFlip() {
    List<Long> times = new ArrayList<Long>( dataMap.keySet() );
    Collections.sort(times);
    int timeIdx = 0;
    for (int i = 0; i < times.size(); ++i) {
      if (times.get(i) <= trimmedStart) {
        if ( times.get(i + 1) > trimmedStart ) {
          timeIdx = i;
        }
      } else {
        break;
      }
    }
    
    long time = times.get(timeIdx);
    List<Number> data = dataMap.get(time);
    
    int startPoint =  
        ( (int) ( trimmedStart - time ) ) / (int) interval;
    
    double max = 0.;
    int maxPoint = startPoint;
    
    for (int i = startPoint; i < ( data.size() - startPoint) / 2; ++i) {
      max = Math.max( max, Math.abs( data.get(i).doubleValue() ) );
      maxPoint = i;
    }
    
    return data.get(maxPoint).doubleValue() > 0;
  }
  
  public List<Number> getData() {
    
    if (!rebuildList) {
      return cachedTimeSeries;
    }
    
    // WARNING: SKIP FACTOR > 0 SHOULD ONLY BE USED TO DISPLAY A REPRESENTATION
    // OF THE DATA, NOT TO PRODUCE DATA TO SEND TO BACKENDS (USE DECIMATION)
    
    List<Number> dataOut = new ArrayList<Number>();
    List<Long> times = new ArrayList<Long>( dataMap.keySet() );
    Collections.sort(times);
    
    long timeCursor = trimmedStart;
    
    // make sure the correct number of points are loaded in to prevent
    // off-by-one errors later on
    int numPoints = 
        (int) ( (trimmedEnd - trimmedStart) / (interval) );
    System.out.println("num. points: " + numPoints);
    
    for (int i = 0; i < times.size(); ++i) {
      
      if ( dataOut.size() == numPoints ) {
        break;
      }
      
      int startIndex;
      long now = times.get(i);
      List<Number> data = dataMap.get(now);
      long next = -1;
      if ( i + 1 < times.size() ) {
        next = times.get(i + 1);
      }
      
      // either we are in the last entry in the map or we need to find the
      // location in the map and corresponding list closest to the current time
      if (timeCursor == now) {
        startIndex = 0;
      } else if ( now < timeCursor && (next > timeCursor || next < 0) ) {
        // get value of time closest to start
        // this is done to deal with the case of differing quantizations
        // between data sets, where all data in a file may be off by less
        // than the interval length, i.e., a millisecond or two
        int closeIdx = (int) ( (timeCursor - now) / interval );
        // this is the index rounded down
        long candidateTime = now + (interval * closeIdx);
        long nextSample = candidateTime + interval;
        if ( Math.abs(timeCursor - candidateTime) 
            <= Math.abs(timeCursor - nextSample) ) {
          startIndex = closeIdx;
        } else {
          startIndex = closeIdx + 1;
        }
        System.out.println("index: " + startIndex);
      } else {
        continue;
      }
      
      System.out.println("Current list length: " + data.size());
      
      if ( startIndex < data.size() ) {
        // make sure we are not in a gap to start with
        int end = startIndex + ( numPoints - dataOut.size() );
        // copy either up to our current end point, or the limit of the block
        end = Math.min( data.size(), end );
        dataOut.addAll( data.subList(startIndex, end) );
      }
      
      timeCursor = trimmedStart + ( interval * dataOut.size() );
      if ( next - timeCursor > (interval * 2) ) {
        // deal with any gaps between two parts of the list
        while (timeCursor < next && dataOut.size() < numPoints) {
          dataOut.add(0.);
          timeCursor += interval;
        }
      }
    }
    
    while (dataOut.size() < numPoints) {
      dataOut.add(0.);
    }
    
    if (interval != targetInterval) {
      dataOut = TimeSeriesUtils.decimate(dataOut, interval, targetInterval);
    }
    
    cachedTimeSeries = dataOut;
    rebuildList = false;
    
    return dataOut;
    
  }
  
  public int size() {
    long timeDiff = trimmedEnd - trimmedStart;
    return (int) (timeDiff / targetInterval);
  }
  
  /**
   * Returns the end time of the data, used mainly in getting range bounds
   * for things like setting end markers, etc.
   * @return The time between two samples of data in microseconds
   */
  public long getEndTime() {
    return trimmedEnd;
  }
  
  /**
   * Get the interval of the data. The timestamp for a given data point in the
   * block can be calculated by startTime + (index * interval).
   * @return The time between two samples of data in microseconds
   */
  public long getInterval() {
    return targetInterval;
  }
  
  public long getInitialInterval() {
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
    return (double) TimeSeriesUtils.ONE_HZ_INTERVAL / (double) targetInterval;
  }

  /**
   * Gives the start timestamp of the miniSEED data. This is a long compatible
   * with the Java System Library's Date and Calendar objects and expressed
   * as microseconds from the UTC epoch
   * @return When the miniSEED data logging started in microseconds 
   */
  public long getStartTime() {
    return trimmedStart;
  }

  public long getInitialStartTime() {
    return startTime;
  }
  
  public long getInitialEndTime() {
    return endTime;
  }

  /**
   * Adjust the target interval of the produced data. This will be used when
   * generating a new series of data from the time series map this object holds.
   * @param newInterval The new interval (time between samples in microseconds)
   */
  public void resample(long newInterval) {
    targetInterval = newInterval;
    rebuildList = true;
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
   * Get start time of data series as a Java calendar object
   * @return Calendar object representing start time (UTC time zone)
   */
  public Calendar getStartCalendar() {
    Calendar cCal = Calendar.getInstance( TimeZone.getTimeZone("UTC") );
    cCal.setTimeInMillis( startTime / 1000 );
    return cCal;
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
    List<Number> data = getData();
    int skipFactor = data.size() / MAX_POINTS + 1; // must be >= 1
    
    // 1000 milliseconds in a microsecond
    long divisor = TimeSeriesUtils.ONE_HZ_INTERVAL / 1000;
    
    XYSeries out = new XYSeries(name);
    long thisTime = trimmedStart;
    for (int i = 0; i < data.size(); i+=skipFactor) {
      Number point = data.get(i);
      double xTime = (double) thisTime / divisor;
      out.add(xTime, point);
      thisTime += skipFactor*targetInterval;
    }
    
    return out;
  }
 
  /**
   * Trim data to a given range (start time, end time)
   * @param start Start time to trim to in microseconds from epoch
   * @param end End time to trim to in microseconds from epoch
   */
  public void trim(long start, long end) {
    
    trimmedStart = Math.max(startTime, start);
    trimmedEnd = Math.min(endTime, end);
    rebuildList = true;
    
  }
  
  private void mergeContiguousTimes() {
    
    // for blocks that start and end at the same point
    List<Long> startTimes = new ArrayList<Long>( dataMap.keySet() );
    Collections.sort(startTimes);
    
    Map<Long, List<Number>> mergedMap = new HashMap<Long, List<Number>>();
    
    int startingPoint = 0;
    int cursor;
    while ( startingPoint < startTimes.size() ) {
      long currentTime = startTimes.get(startingPoint);
      List<Number> currentSeries = dataMap.get(currentTime);
      long lastTimeInList = 
          currentTime + ( currentSeries.size() * interval );
      
      cursor = startingPoint + 1;
      long nextTime = startTimes.get(cursor);
      long difference = Math.abs(nextTime - lastTimeInList);
      
      while ( difference < (interval / 4) ) {
        List<Number> nextGroup = dataMap.get(nextTime);
        currentSeries.addAll(nextGroup);
        lastTimeInList = nextTime + ( nextGroup.size() * interval );
        ++cursor;
        if ( cursor >= startTimes.size() ) {
          break;
        }
        nextTime = startTimes.get(cursor);
        difference = Math.abs(nextTime - lastTimeInList);
      }
      
      mergedMap.put(currentTime, currentSeries);
      startingPoint = cursor + 1;
      
    }
    
    dataMap = mergedMap;
    
  }
  
}
