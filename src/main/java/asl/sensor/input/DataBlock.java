package asl.sensor.input;

import java.io.IOException;
import java.time.Instant;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import org.apache.commons.math3.util.Pair;
import org.jfree.data.xy.XYSeries;
import asl.sensor.utils.TimeSeriesUtils;
import edu.iris.dmc.seedcodec.CodecException;
import edu.sc.seis.seisFile.mseed.SeedFormatException;

/**
 * Holds the time series and metadata for a miniSEED file loaded in by the user.
 * Includes methods for resampling data (backed by the TimeSeriesUtils decimate
 * function) and for trimming to specific time regions.
 *
 * Data is stored as a series of contiguous data blocks, stored in a map
 * according to their start time (given in milliseconds from epoch). The data
 * can be set to be windowed to a given region of the valid data, from which
 * the timeseries data can be extracted. The timeseries data is then returned
 * as a list of samples (as double values) over that time range.
 *
 * The contiguous blocks taken from a given miniSEED file should have the same
 * length and start times as given by reading in the same file into a program
 * such as ObsPy. Sample rate is stored as the length of milliseconds between
 * samples (as a long) but is also displayable as the
 *
 * Data names are given as a single string using SNCL (sensor, name, channel,
 * location) metadata, and are generally used within this program as keys
 * for data since experiments generally expect data to be taken over the
 * same time range. Thus two sensors' data fed into an experiment will be
 * considered equal if their names are the same.
 *
 * While it is possible to get data from an instance of this object that
 * includes a gapped range, it is not recommended as the behavior of the time
 * series has not been fully tested. In general most experiments are
 * best performed over time ranges that are fully inside a given contiguous
 * block.
 *
 * @author akearns
 */
public class DataBlock {

  private static final int MAX_POINTS = 5000;

  private final long interval;
  private final String name;
  private long targetInterval;
  private long startTime, endTime;
  private Map<Long, double[]> dataMap;
  private long trimmedStart, trimmedEnd;
  private double[] cachedTimeSeries;
  private boolean rebuildList;

  /**
   * Creates a copy of a given DataBlock, which has the same parameters
   *
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

    cachedTimeSeries = in.getData().clone();
    //Set to false, as it was just rebuilt in getData()
    rebuildList = false;

  }

  /**
   * Create a new DataBlock using a single contiguous block of data; the
   * resulting datamap will have a single entry at the specified start time
   * which points to the passed-in array of data.
   *
   * @param dataIn Timeseries data to be loaded into this datablock
   * @param intervalIn Sampling interval of data
   * @param nameIn SNCL metadata of data source
   * @param start Time (in ms from epoch) that the data begins on
   */
  public DataBlock(double[] dataIn, long intervalIn, String nameIn, long start) {
    interval = intervalIn;
    targetInterval = intervalIn;
    startTime = start;
    dataMap = new LinkedHashMap<>();
    dataMap.put(startTime, dataIn);

    trimmedStart = startTime;
    // subt. 1 from length as 'end time' should be identical to time of last point, not
    // when the point after it is supposed to start
    endTime = startTime + (interval * (dataIn.length - 1));
    trimmedEnd = endTime;

    name = nameIn;
    cachedTimeSeries = dataIn;
    rebuildList = false;
  }

  /**
   * Create a new datablock by feeding in specific parameters
   *
   * @param dataIn Map of contiguous data blocks (which may require merging)
   * @param intervalIn Sampling interval of data in ms
   * @param nameIn SNCL metadata of data source
   */
  public DataBlock(Map<Long, double[]> dataIn, long intervalIn, String nameIn) {
    interval = intervalIn;
    targetInterval = intervalIn;

    name = nameIn;
    dataMap = dataIn;

    mergeContiguousTimes();
    recalculateTimes();
    rebuildList = true;
  }

  /**
   * If necessary, construct an array representing the data in the given window
   * from specified start and end times and return it. The given array is
   * cached until either it, the underlying data map, the time range of the
   * window, or desired output sample rate (requiring decimation) change.
   * If the current time window includes gaps, these will be populated by zeros.
   * It is not recommended to trim data to a region including these gaps because
   * they can produce undesired behavior in the results of experiments. If the
   * data time range needs to be reduce it will also perform decimation.
   *
   * @return Array representing the data found within a given time range
   */
  public double[] getData() {

    if (!rebuildList) {
      return cachedTimeSeries;
    }

    List<Long> times = new ArrayList<>(dataMap.keySet());
    Collections.sort(times);

    long timeCursor = trimmedStart;

    int numPoints = sizeNativeSampleRate();
    cachedTimeSeries = new double[numPoints];
    int lastFilledIndex = 0;

    for (int i = 0; i < times.size(); ++i) {

      if (lastFilledIndex == numPoints) {
        break;
      }

      int startIndex;
      long now = times.get(i);
      double[] data = dataMap.get(now);
      long next = -1;
      if (i + 1 < times.size()) {
        next = times.get(i + 1);
      }

      // either we are in the last entry in the map or we need to find the
      // location in the map and corresponding list closest to the current time
      if (timeCursor == now) {
        startIndex = 0;
      } else if (now < timeCursor && (next > timeCursor || next < 0)) {
        // get value of time closest to start
        // this is done to deal with the case of differing quantizations
        // between data sets, where all data in a file may be off by less
        // than the interval length, i.e., a millisecond or two
        int closeIdx = (int) ((timeCursor - now) / interval);
        // this is the index rounded down
        long candidateTime = now + (closeIdx * interval);
        long nextSample = candidateTime + interval;
        if (Math.abs(timeCursor - candidateTime)
            <= Math.abs(timeCursor - nextSample)) {
          startIndex = closeIdx;
        } else {
          startIndex = closeIdx + 1;
        }
      } else {
        continue;
      }

      if (startIndex < data.length) {
        // make sure we are not in a gap to start with
        int end = startIndex + (numPoints - lastFilledIndex);
        // copy either up to our current end point, or the limit of the block
        end = Math.min(data.length, end);
        double[] sublist = Arrays.copyOfRange(data, startIndex, end);
        for (int j = 0; j < sublist.length; ++j, ++lastFilledIndex) {
          cachedTimeSeries[lastFilledIndex] = sublist[j];
        }

      }

      timeCursor = trimmedStart + (interval * lastFilledIndex);
      if (next - timeCursor > (interval * 2)) {
        // deal with any gaps between two parts of the list
        while (timeCursor < next && lastFilledIndex < numPoints) {
          cachedTimeSeries[lastFilledIndex] = 0.;
          ++lastFilledIndex;
          timeCursor += interval;
        }
      }
    }

    while (lastFilledIndex < numPoints) {
      cachedTimeSeries[lastFilledIndex] = 0.;
      ++lastFilledIndex;
    }

    if (interval != targetInterval) {
      cachedTimeSeries =
          TimeSeriesUtils.decimate(cachedTimeSeries, interval, targetInterval);
      // prevent issues cause by rounding of time series length after downsampling
      int finalNumPoints =
          (int) Math.ceil((trimmedEnd - trimmedStart) / ((double) targetInterval));
      if (finalNumPoints < cachedTimeSeries.length) {
        cachedTimeSeries = Arrays.copyOfRange(cachedTimeSeries, 0, finalNumPoints);
      }
    }

    rebuildList = false;
    return cachedTimeSeries;
  }

  /**
   * Return a copy of the underlying data structure of this object. It is not
   * returned directly because modification of the underlying object would
   * invalidate any cached timeseries array.
   * To modify the datamap inside this object, call this function and then
   * make a call to the corresponding setter to replace it. @see #setDataMap
   * The data structure is a map of contiguous block start times to the
   * timeseries data of that block (as an array)
   *
   * @return copy of this datablock's underlying contiguous block map
   */
  public Map<Long, double[]> getDataMap() {
    return new LinkedHashMap<>(dataMap);
  }

  /**
   * Returns the end time of the data, used mainly in getting range bounds
   * for things like setting end markers, etc.
   *
   * @return The time between two samples of data in milliseconds
   */
  public long getEndTime() {
    return trimmedEnd;
  }

  /**
   * Output the boundaries between contiguous blocks inside the currently
   * selected data window. Used to display the input panel's orange regions.
   * This is given as a list of paired start and end times specified in ms
   * from epoch.
   *
   * @return List of time ranges specifying lengths of time with no samples
   */
  public List<Pair<Long, Long>> getGapBoundaries() {

    List<Pair<Long, Long>> gapList = new ArrayList<>();

    List<Long> times = new ArrayList<>(dataMap.keySet());
    Collections.sort(times);
    // contiguous blocks must have been merged for this to work correctly!
    for (int i = 0; i < times.size(); ++i) {
      long timeNow = times.get(i);
      long blockEnd = dataMap.get(timeNow).length * interval + timeNow;
      boolean hasNext = (i + 1) < times.size();

      if (blockEnd < trimmedStart) {
        // does data (re-)start before our trimmed region does?
        // if not, data begins with a gap
        if (hasNext && times.get(i + 1) > trimmedStart) {
          // does the next data point start before our region of interest ends?
          long gapEnd = Math.min(times.get(i + 1), trimmedEnd);
          gapList.add(new Pair<>(trimmedStart, gapEnd));
        } else if (!hasNext) {
          gapList.add(new Pair<>(trimmedStart, trimmedEnd));
        }
        continue;
      }
      if (timeNow > trimmedEnd) {
        break;
      }
      // check if a gap exists completely inside our selection window
      if ((i + 1) < times.size()) {
        long timeNext = times.get(i + 1);
        // is there a discrepancy, and is it big enough to be a gap?
        if (timeNext - blockEnd > (3 * interval) / 2) {
          long gapEnd = Math.min(timeNext, trimmedEnd);
          gapList.add(new Pair<>(blockEnd, gapEnd));
        }
      }
    }

    return gapList;
  }

  /**
   * Gives the end timestamp of the miniSEED data. This is a long compatible
   * with the Java System Library's Date and Calendar objects and expressed
   * as milliseconds from the UTC epoch
   *
   * @return Time after the last sample of the latest contiguous block in ms
   */
  public long getInitialEndTime() {
    return endTime;
  }

  /**
   * Get the sampling interval of the internal representation of the data.
   * Mainly useful in terms of determining when data needs to be downsampled.
   *
   * @return Interval between samples in ms
   */
  public long getInitialInterval() {
    return interval;
  }

  /**
   * Gives the start timestamp of the miniSEED data. This is a long compatible
   * with the Java System Library's Date and Calendar objects and expressed
   * as milliseconds from the UTC epoch
   *
   * @return When the miniSEED data logging started in milliseconds
   */
  public long getInitialStartTime() {
    return startTime;
  }

  /**
   * Get the interval of the output data. The timestamp for a given data point
   * in the block can be calculated by startTime + (index * interval).
   *
   * @return The time between two samples of data in milliseconds
   */
  public long getInterval() {
    return targetInterval;
  }

  /**
   * Get the name of the dataset (used as a key in dataseries, etc.)
   * The name is defined by station-location-channel.
   *
   * @return Name of dataset
   */
  public String getName() {
    return name;
  }

  /**
   * Return the sample rate of the data, in Hz. This is the inverse of the
   * interval, scaled to be in seconds rather than milliseconds.
   *
   * @return Sample rate in Hz.
   */
  public double getSampleRate() {
    return (double) TimeSeriesUtils.ONE_HZ_INTERVAL / (double) interval;
  }

  /**
   * Get (untrimmed) start time of the data
   *
   * @return DateTime object representing start time in UTC time zone
   */
  public Instant getStartInstant() {
    return Instant.ofEpochMilli(startTime);
  }

  /**
   * Gives the start timestamp of the trim window. This is a long compatible
   * with the Java System Library's Date and Calendar objects and expressed
   * as milliseconds from the UTC epoch
   *
   * @return When the miniSEED data logging started in milliseconds
   */
  public long getStartTime() {
    return trimmedStart;
  }

  /**
   * Find contiguous blocks of data and merge into a single series. If there are
   * duplicated data points, ignore them.
   */
  private void mergeContiguousTimes() {

    // for blocks that start and end at the same point
    List<Long> startTimes = new ArrayList<>(dataMap.keySet());
    Collections.sort(startTimes);

    Map<Long, double[]> mergedMap = new LinkedHashMap<>();

    int startingPoint = 0;
    int cursor;
    while (startingPoint < startTimes.size()) {
      List<double[]> toMerge = new ArrayList<>();
      long currentTime = startTimes.get(startingPoint);

      double[] currentSeries = dataMap.get(currentTime);
      toMerge.add(currentSeries);
      long timeAtSublistEnd =
          currentTime + (currentSeries.length * interval);

      cursor = startingPoint + 1;

      if (cursor >= startTimes.size()) {
        mergedMap.put(currentTime, currentSeries);
        dataMap = mergedMap;
        return;
      }

      long nextTime = startTimes.get(cursor);

      long difference = nextTime - timeAtSublistEnd;
      while (difference < (interval / 4)) {

        if (difference < 0) {
          // duplicated data check
          long diff = timeAtSublistEnd - nextTime; // take the difference
          long mod = diff % interval; // round up if greater than 75% interval

          // division truncates; if the result would be, say, 1.9 as a decimal,
          // then we should round up to 2 instead of start at 1 since this is
          // common enough with timing differences between gaps
          // whereas if this conditional was triggered beacuse the next block
          // occurs at, say, 90% of the interval, then the difference is only
          // 0.1 from the expected start and we copy everything anyway
          int fstUndupIdx = (int) (diff / interval);
          if (mod > 3 * interval / 4) {
            ++fstUndupIdx;
          }

          double[] next = dataMap.get(nextTime);
          if (fstUndupIdx < next.length) {
            // is there data between list end and current time cursor?
            // there is some data where records are 100% redundant and
            // have no data that can be loaded in from that record
            // (the unduplicated index is past the length of the data)
            // the conditional would be false -- and this would trigger an error
            double[] truncated =
                Arrays.copyOfRange(next, fstUndupIdx, next.length);
            toMerge.add(truncated);
            timeAtSublistEnd = timeAtSublistEnd + (truncated.length * interval);
          }

          ++cursor;
        } else {
          // data not duplicated, so copy it all
          double[] nextGroup = dataMap.get(nextTime);
          toMerge.add(nextGroup);
          // currentSeries = TimeSeriesUtils.addAll(currentSeries, nextGroup);
          timeAtSublistEnd = nextTime + (nextGroup.length * interval);
          ++cursor;
        }

        if (cursor >= startTimes.size()) {
          double[] contiguousSeries = TimeSeriesUtils.concatAll(toMerge);
          mergedMap.put(currentTime, contiguousSeries);
          dataMap = mergedMap;
          return;
        }

        nextTime = startTimes.get(cursor);
        difference = nextTime - timeAtSublistEnd;

      }

      // end of the contiguous block. merge data and iterate
      double[] contiguousSeries = TimeSeriesUtils.concatAll(toMerge);
      mergedMap.put(currentTime, contiguousSeries);
      startingPoint = cursor;

    }
    dataMap = mergedMap;
  }

  /**
   * Identifies whether or not input of signal starts positive. Used
   * in step calibration solver to figure out if the data's signs are inverted.
   *
   * @return True if initial step response is negative
   */
  public boolean needsSignFlip() {
    return TimeSeriesUtils.needsSignFlip(getData());
  }

  /**
   * Adjust the target interval of the produced data. This will be used when
   * generating a new series of data from the time series map this object holds.
   * The specified interval will be ignored if it represents a higher sample
   * rate than the current sample; only decimation is performed, not upsampling.
   *
   * @param newInterval The new interval (time between samples in milliseconds)
   */
  public void resample(long newInterval) {
    targetInterval = Math.max(interval, newInterval);
    rebuildList = rebuildList || (targetInterval != interval);
  }

  /**
   * Return the length of datapoints to be contained in the continuous data
   * the program will return given by the currently-specified time window
   *
   * @return length of the data region to be returned, including filled-in gaps
   */
  public int size() {
    long timeDiff = trimmedEnd - trimmedStart;
    return (int) (timeDiff / targetInterval);
  }

  /**
   * Return the number of points in the given data range assuming the data does not require
   * any level of decimation in order to be read in
   *
   * @return length of the data region to be returned, including filled-in gaps
   */
  int sizeNativeSampleRate() {
    long timeDiff = trimmedEnd - trimmedStart;
    return (int) (timeDiff / interval);
  }


  /**
   * Converts this object's time series data into a form plottable by a chart.
   * The format is a pair of data: the time of a sample and that sample's value.
   * The data is truncated (without smoothing) to prevent the plot from slowing
   * down the performance of the program as much as possible.
   *
   * @return JFreeChart XYSeries representation of the data
   */
  public XYSeries toXYSeries() {
    // doing a quick decimation here on the displays for datapanel
    // so that we can do the sliding/zooming operations relatively expediently
    // trying to draw the charts with too much data slows it down terribly
    double[] data = getData();

    int skipFactor = data.length / MAX_POINTS + 1; // must be >= 1

    XYSeries out = new XYSeries(name);
    long thisTime = trimmedStart;
    for (int i = 0; i < data.length; i += skipFactor) {
      double point = data[i];
      double xTime = thisTime;
      out.add(xTime, point);
      thisTime += skipFactor * targetInterval;
    }

    return out;
  }

  /**
   * Adjust the window of data to collect samples from when getting the data
   *
   * @param start Start time to trim window to in milliseconds from epoch
   * @param end End time to trim window to in milliseconds from epoch
   */
  public void trim(long start, long end) {

    long temp = Math.min(start, end);
    end = Math.max(start, end);
    start = temp;

    trimmedStart = Math.max(startTime, start);
    trimmedEnd = Math.min(endTime, end);
    rebuildList = rebuildList ||
        (startTime != trimmedStart) || (endTime != trimmedEnd);

  }

  /**
   * Reset data window to entire region specified by underlying datamap.
   */
  public void untrim() {
    boolean regen = (trimmedStart != startTime) || (trimmedEnd != endTime);
    trimmedStart = startTime;
    trimmedEnd = endTime;
    rebuildList = rebuildList || regen;

  }

  private void recalculateTimes() {
    List<Long> times = new ArrayList<>(dataMap.keySet());
    Collections.sort(times);
    startTime = times.get(0);
    trimmedStart = startTime;
    long lastListStart = times.get(times.size() - 1);
    // as previously seen, subtract 1 to get time of last point, not expected start of next point
    int pointsToEnd = dataMap.get(lastListStart).length - 1;
    endTime = lastListStart + (pointsToEnd * interval);
    trimmedEnd = endTime;
    rebuildList = true;
  }

  public void appendTimeSeries(String filepath)
      throws IOException, SeedFormatException, CodecException {
    Map<Long, double[]> toAppend =
        TimeSeriesUtils.getTimeSeriesMap(filepath, name).getSecond();
    dataMap.putAll(toAppend);
    mergeContiguousTimes();
    recalculateTimes();
  }

}
