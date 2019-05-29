package asl.sensor.utils;

import edu.iris.dmc.seedcodec.Codec;
import edu.sc.seis.seisFile.SeisFileException;
import edu.sc.seis.seisFile.fdsnws.FDSNDataSelectQuerier;
import edu.sc.seis.seisFile.fdsnws.FDSNDataSelectQueryParams;
import edu.sc.seis.seisFile.mseed.DataRecordIterator;
import java.io.DataInputStream;
import java.io.EOFException;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.time.Instant;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TimeZone;
import org.apache.commons.math3.util.Pair;
import asl.sensor.input.DataBlock;
import edu.iris.dmc.seedcodec.CodecException;
import edu.sc.seis.seisFile.mseed.Blockette;
import edu.sc.seis.seisFile.mseed.Blockette1000;
import edu.sc.seis.seisFile.mseed.Btime;
import edu.sc.seis.seisFile.mseed.DataHeader;
import edu.sc.seis.seisFile.mseed.DataRecord;
import edu.sc.seis.seisFile.mseed.SeedFormatException;
import edu.sc.seis.seisFile.mseed.SeedRecord;

/**
 * Contains static methods for grabbing data from miniSEED files
 * and some very basic timeseries processing tools (i.e., decimation)
 *
 * @author akearns
 */
public class TimeSeriesUtils {


  /**
   * Interval for data that has been sampled at 1 Hz in milliseconds
   */
  public final static long ONE_HZ_INTERVAL = 1000L;

  /**
   * Sample rate of a 1 Hz sample, in Hz, as a double (that is, 1.0)
   */
  public final static double ONE_HZ = 1.0;

  // divide by this to go from nanoseconds to milliseconds
  public static final int TO_MILLI_FACTOR = 1000000;

  public static final ThreadLocal<SimpleDateFormat> DATE_TIME_FORMAT =
      ThreadLocal.withInitial(() -> {
        SimpleDateFormat format = new SimpleDateFormat("YYYY.DDD.HH:mm:ss.SSS");
        format.setTimeZone(TimeZone.getTimeZone("UTC"));
        return format;
      });

  public static String formatEpochMillis(long millis) {
    return DATE_TIME_FORMAT.get().format(Date.from(Instant.ofEpochMilli(millis)));
  }

  /**
   * Merge arrays from multiple timeseries into a single object
   *
   * @param arrays Series of arrays
   * @return All inputted arrays concatenated into a single array (can be specified as a 2D array)
   */
  public static double[] concatAll(double[]... arrays) {

    if (arrays.length == 0) {
      return new double[]{};
    }

    if (arrays.length == 1) {
      return arrays[0];
    }

    int len = 0;
    for (double[] arr : arrays) {
      len += arr.length;
    }

    double[] result = new double[len];
    int start = 0;
    for (double[] arr : arrays) {
      if (arr.length == 0) {
        continue;
      }
      int end = arr.length;
      System.arraycopy(arr, 0, result, start, end);
      start += end;
    }

    return result;
  }

  /**
   * Merge a series of arrays into a single array in order. Used to concatenate
   * all contiguous arrays in a data map. If the list is empty, an empty array
   * is returned. If the list has one object, that object is returned. Otherwise
   * a new array is constructed holding all the inputted arrays at once. The
   * size of the combined array is equal to the total length of all arrays in
   * the list.
   *
   * @param list List of arrays to merge all together into a single array
   * @return Single array holding each array's data in sequence.
   */
  public static double[] concatAll(List<double[]> list) {
    return concatAll(list.toArray(new double[list.size()][]));
  }

  /**
   * Initial driver for the decimation utility
   * which takes a timeseries of unknown rate and
   * runs downsampling to convert it to a target
   * frequency of a 1Hz interval.
   *
   * @param data The timeseries to be decimated
   * @param source The source frequency as interval between samples (milliseconds)
   * @return A timeseries decimated to the correct frequency
   */
  public static double[] decimate(double[] data, long source, long target) {

    // a sample lower than 1Hz frq has longer time between samples
    // since it's an inverse relationship and all
    if (source >= target) {
      // if data is too low-frequency to decimate, do nothing
      return data;
    }

    // find what the change in size is going to be
    long gcd = euclidGCD(source, target);
    // conversion up- and down-factors
    // (upsample by target, downsample by source)
    // cast is valid because any discrete interval
    // from 1Hz and up is already expressable
    // as an int
    int upf = (int) (source / gcd);
    int dnf = (int) (target / gcd);

    double higherFreq = (1. / source) * upf * ONE_HZ_INTERVAL;
    double lowerFreq = (1. / target) * ONE_HZ_INTERVAL / 2;
    // nyquist rate of downsampled data

    // one valid sample rate for data is 2.5Hz
    // with 1Hz that comes out as a ratio of 5/2, which won't
    // downsample neatly in some cases so we would first upsample,
    // filter out any noise terms, then downsample

    double[] upped = upsample(data, upf);
    double[] lpfed = FFTResult.lowPassFilter(upped, higherFreq, lowerFreq);

    return downsample(lpfed, dnf);

  }

  /**
   * Remove mean (constant value) from a dataset and include
   *
   * @return timeseries as numeric list with previous mean subtracted
   */
  public static double[] demean(double[] dataSet) {
    double[] dataOut = dataSet.clone();
    TimeSeriesUtils.demeanInPlace(dataOut);
    return dataOut;
  }

  /**
   * In-place subtraction of mean from each point in an incoming data set.
   * This is a necessary step in calculating the power-spectral density.
   *
   * @param dataSet The data to have the mean removed from.
   */
  public static void demeanInPlace(double[] dataSet) {

    // I'm always getting the demeaning tasks, huh?

    if (dataSet.length == 0) {
      return; // shouldn't happen but just in case
    }

    double mean = getMean(dataSet);
    // mean /= dataSet.length;

    for (int i = 0; i < dataSet.length; ++i) {
      // iterate over index rather than for-each cuz we must replace data
      dataSet[i] -= mean;
    }

    // test shows this works as in-place method
  }

  /**
   * Linear detrend applied to an array of doubles rather than a list.
   * This operation is not done in-place.
   *
   * @param dataSet The double array to be detrended
   * @return Array of doubles with linear detrend removed
   */
  public static double[] detrend(double[] dataSet) {
    double sumX = 0.0;
    double sumY = 0.0;
    double sumXSqd = 0.0;
    double sumXY = 0.0;

    for (int i = 0; i < dataSet.length; ++i) {
      sumX += i;
      sumXSqd += (double) i * (double) i;
      double value = dataSet[i];
      sumXY += value * i;
      sumY += value;
    }

    // brackets here so you don't get confused thinking this should be
    // algebraic division (in which case we'd just factor out the size term)
    //

    double del = sumXSqd - (sumX * sumX / dataSet.length);

    double slope = sumXY - (sumX * sumY / dataSet.length);
    slope /= del;

    double yOffset = (sumXSqd * sumY) - (sumX * sumXY);
    yOffset /= del * dataSet.length;

    double[] detrended = new double[dataSet.length];

    for (int i = 0; i < dataSet.length; ++i) {
      detrended[i] = dataSet[i] - ((slope * i) + yOffset);
    }

    return detrended;
  }

  public static double[] detrendEnds(double[] data) {
    int lastIdx = data.length - 1;
    double start = data[0];
    double end = data[lastIdx];
    double diff = end - start;
    double delta = diff / data.length;
    double[] dataCopy = data.clone();
    for (int i = 0; i < data.length; ++i) {
      dataCopy[i] -= start + ((delta * i) + 1) * diff / (lastIdx);
    }

    return dataCopy;
  }

  /**
   * Downsamples data by a multiple of passed factor.
   * Result is data.length/factor cells in size, data.length is expected to be even multiple of factor
   * Requires previous use of a low-pass filter to avoid aliasing
   *
   * @param data The timeseries to be downsampled
   * @param factor The factor to decrease the size by
   * @return The downsampled series
   */
  public static double[] downsample(double[] data, int factor) {

    double[] downsamp = new double[data.length / factor];
    for (int i = 0; i < downsamp.length; i++) {
      downsamp[i] = data[i * factor];
    }

    return downsamp;
  }

  /**
   * Implements Euclid's algorithm for finding GCD
   * used to find common divisors to give us upsample
   * and downsample rates by dividing the timeseries intervals
   * by this value
   *
   * @param source Initially, one of two frequencies to calculate
   * @param target Initially, one of two frequencies to calculate
   * @return The GCD of the two frequencies
   */
  public static long euclidGCD(long source, long target) {

    // take remainders until we hit 0
    // which means the divisor is the gcd
    long rem = source % target;
    if (rem == 0) {
      return target;
    }

    return euclidGCD(target, rem);
  }

  /**
   * Extract SNCL data from a SEED data header
   *
   * @param dataHeader found in a seed file
   * @return String containing the SNCL identifier of the data
   */
  private static String extractName(DataHeader dataHeader) {
    StringBuilder fileID = new StringBuilder();
    String station = dataHeader.getStationIdentifier();
    // remove all whitespace from station name
    station = station.replaceAll("\\s+", "");
    fileID.append(dataHeader.getNetworkCode()).append("_");
    fileID.append(station).append("_");
    fileID.append(dataHeader.getLocationIdentifier()).append("_");
    fileID.append(dataHeader.getChannelIdentifier());
    return fileID.toString();
  }

  /**
   * Returns an int representing the number of bytes in a record for
   * a miniSEED file
   *
   * @param filename Full path to the file to be read in
   * @return number of bytes of a record (i.e., 512, 40096)
   * @throws FileNotFoundException If file does not exist
   */
  private static int getByteSize(String filename) throws IOException, SeedFormatException {
    try (DataInputStream dis = new DataInputStream(new FileInputStream(filename))) {
      while (true) {
        SeedRecord sr = SeedRecord.read(dis, 4096);

        Blockette[] blockettes = sr.getBlockettes();

        for (Blockette blockette : blockettes) {
          if (blockette.getType() == 1000) {
            Blockette1000 b1000 = (Blockette1000) blockette;
            return b1000.getDataRecordLength();
          }
        }
      }
    }
  }

  /**
   * Used to quickly get the first data in a file. This is useful if loading in
   * data from a file that is known to not be multiplexed (i.e., containing
   * only the data from a single channel).
   *
   * @param filename Filename of miniSEED data to load in
   * @return Datablock representing the data inside the miniSEED
   */
  public static DataBlock getFirstTimeSeries(String filename)
      throws SeedFormatException, CodecException, IOException {
    String filter = getMplexNameList(filename).get(0);
    return getTimeSeries(filename, filter);
  }

  /**
   * Used to quickly get the data in a list of files, where the first
   * file's first data is the SNCL to filter on. This is useful if loading in
   * data from files that are known to not be multiplexed (i.e., containing
   * only the data from a single channel). Used mainly to load calibration spanning multiple days.
   *
   * @param filenames Filenames of miniSEED data to load in
   * @return Datablock representing the data inside the miniSEEDs
   */
  public static DataBlock getFirstTimeSeries(String[] filenames)
      throws SeedFormatException, CodecException, IOException {
    String filter = getMplexNameList(filenames[0]).get(0);
    return getTimeSeries(filenames, filter);
  }

  /**
   * Return the calculation of the arithmetic mean (using a recursive definition for stability)
   *
   * @param dataSet Range of data to get the mean value from
   * @return the arithmetic mean
   */
  public static double getMean(double[] dataSet) {
    double mean = 0.0;
    double inc = 1;

    for (double data : dataSet) {
      mean = mean + ((data - mean) / inc);
      ++inc;
    }
    return mean;
  }

  /**
   * Get the set of available data series in a multiplexed miniseed file.
   * Because the list is derived from the set, the result of this function
   * should have no duplicate entries.
   * Because Java Sets do not specify ordering, conversion to list (which does) allows for indexing
   * operations to be performed on the returned data more easily.
   *
   * @param filename Name of file to read in
   * @return List of strings corresponding to metadata of each data series
   * in the given miniseed file
   */
  public static List<String> getMplexNameList(String filename)
      throws SeedFormatException, IOException {
    return new ArrayList<>(getMplexNameSet(filename));
  }

  /**
   * Returns list of SNCL (station, network, channel, location) data for
   * a multiplexed miniseed file as a set of strings
   *
   * @param filename miniseed file to be read in
   * @return set of all (unique) SNCL strings
   */
  public static Set<String> getMplexNameSet(String filename)
      throws SeedFormatException, IOException {
    Set<String> dataNames = new HashSet<>();

    int byteSize;
    byteSize = getByteSize(filename);

    DataInputStream dis;

    dis = new DataInputStream(new FileInputStream(filename));

    while (true) {
      // read in the current record, then get the SNCL data
      try {
        SeedRecord sr = SeedRecord.read(dis, byteSize);

        if (sr instanceof DataRecord) {
          DataRecord dr = (DataRecord) sr;
          DataHeader dh = dr.getHeader();

          String fileID = extractName(dh);
          dataNames.add(fileID);
        }

      } catch (EOFException e) {
        // just break out of loop, this means we reached the file end
        break;
      }

    } // end loop until EOF exception

    return dataNames;
  }

  /**
   * Reads in the time series data from a miniSEED file and produces it as a
   * list of doubles according to a given filter (to handle multiplexed data).
   * This is packaged into a data structure that also includes the file's
   * metadata (station, channel, etc.) and the start time and period between
   * samples.
   * Some of this code is based on the miniseed to float array example given
   * in the repository for the included seisFile miniSEED parser library;
   * see the src/.../examples folder under
   * https://github.com/crotwell/seisFile/ for more
   *
   * @param filename The full path to the file to be loaded in
   * @param filter Specifies which data to load in, for multiplexed files
   * @return A structure containing the time series and metadata for the file
   * @throws FileNotFoundException If file cannot be read in
   */
  public static DataBlock getTimeSeries(String filename, String filter)
      throws IOException, SeedFormatException, CodecException {
    Pair<Long, Map<Long, double[]>> intervalSeriesMapPair =
        getTimeSeriesMap(filename, filter);
    return mapToTimeSeries(intervalSeriesMapPair, filter);
  }

  /**
   * Reads in the time series data from miniSEED files and concatenates it as a
   * list of doubles according to a given filter (to handle multiplexed data).
   * This is packaged into a data structure that also includes the file's
   * metadata (station, channel, etc.) and the start time and period between
   * samples.
   * This method is used to handle data that may span multiple days (such as some calibrations).
   * Some of this code is based on the miniseed to float array example given
   * in the repository for the included seisFile miniSEED parser library;
   * see the src/.../examples folder under
   * https://github.com/crotwell/seisFile/ for more
   *
   * @param filenames Each entry is full path of each file to be loaded in
   * @param filter Specifies which data to load in, for multiplexed files
   * @return A structure containing the time series and metadata for the file
   * @throws FileNotFoundException If file cannot be read in
   */
  private static DataBlock getTimeSeries(String[] filenames, String filter)
      throws IOException, SeedFormatException, CodecException {

    Pair<Long, Map<Long, double[]>> intervalSeriesMapPair =
        getTimeSeriesMap(filenames, filter);
    return mapToTimeSeries(intervalSeriesMapPair, filter);

  }

  /**
   * This code is taken from an example in Crotwell's seisFile example code for reading in data
   * from the FDSN in the standard SEED format, which is then read in like SEED data is.
   * @param network Network code (2-character)
   * @param station Station code (4- or 5-character expected)
   * @param location Location code (2-digit numeric or empty expected)
   * @param channel Channel code (3-character alphanumeric expected)
   * @param startTime UTC-zoned date-time representing start time of data from query
   * @param endTime UTC-zoned date-time representing end time of data from query
   * @return DataBlock
   * @throws SeisFileException If there was a problem executing the query
   * @throws IOException If the queried data cannot be parsed as expected
   * @throws CodecException If the data produced in the query has a compression error
   */
  public static DataBlock getTimeSeriesFromFDSNQuery(String network, String station,
      String location, String channel, Date startTime, Date endTime)
      throws SeisFileException, IOException, CodecException {
    // see https://github.com/crotwell/seisFile/blob/master/src/example/.../FDSNDataSelect.java
    FDSNDataSelectQueryParams params = new FDSNDataSelectQueryParams();
    params.setStartTime(startTime).setEndTime(endTime).appendToNetwork(network)
        .appendToStation(station).appendToLocation(location).appendToChannel(channel);
    FDSNDataSelectQuerier querier = new FDSNDataSelectQuerier(params);
    DataRecordIterator iterator = querier.getDataRecordIterator();

    List<DataRecord> dataRecordList = new LinkedList<>();
    while (iterator.hasNext()) {
      dataRecordList.add(iterator.next());
    }

    String filter = network + "_" + station + "_" + location + "_" + channel;

    return mapToTimeSeries(getTimeSeriesMap(dataRecordList), filter);
  }


  /**
   * Extract data from records in a miniseed file and return them as a map
   * of sampled data points at various times
   *
   * @param filename Name of miniseed file to read in
   * @param filter SNCL data of relevant channel to get data from
   * @return Paired value, first entry of which is the interval between points
   * given as a long and second of which is a map from sample times to data
   * points from each given time value in the miniseed records
   * @throws FileNotFoundException if given file from filename cannot be read
   */
  public static Pair<Long, Map<Long, double[]>>
  getTimeSeriesMap(String filename, String filter)
      throws IOException, SeedFormatException, CodecException {
    return getTimeSeriesMap(new String[]{filename}, filter);
  }

  /**
   * Read in multiple miniseed files and concatenate data as long as the data
   * refers to the SNCL data according to the provided filter. If a file does
   * not have any data matching the filter, then its contents will not be added.
   * This is useful for concatenating data automatically from day-crossing
   * calibration data.
   *
   * @param filenames List of miniseed data
   * @param filter SNCL data of relevant channel to get data from
   * @return Paired value, first entry of which is the interval between points
   * given as a long and second of which is a map from sample times to data
   * points from each given time value in the miniseed records
   */
  private static Pair<Long, Map<Long, double[]>>
  getTimeSeriesMap(String[] filenames, String filter)
      throws IOException, SeedFormatException,
      CodecException {

    List<DataRecord> dataRecordList = new LinkedList<>();

    for (String filename : filenames) {
      int byteSize = getByteSize(filename);

      try (DataInputStream dis = new DataInputStream(new FileInputStream(filename))) {
        while (true) {
          try {
            SeedRecord sr = SeedRecord.read(dis, byteSize);
            if (sr instanceof DataRecord) {
              DataRecord dr = (DataRecord) sr;
              DataHeader dh = dr.getHeader();
              String seriesID = extractName(dh);

              if (seriesID.equals(filter)) {
                dataRecordList.add((DataRecord) sr);
              }
            }
          } catch(EOFException e){
              break;
          } // end infinite while loop (read until EOF)
        }
      }
    }

    return getTimeSeriesMap(dataRecordList);
  }


  private static Pair<Long, Map<Long, double[]>>
  getTimeSeriesMap(List<DataRecord> dataRecords)
      throws SeedFormatException, CodecException {

    long interval = 0L;
    Map<Long, double[]> timeListMap = new LinkedHashMap<>();

    for (DataRecord dr : dataRecords) {
      DataHeader dh = dr.getHeader();

      byte af = dh.getActivityFlags();
      byte correctionFlag = 0b00000010; // is there a time correction?
      int correction = 0;
      if ((af & correctionFlag) != 0) {
        correction = dh.getTimeCorrection();
      }
      if (correction > 0) {
        System.out.println("Time correction? " + correction);
      }
      Btime bt = dh.getStartBtime();

      // convert Btime to milliseconds
      long start = bt.convertToCalendar().getTimeInMillis();
      start += correction / 10; // correction in tenths of millis
      //start = (start * 10) + bt.getTenthMilli();

      int fact = dh.getSampleRateFactor();
      int mult = dh.getSampleRateMultiplier();

      // we can assume interval is consistent through a file
      if (fact > 0 && mult > 0) {
        interval = ONE_HZ_INTERVAL / (fact * mult);
      } else if (fact > 0 && mult < 0) {
        interval = Math.abs((ONE_HZ_INTERVAL * mult) / fact);
      } else if (fact < 0 && mult > 0) {
        interval = Math.abs((ONE_HZ_INTERVAL * fact) / mult);
      } else {
        interval = ONE_HZ_INTERVAL * fact * mult;
      }

      timeListMap.put(start, dr.decompress().getAsDouble());
    }

    return new Pair<>(interval, timeListMap);
  }

  /**
   * Convert map from getTimeSeriesMap to a datablock format. This process
   * removes the DC offset and fills in any data gaps with zeros.
   *
   * @param data Pair, first value is long representing interval and the
   * second value is a map from longs representing time of a data sample
   * to a numeric type representing the recorded value at that time
   * @param filter String with SNCL (station, network, channel, location) data
   * @return DataBlock with the given timeseries and metadata
   */
  private static DataBlock
  mapToTimeSeries(Pair<Long, Map<Long, double[]>> data, String filter) {

    long interval = data.getFirst();
    Map<Long, double[]> timeMap = data.getSecond();
    DataBlock db;

    db = new DataBlock(timeMap, interval, filter);
    return db;

  }

  /**
   * Used to determine polarity of data in step calibrations in particular. If the first step
   * response goes negative, the wiring of the sensor has been done effectively backwards.
   * That is, step expects first response to be positive, and so this informs sensor that the
   * input should be multiplied by -1 before calibration params are solved for.
   */
  public static boolean needsSignFlip(double[] data) {
    double min = data[0]; // most negative value in data
    double max = data[0]; // most positive value in data
    int minIdx = 0; // where in list min occurs
    int maxIdx = 0; // where in list max occurs

    for (int i = 1; i < data.length; ++i) {
      if (data[i] > max) {
        max = data[i];
        maxIdx = i;
      } else if (data[i] < min) {
        min = data[i];
        minIdx = i;
      }
    }

    // if the negative peak is first, we need the flip
    return minIdx < maxIdx;
  }

  /**
   * Normalize result according to value of absolute maximum of the data.
   * This is intended to replicate the normalization behavior of Obspy.
   *
   * @param data Time series data to be normalized
   * @return The data normalized by its maximum absolute value
   */
  public static double[] normalize(double[] data) {
    double absMax = Math.abs(data[0]); // initialize with first value in array
    // first get the absolute max
    for (double point : data) {
      absMax = Math.max(Math.abs(point), absMax);
    }

    if (absMax == 0) {
      return data.clone();
    }

    double[] normData = new double[data.length];
    // now scale the data accordingly
    for (int i = 0; i < data.length; ++i) {
      normData[i] = data[i] / absMax;
    }
    return normData;
  }

  /**
   * Rotates a north and east (known orthogonal) set of data and produces a new
   * DataBlock along the north axis in the rotated coordinate system from
   * the given angle, clockwise (y' = y cos theta - x sin theta)
   *
   * @param north Sensor assumed to point north
   * @param east Sensor assumed to point east, orthogonal to north sensor
   * @param angle Angle to rotate the data along
   * @return New DataBlock whose time series is the rotated data along the
   * given angle, facing north
   */
  public static DataBlock rotate(DataBlock north, DataBlock east, double angle) {
    long start = north.getStartTime();
    long interval = north.getInterval();
    String name = north.getName();

    double[] northData = north.getData();
    double[] eastData = east.getData();

    return new DataBlock(rotate(northData, eastData, angle), interval, name, start);
  }

  /**
   * Rotates a north and east (known orthogonal) timeseries and produces a new
   * timeseries along the north axis in the rotated coordinate system from
   * the given angle, clockwise (y' = y cos theta - x sin theta). Though the
   * data given as input are not necessarily from a north-facing and
   * east-facing sensor, they are presumed to be orthogonal to each other.
   *
   * @param northData Timeseries data expected to point north
   * @param eastData Timeseries assumed to point east,
   * orthogonal to north sensor
   * @param angle Angle to rotate the data along
   * @return New timeseries data rotated data along the
   * given angle, facing north
   */
  public static double[] rotate(double[] northData, double[] eastData, double angle) {
    double[] rotatedData = new double[northData.length];

    // clockwise rotation matrix!! That's why things are so screwy
    double sinTheta = Math.sin(angle);
    double cosTheta = Math.cos(angle);

    for (int i = 0; i < northData.length; ++i) {
      rotatedData[i] =
          northData[i] * cosTheta -
              eastData[i] * sinTheta;
    }

    return rotatedData;
  }

  /**
   * Rotates a north and east (known orthogonal) set of data and produces a new
   * DataBlock along the east axis in the rotated coordinate system from
   * the given angle, clockwise (x' = x cos theta + y sin theta)
   *
   * @param north Sensor assumed to point north
   * @param east Sensor assumed to point east, orthogonal to north sensor
   * @param angle Angle to rotate the data along
   * @return New DataBlock whose time series is the rotated data along the
   * given angle, facing east.
   */
  public static DataBlock rotateX(DataBlock north, DataBlock east, double angle) {
    long start = east.getStartTime();
    long interval = east.getInterval();
    String name = east.getName();

    double[] northData = north.getData();
    double[] eastData = east.getData();

    return new DataBlock(rotateX(northData, eastData, angle), interval, name, start);
  }

  /**
   * Rotates a north and east (known orthogonal) timeseries and produces a new
   * timeseries along the east axis in the rotated coordinate system from
   * the given angle, clockwise (x' = x cos theta + y sin theta). Though the
   * data given as input are not necessarily from a north-facing and
   * east-facing sensor, they are presumed to be orthogonal to each other.
   *
   * @param northData Timeseries data expected to point north
   * @param eastData Timeseries assumed to point east,
   * orthogonal to north sensor
   * @param angle Angle to rotate the data along
   * @return New timeseries data rotated data along the
   * given angle, facing east.
   */
  private static double[] rotateX(double[] northData, double[] eastData, double angle) {
    double[] rotatedData = new double[northData.length];

    double sinTheta = Math.sin(angle);
    double cosTheta = Math.cos(angle);

    for (int i = 0; i < northData.length; ++i) {
      rotatedData[i] =
          eastData[i] * cosTheta +
              northData[i] * sinTheta;
    }

    return rotatedData;
  }

  /**
   * Upsamples data by a multiple of passed factor, placing zeros
   * between each data point. Result is data.length*factor cells in size.
   * Requires use of a low-pass filter to remove discontinuities.
   *
   * @param data The timeseries to be upsampled
   * @param factor The factor to increase the size by
   * @return The upsampled series
   */
  private static double[] upsample(double[] data, int factor) {

    int newLength = data.length * factor;

    double[] upsamp = new double[newLength];

    for (int i = 0; i < data.length; ++i) {
      upsamp[i * factor] = data[i]; // index, element
    }

    return upsamp;
  }

}

