package asl.sensor.utils;

import java.io.DataInputStream;
import java.io.EOFException;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math3.util.Pair;

import asl.sensor.input.DataBlock;
import edu.iris.dmc.seedcodec.B1000Types;
import edu.iris.dmc.seedcodec.CodecException;
import edu.iris.dmc.seedcodec.DecompressedData;
import edu.iris.dmc.seedcodec.UnsupportedCompressionType;
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
 * @author akearns
 *
 */
public class TimeSeriesUtils {

  /**
   * Interval for data that has been sampled at 1 Hz in microseconds
   */
  public final static long ONE_HZ_INTERVAL = 1000000L;
  /**
   * Sample rate of a 1 Hz sample, in Hz, as a double (that is, 1.0)
   */
  public final static double ONE_HZ = 1.0;

  /**
   * Initial driver for the decimation utility
   * which takes a timeseries of unknown rate and
   * runs downsampling to convert it to a target
   * frequency of a 1Hz interval.
   * @param data The timeseries to be decimated
   * @param src The source frequency as interval between samples (microseconds)
   * @return A timeseries decimated to the correct frequency
   */
  public static List<Number> decimate(List<Number> data, long src, long tgt){

    // long tgt = ONE_HZ_INTERVAL; // target frequency
    // a sample lower than 1Hz frq has longer time between samples
    // since it's an inverse relationship and all
    if(src >= tgt){
      // if data is too low-frequency to decimate, do nothing
      return data;
    }

    // find what the change in size is going to be
    long gcd = euclidGCD(src, tgt);
    // conversion up- and down-factors
    // (upsample by target, downsample by source)
    // cast is valid because any discrete interval
    // from 1Hz and up is already expressable
    // as an int
    int upf = (int)(src/gcd);
    int dnf = (int)(tgt/gcd);

    double higherFreq = (1. / src) * upf * ONE_HZ_INTERVAL;
    double lowerFreq = (1. / tgt) * ONE_HZ_INTERVAL / 2; 
    // nyquist rate of downsampled data

    // one valid sample rate for data is 2.5Hz
    // with 1Hz that comes out as a ratio of 5/2, which won't
    // downsample neatly in some cases so we would first upsample,
    // filter out any noise terms, then downsample
    List<Number> upped = upsample(data,upf);
    List<Number> lpfed = lowPassFilter(upped, higherFreq, lowerFreq);
    List<Number> down = downsample(lpfed,dnf);

    return down;

  }

  /**
   * Remove the mean from some data
   * @param dataSet timeseries array, double values
   * @return same array but with the mean removed
   */
  public static double[] demean(double[] dataSet) {
    
    List<Number> dataToProcess = new ArrayList<Number>();
    
    for (double number : dataSet) {
      dataToProcess.add(number);
    }
    
    TimeSeriesUtils.demeanInPlace(dataToProcess);
    
    double[] out = new double[dataSet.length];
    for (int i = 0; i < out.length; ++i) {
      out[i] = dataToProcess.get(i).doubleValue();
    }
    return out;
    
  }

  /**
   * Remove mean (constant value) from a dataset and include
   * @param dataSet
   * @return timeseries as numeric list with previous mean subtracted
   */
  public static List<Number> demean(List<Number> dataSet) {
    List<Number> dataOut = new ArrayList<Number>(dataSet);
    TimeSeriesUtils.demeanInPlace(dataOut);
    return dataOut;
  }

  /**
   * In-place subtraction of mean from each point in an incoming data set.
   * This is a necessary step in calculating the power-spectral density.
   * @param dataSet The data to have the mean removed from.
   */
  public static void demeanInPlace(List<Number> dataSet) {
    
    // I'm always getting the demeaning tasks, huh?
    
    if(dataSet.size() == 0) return; // shouldn't happen but just in case
    
    double mean = 0.0;
    
    for(Number data : dataSet) {
      mean += data.doubleValue();
    }
    
    mean /= dataSet.size();
    
    for(int i = 0; i < dataSet.size(); ++i) {
      // iterate over index rather than for-each cuz we must replace data
      // ugly syntax because numeric data types are immutable
      dataSet.set(i, dataSet.get(i).doubleValue() - mean);
    }
    
    // test shows this works as in-place method
  }

  /**
   * Linear detrend applied to an array of doubles rather than a list.
   * This operation is not done in-place.
   * @param dataSet The double array to be detrended
   * @return Array of doubles with linear detrend removed
   */
  public static double[] detrend (double[] dataSet) {
    double sumX = 0.0;
    double sumY = 0.0;
    double sumXSqd = 0.0;
    double sumXY = 0.0;
    
    for (int i = 0; i < dataSet.length; ++i) {
      sumX += (double) i;
      sumXSqd += (double) i * (double) i;
      double value = dataSet[i];
      sumXY += value * (double) i;
      sumY += value;
    }
    
    // brackets here so you don't get confused thinking this should be
    // algebraic division (in which case we'd just factor out the size term)
    // 
    
    double del = sumXSqd - ( sumX * sumX / dataSet.length );
    
    double slope = sumXY - ( sumX * sumY / dataSet.length );
    slope /= del;
    
    double yOffset = (sumXSqd * sumY) - (sumX * sumXY);
    yOffset /= del * dataSet.length;
    
    double[] detrended = new double[dataSet.length];
    
    for (int i = 0; i < dataSet.length; ++i) {
      detrended[i] = dataSet[i] - ( (slope * i) + yOffset);
    }
    
    return detrended;
  }

  /**
   * In-place subtraction of trend from each point in an incoming data set.
   * This is a necessary step in calculating the power-spectral density.
   * @param dataSet The data to have the trend removed from.
   */
  public static void detrend(List<Number> dataSet) {
    
    double sumX = 0.0;
    double sumY = 0.0;
    double sumXSqd = 0.0;
    double sumXY = 0.0;
    
    for (int i = 0; i < dataSet.size(); ++i) {
      sumX += (double) i;
      sumXSqd += (double) i * (double) i;
      double value = dataSet.get(i).doubleValue();
      sumXY += value * (double) i;
      sumY += value;
    }
    
    // brackets here so you don't get confused thinking this should be
    // algebraic division (in which case we'd just factor out the size term)
    // 
    
    double del = sumXSqd - ( sumX * sumX / dataSet.size() );
    
    double slope = sumXY - ( sumX * sumY / dataSet.size() );
    slope /= del;
    
    double yOffset = (sumXSqd * sumY) - (sumX * sumXY);
    yOffset /= del * dataSet.size();
    
    for (int i = 0; i < dataSet.size(); ++i) {
      dataSet.set(i, dataSet.get(i).doubleValue() - ( (slope * i) + yOffset) );
    }
    
  }

  /**
   * Downsamples data by a multiple of passed factor. Result is
   * data.length/factor cells in size
   * Requires previous use of a low-pass filter to avoid aliasing
   * @param data The timeseries to be downsampled
   * @param factor The factor to decrease the size by
   * @return The downsampled series
   */
  public static List<Number> downsample(List<Number> data, int factor){

    List<Number> downsamp = Arrays.asList(new Number[data.size()/factor]);
    for(int i=0; i < downsamp.size(); i++){
      downsamp.set( i, data.get(i*factor) ); 
    }

    return downsamp;
  }
  
  /**
   * Implements Euclid's algorithm for finding GCD
   * used to find common divisors to give us upsample
   * and downsample rates by dividing the timeseries intervals
   * by this value
   * @param src Initially, one of two frequencies to calculate
   * @param tgt Initially, one of two frequencies to calculate
   * @return The GCD of the two frequencies
   */
  public static long euclidGCD(long src,long tgt){

    // take remainders until we hit 0
    // which means the divisor is the gcd
    long rem = src % tgt;
    if(rem == 0){
      return tgt;
    }

    return euclidGCD(tgt, rem);
  }
  
  /**
   * Extract SNCL data from a SEED data header
   * @param dh found in a seed file
   * @return String containing the SNCL identifier of the data
   */
  private static String extractName(DataHeader dh) {
    StringBuilder fileID = new StringBuilder();
    String station = dh.getStationIdentifier();
    // remove all whitespace from station name
    station = station.replaceAll("\\s+","");;
    fileID.append(dh.getNetworkCode() + "_");
    fileID.append(station + "_");
    fileID.append(dh.getLocationIdentifier() + "_");
    fileID.append(dh.getChannelIdentifier());
    return fileID.toString();
  }
  
  /**
   * Returns an int representing the number of bytes in a record for
   * a miniSEED file
   * @param filename Full path to the file to be read in
   * @return number of bytes of a record (i.e., 512, 40096)
   * @throws FileNotFoundException If file does not exist
   */
  public static int getByteSize(String filename) throws FileNotFoundException {

    DataInputStream dis; // used to read in input to get b1000
    int byteSize;
    try {
      dis = new DataInputStream( new FileInputStream(filename) );

      while (true) {

        try {
          SeedRecord sr = SeedRecord.read(dis, 4096);

          Blockette[] blockettes = sr.getBlockettes();

          for (Blockette blockette : blockettes) {
            if ( blockette.getType() == 1000 ) {
              Blockette1000 b1000 = (Blockette1000) blockette;
              byteSize = b1000.getDataRecordLength(); // expect either 9 or 12
              return byteSize;
            }
          } // end of loop over blockettes

        } catch (SeedFormatException e) {
          e.printStackTrace();
        } catch (IOException e) {
          e.printStackTrace();
        } // end of try-catch blocks for parsing an individual record

      } // end of while loop for gotByteSize

    } catch (FileNotFoundException e) {
      throw e;
    } // end of try block for creating DataInputStream

  }

  /**
   * Return the locations of gaps in loaded data
   * @param data Pair, first value is long representing interval and the
   * second value is a map from longs representing time of a data sample
   * to a numeric type representing the recorded value at that time
   * @return Paired values representing the start and end points of time ranges
   * where no data was available in the map
   */
  public static 
  List<Pair<Long, Long>> getGaps(Pair<Long, Map<Long, Number>> data) {
    List<Pair<Long, Long>> gapList = new ArrayList<Pair<Long, Long>>();
    
    long interval = data.getFirst();
    List<Long> times = new ArrayList<Long>( data.getSecond().keySet() );
    Collections.sort(times);

    for (int i = 0; i < times.size(); ++i) {
      long timeNow = times.get(i);
      if ( (i + 1) < times.size() ) {
        long timeNext = times.get(i + 1);
        // is there a discrepancy, and is it big enough for us to care?
        if (timeNext - timeNow != interval) {
          
          if (timeNext - timeNow < interval / 2) {
            System.out.println("SKIPPING POSSIBLE REPEATED VALUE");
            continue;
          }
          
          // it's a gap if it's more than twice the interval size
          if (timeNext - timeNow  > interval * 2) {
            gapList.add( new Pair<Long, Long>(timeNow, timeNext) );
          }

          // long gap = timeNext - timeNow;
          // System.out.println("FOUND GAP: " + timeNow + ", " + timeNext);
          // System.out.println("(Itvl: " + interval + "; gap: " + gap + ")");
          while (timeNext - timeNow > interval * 2) {
            timeNow += interval;
          }
        }
      }
    }
    
    return gapList;
  }
  
  /**
   * Get the set of available data series in a multiplexed miniseed file.
   * Because the list is derived from the set, the result of this function
   * should have no duplicate entries.
   * Because Java Sets do not specify ordering, this allows for indexing
   * operations to be performed on the returned data more easily.
   * @param filename Name of file to read in
   * @return List of strings corresponding to metadata of each data series
   * in the given miniseed file
   * @throws FileNotFoundException If given file from filename cannot be read
   */
  public static List<String> getMplexNameList(String filename)
      throws FileNotFoundException {
    return new ArrayList<String>( getMplexNameSet(filename) );
  }
  
  /**
   * Returns list of SNCL (station, network, channel, location) data for
   * a multiplexed miniseed file as a set of strings
   * @param filename miniseed file to be read in
   * @return set of all (unique) SNCL strings
   * @throws FileNotFoundException if file cannot be read
   */
  public static Set<String> getMplexNameSet(String filename) 
      throws FileNotFoundException {
    Set<String> dataNames = new HashSet<String>();

    int byteSize;
    try {
      byteSize = getByteSize(filename);
    } catch (FileNotFoundException e1) {
      throw e1;
    }

    DataInputStream dis;

    try {
      dis = new DataInputStream( new FileInputStream(filename) );

      while (true) {
        // read in the current record, then get the SNCL data
        try {
          SeedRecord sr = SeedRecord.read(dis, byteSize);

          if (sr instanceof DataRecord) {
            DataRecord dr = (DataRecord)sr;
            DataHeader dh = dr.getHeader();

            String fileID = extractName(dh);
            dataNames.add(fileID);
          }

        } catch (EOFException e) {
          // just break out of loop, this means we reached the file end
          break;
        }

      } // end loop until EOF exception

    } catch (FileNotFoundException e) {
      // Auto-generated catch block
      e.printStackTrace();
    } catch (IOException e) {
      // Auto-generated catch block
      e.printStackTrace();
    } catch (SeedFormatException e) {
      // Auto-generated catch block
      e.printStackTrace();
    }

    return dataNames;
  }
  
  /**
   * Get all contiguous ranges of data from a read-in time-value map
   * @param data Pair, first value is long representing interval and the
   * second value is a map from longs representing time of a data sample
   * to a numeric type representing the recorded value at that time
   * @return List of paired longs, representing the start and end times
   * in microseconds of each contiguous range in the map
   */
  public static List<Pair<Long, Long>> 
  getRanges(Pair<Long, Map<Long, Number>> data) {
    
    List<Pair<Long, Long>> gapList = getGaps(data);
    
    List<Long> times = new ArrayList<Long>( data.getSecond().keySet() );
    Collections.sort(times);
    
    List<Pair<Long, Long>> contiguous = new ArrayList<Pair<Long, Long>>();
    
    // PairComparator sorts the gaps by start time
    Collections.sort(gapList, new StartTimeComparator() );
    long start = times.get(0);
    long end = times.get( times.size() - 1 );
    
    if (gapList.size() == 0) {
      contiguous.add( new Pair<Long, Long>(start, end) );
      return contiguous;
    }
    
    contiguous.add( new Pair<Long, Long>(start, gapList.get(0).getFirst() ) );
    
    for (int i = 0; i < gapList.size() - 1; ++i) {
      contiguous.add( 
          new Pair<Long, Long>( 
              gapList.get(i).getSecond(), 
              gapList.get(i + 1).getFirst() ) );
    }
    
    Long last = gapList.get( gapList.size() - 1).getSecond();
    contiguous.add( new Pair<Long, Long>(last, end) );
    Collections.sort(contiguous, new StartTimeComparator() );
    return contiguous;
  }

  /**
   * Get range of contiguous data blocks in a file
   * @param file Seed file to read in
   * @param filter SNCL data representing relevant timeseries
   * @return List of regions of data within which no time gaps exist
   * @throws FileNotFoundException If file does not exist
   */
  public static List<Pair<Long, Long>> getRanges(String file, String filter) 
      throws FileNotFoundException {
    return getRanges( getTimeSeriesMap(file, filter) );
  }

  /**
   * Reads in the time series data from a miniSEED file and produces it as a
   * list of Java numerics, which can be shorts, floats, doubles, or longs,
   * reflecting the format of the data in the file which can be any of these.
   * This is packaged into a data structure that also includes the file's
   * metadata (station, channel, etc.) and the start time and period between
   * samples.
   * Some of this code is based on the miniseed to float array example given
   * in the repository for the included seisFile miniSEED parser library;
   * see the src/.../examples folder under
   * https://github.com/crotwell/seisFile/ for more
   * @param filename The full path to the file to be loaded in
   * @param filter Specifies which data to load in, for multiplexed files
   * @return A structure containing the time series and metadata for the file
   * @throws FileNotFoundException If file cannot be read in
   */
  public static DataBlock getTimeSeries(String filename, String filter)
      throws FileNotFoundException {

    // XYSeries xys = null;
    DataBlock db = null;
    Pair<Long, Map<Long, Number>> intervalSeriesMapPair = 
        getTimeSeriesMap(filename, filter);
    List<Pair<Long, Long>> gaps = getGaps(intervalSeriesMapPair);
    db = mapToTimeSeries(intervalSeriesMapPair, filter);
    db.setGapLocations(gaps);
    return db;

  }

  /**
   * Reads in the time series data from a miniSEED file and produces it as a
   * list of Java numerics, which can be shorts, floats, doubles, or longs,
   * reflecting the format of the data in the file which can be any of these.
   * This is packaged into a data structure that also includes the file's
   * metadata (station, channel, etc.) and the start time and period between
   * samples. The data loaded in is pre-trimmed according a specified range,
   * mainly to be used for loading in contiguous regions of data based on
   * gap locations
   * @param filename The full path to the file to be loaded in
   * @param filter Specifies which data to load in, for multiplexed files
   * @param range Range of time to trim data to before loading
   * @return A structure containing the time series and metadata for the file
   * @throws FileNotFoundException If file cannot be read in
   */
  public static DataBlock 
  getTimeSeries(String filename, String filter, Pair<Long, Long> range) 
      throws FileNotFoundException {
    
    DataBlock db = null;
    Pair<Long, Map<Long, Number>> intervalSeriesMapPair = 
        getTimeSeriesMap(filename, filter);
    List<Pair<Long, Long>> gaps = getGaps(intervalSeriesMapPair);
    db = mapToTimeSeries(intervalSeriesMapPair, filter, range);
    db.setGapLocations(gaps);
    return db;
  }

  /**
   * Extract data from records in a miniseed file and return them as a map
   * of sampled data points at various times
   * @param filename Name of miniseed file to read in
   * @param filter SNCL data of relevant channel to get data from
   * @return Paired value, first entry of which is the interval between points
   * given as a long and second of which is a map from sample times to data 
   * points from each given time value in the miniseed records
   * @throws FileNotFoundException if given file from filename cannot be read
   */
  public static Pair<Long, Map<Long, Number>>
   getTimeSeriesMap(String filename, String filter) 
      throws FileNotFoundException {

    long interval = 0L;
    DataInputStream dis;
    Map<Long, Number> timeMap = new HashMap<Long, Number>();
    
    int byteSize = 512;
    try {
      byteSize = getByteSize(filename);
    } catch (FileNotFoundException e1) {
      throw e1;
    }

    try {
      dis = new DataInputStream(  new FileInputStream(filename) );

      while ( true ) {

        try {
          SeedRecord sr = SeedRecord.read(dis, byteSize);
          if (sr instanceof DataRecord) {
            DataRecord dr = (DataRecord)sr;
            DataHeader dh = dr.getHeader();
            String seriesID = extractName(dh);

            if ( !seriesID.equals(filter) ){
              // System.out.println(seriesID);
              continue; // skip to next seedRecord
            }

            // byte af = dh.getActivityFlags();
            // byte correctionFlag = 0b00000010; // is there a time correction?
            // int correction = 0;
            // if ( (af & correctionFlag) != 0 ) {
            //   correction = dh.getTimeCorrection();
            // }

            Btime bt = dh.getStartBtime();

            // convert Btime to microseconds first as milliseconds
            long start = bt.convertToCalendar().getTimeInMillis();

            // .1 ms = 100 microseconds
            start *= 1000;
            // start += correction;

            int fact = dh.getSampleRateFactor();
            int mult = dh.getSampleRateMultiplier();

            // we can assume interval is consistent through a file
            if( fact > 0 && mult > 0) {
              interval = ONE_HZ_INTERVAL / (fact * mult);
            } else if (fact > 0 && mult < 0) {
              interval = Math.abs( (ONE_HZ_INTERVAL * mult) / fact);
            } else if (fact < 0 && mult > 0) {
              interval = Math.abs( (ONE_HZ_INTERVAL * fact) / mult);
            } else {
              interval = ONE_HZ_INTERVAL * fact * mult;
            }

            DecompressedData decomp = dr.decompress();

            // get the original datatype of the series (loads data faster)
            // otherwise the decompressed data gets converted (cloned) as
            // the other type instead
            int dataType = decomp.getType();
            long timeOfData = start;

            switch (dataType) {
            case B1000Types.INTEGER:
              int[] decomArrayInt = decomp.getAsInt();
              for (int dataPoint : decomArrayInt ) {
                timeMap.put(timeOfData,dataPoint);
                timeOfData += interval;
              }
              break;
            case B1000Types.FLOAT:
              float[] decomArrayFlt = decomp.getAsFloat();
              for (float dataPoint : decomArrayFlt ) {
                timeMap.put(timeOfData,dataPoint);
                timeOfData += interval;
              }
              break;
            case B1000Types.SHORT:
              short[] decomArrayShr = decomp.getAsShort();
              for (short dataPoint : decomArrayShr ) {
                timeMap.put(timeOfData,dataPoint);
                timeOfData += interval;
              }
              break;
            default:
              double[] decomArrayDbl = decomp.getAsDouble();
              for (double dataPoint : decomArrayDbl ) {
                timeMap.put(timeOfData,dataPoint);
                timeOfData += interval;
              }
              break;
            }

          }
        } catch(EOFException e) {
          break;
        }

      } // end infinite while loop (read until EOF)

    } catch (FileNotFoundException e) {
      // Auto-generated catch block
      e.printStackTrace();
    } catch (IOException e) {
      // Auto-generated catch block
      e.printStackTrace();
    } catch (SeedFormatException e) {
      // Auto-generated catch block
      e.printStackTrace();
    } catch (UnsupportedCompressionType e) {
      // Auto-generated catch block
      e.printStackTrace();
    } catch (CodecException e) {
      // Auto-generated catch block
      e.printStackTrace();
    }

    return new Pair<Long, Map<Long, Number>>(interval, timeMap);
  }

  /**
   * Implements low pass band filter
   * @param timeseries  The data to be filtered
   * @param sps         Samples per second
   * @return            The filtered data
   */
  public static List<Number> 
  lowPassFilter(List<Number> timeseries, double sps, double corner)
  {

    double fl = 0.;
    double fh = corner;

    return FFTResult.bandFilter(timeseries, sps, fl, fh);

    /*
    List<Number> timeseriesOut = new ArrayList<Number>();

    for (int i = 0; i < timeseries.size(); ++i) {
      double point = timeseriesFilter.get(i);
      // System.out.println(point);
      timeseriesOut.add(point);
    }

    return timeseriesOut;
     */
  }

  /**
   * Convert map from getTimeSeriesMap to a datablock format. This process
   * removes the DC offset and fills in any data gaps with zeros.
   * @param data Pair, first value is long representing interval and the
   * second value is a map from longs representing time of a data sample
   * to a numeric type representing the recorded value at that time
   * @param filter String with SNCL (station, network, channel, location) data
   * @return DataBlock with the given timeseries and metadata
   */
  public static 
  DataBlock mapToTimeSeries(Pair<Long, Map<Long, Number>> data, String filter) {
    
    // use max range instead of a trim region
    // TODO: likely will need to change range in year 2038
    // or to use data from before 1970
    Pair<Long, Long> range = 
        new Pair<Long, Long>( 0L, Long.MAX_VALUE );
    
    return mapToTimeSeries(data, filter, range);
    
  }

  /**
   * Convert loaded in time-value map to contiguous timeseries over a given
   * range.
   * @param data Pair, first value is long representing interval and the
   * second value is a map from longs representing time of a data sample
   * to a numeric type representing the recorded value at that time
   * @param filter SNCL data taken from miniseed data
   * @param range Range to pre-trim data to
   * @return DataBlock consisting of timeseries data within given range
   */
  public static DataBlock 
  mapToTimeSeries(Pair<Long, Map<Long, Number>> data, String filter, 
      Pair<Long, Long> range) {
    
    long start = Math.min( range.getFirst(), range.getSecond() );
    long end = Math.max( range.getFirst(), range.getSecond() );
    
    long interval = data.getFirst();
    Map<Long, Number> timeMap = data.getSecond();
    DataBlock db;
    
    // get data mean
    double mean = 0.;
    for ( Number sample : timeMap.values() ) {
      mean += sample.doubleValue();
    }
    
    mean /= timeMap.values().size();

    // now we have all the data in a convenient map timestamp -> value
    // which we can then convert into an easy array
    List<Long> times = new ArrayList<Long>( timeMap.keySet() );
    Collections.sort(times);
    
    // get the min value in the set, the start time for the series
    long startTime = times.get(0);
    
    start = Math.max(start, startTime); // data can't start before first sample
    // this is necessary since we use this value as start time of data block
   
    // read in data from the records as long as they exist
    // if no data exists (there's a gap in the record), set value to 0
    // this is done to handle cases where multiplexed files have non-matching
    // gaps and similar issues that previous code was not able to handle

    List<Number> timeList = new ArrayList<Number>();
    
    for (int i = 0; i < times.size(); ++i) {
      long timeNow = times.get(i);
      
      if (timeNow < start) {
        continue;
      } else if (timeNow > end) {
        break;
      }
      
      // do a demean here so that we can add 0-values to empty points
      // without those values affecting the removal of a DC offset later
      timeList.add( timeMap.get(timeNow).doubleValue() );

      if ( (i + 1) < times.size() ) {
        long timeNext = times.get(i + 1);
        // is there a discrepancy, and is it big enough for us to care?
        if (timeNext - timeNow != interval) {
          
          if (timeNext - timeNow < interval / 2) {
            System.out.println("SKIPPING POSSIBLE REPEATED VALUE");
            continue;
          }
          
          if (timeNext - timeNow > interval && 
              timeNext - timeNow < interval * 1.75) {
            System.out.println("CORRECTING SMALL TIME ROUNDING ERROR");
            continue;
          }
          
          // long gap = timeNext - timeNow;
          // System.out.println("FOUND GAP: " + timeNow + ", " + timeNext);
          // System.out.println("(Itvl: " + interval + "; gap: " + gap + ")");
          while (timeNext - timeNow > interval) {
            timeList.add(mean);
            // why use mean value? puts the value in gaps as the center of
            // data -- useful if the data has a large DC offset
            timeNow += interval;
          }
        }
      }

    }

    
    // demean the input to remove DC offset before adding it to the data
    // List<Number> listOut = TimeSeriesUtils.demean( timeList );
    // since we've demeaned the data while adding it in, don't need to do that
    db = new DataBlock(timeList, interval, filter, start);
    return db;
    
  }

  /** 
   * Scales data of an arbitrary range to lie within a [-1, 1] range
   * @param data Timeseries data
   * @return Same data, over the range [-1, 1], linearly scaled
   */
  public static double[] normalize(double[] data) {
    double max = Double.NEGATIVE_INFINITY;
    double min = Double.POSITIVE_INFINITY;

    for (double point : data) {
      if (point < min) {
        min = point;
      }
      if (point > max) {
        max = point;
      }
    }

    for (int i = 0; i < data.length; ++i) {
      // scale to range (0,2) then to (-1, 1)
      data[i] = ( 2 * (data[i] - min) / (max - min) ) - 1;
    }

    return data;

  }

  public static List<Number> normalize(List<Number> data) {
    double max = Double.NEGATIVE_INFINITY;
    double min = Double.POSITIVE_INFINITY;

    for (Number point : data) {
      if (point.doubleValue() < min) {
        min = point.doubleValue();
      }
      if (point.doubleValue() > max) {
        max = point.doubleValue();
      }
    }

    for (int i = 0; i < data.size(); ++i) {
      // scale to range (0,2) then to (-1, 1)
      Double previous = data.get(i).doubleValue();
      data.set(i, 2 * ( (previous - min) / (max-min) ) - 1 );
    }

    return data;

  }

  /**
   * Rotates a north and east (known orthognal) set of data and produces a new
   * DataBlock along the north axis in the rotated coordinate system from
   * the given angle, clockwise (y' = y cos theta - x sin theta)
   * @param north Sensor assumed to point north
   * @param east Sensor assumed to point east, orthogonal to north sensor
   * @param ang Angle to rotate the data along
   * @return New DataBlock whose time series is the rotated data along the
   * given angle, facing north
   */
  public static DataBlock rotate(DataBlock north, DataBlock east, double ang) {
    DataBlock rotated = new DataBlock(north);
    List<Number> northData = north.getData();
    List<Number> eastData = east.getData();
    List<Number> rotatedData = new ArrayList<Number>();

    // clockwise rotation matrix!! That's why things are so screwy
    double sinTheta = Math.sin(ang);
    double cosTheta = Math.cos(ang);

    for (int i = 0; i < northData.size(); ++i) {
      rotatedData.add( 
          northData.get(i).doubleValue() * cosTheta - 
          eastData.get(i).doubleValue() * sinTheta );
    }

    rotated.setData(rotatedData);

    return rotated;
  }

  /**
   * Rotates a north and east (known orthognal) set of data and produces a new
   * DataBlock along the east axis in the rotated coordinate system from
   * the given angle, clockwise (x' = x cos theta + y sin theta)
   * @param north Sensor assumed to point north
   * @param east Sensor assumed to point east, orthogonal to north sensor
   * @param ang Angle to rotate the data along
   * @return New DataBlock whose time series is the rotated data along the
   * given angle, facing east.
   */
  public static DataBlock rotateX(DataBlock north, DataBlock east, double ang) {
    DataBlock rotated = new DataBlock(east); // use east component metadata
    List<Number> northData = north.getData();
    List<Number> eastData = east.getData();
    List<Number> rotatedData = new ArrayList<Number>();

    double sinTheta = Math.sin(ang);
    double cosTheta = Math.cos(ang);

    for (int i = 0; i < northData.size(); ++i) {
      rotatedData.add( 
          eastData.get(i).doubleValue() * cosTheta + 
          northData.get(i).doubleValue() * sinTheta );
    }

    rotated.setData(rotatedData);

    return rotated;
  }

  /**
   * Upsamples data by a multiple of passed factor, placing zeros
   * between each data point. Result is data.length*factor cells in size.
   * Requires use of a low-pass filter to remove discontinuities.
   * @param data The timeseries to be upsampled
   * @param factor The factor to increase the size by
   * @return The upsampled series
   */
  public static List<Number> upsample(List<Number> data, int factor){

    List<Number> upsamp = new ArrayList<Number>();

    for(int i = 0; i < data.size() * factor; ++i) {
      upsamp.add( new Double(0.) );
    }

    for(int i = 0; i < data.size(); ++i){
      upsamp.set( i*factor, data.get(i) ); // index, element
    }

    return upsamp;
  }

}

