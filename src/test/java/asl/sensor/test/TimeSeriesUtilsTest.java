package asl.sensor.test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;
import java.io.BufferedInputStream;
import java.io.DataInput;
import java.io.DataInputStream;
import java.io.EOFException;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.time.Instant;
import java.time.OffsetDateTime;
import java.time.ZoneOffset;
import java.time.format.DateTimeFormatter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import org.junit.Before;
import org.junit.Test;
import asl.sensor.input.DataBlock;
import asl.sensor.utils.FFTResult;
import asl.sensor.utils.TimeSeriesUtils;
import edu.iris.dmc.seedcodec.B1000Types;
import edu.iris.dmc.seedcodec.CodecException;
import edu.iris.dmc.seedcodec.DecompressedData;
import edu.iris.dmc.seedcodec.UnsupportedCompressionType;
import edu.sc.seis.seisFile.mseed.DataRecord;
import edu.sc.seis.seisFile.mseed.SeedFormatException;
import edu.sc.seis.seisFile.mseed.SeedRecord;

public class TimeSeriesUtilsTest {

  public static String folder = TestUtils.DL_DEST_LOCATION + TestUtils.SUBPAGE;

  public String station = "TST5";
  public String location = "00";
  public String channel = "BH0";
  public String fileID = station+"_"+location+"_"+channel+".512.seed";

  @Test
  public void canGetFile() {
    try{
      String filename1 = folder + "blocktrim/" + fileID;
      FileInputStream fis = new FileInputStream(filename1);
      fis.close();
    } catch (Exception e) {
      assertNull(e);
    }
  }

  @Test
  public void canGetMultiplexDataNames() {
    String filename2 = folder + "multiplex/cat.seed";
    Set<String> names;

    try {
      names = TimeSeriesUtils.getMplexNameSet(filename2);
      assertTrue(names.contains("IU_ANMO_00_LH1"));
      assertTrue(names.contains("IU_ANMO_00_LH2"));
      assertTrue(names.contains("IU_ANMO_00_LHZ"));
      assertEquals( names.size(), 3 );
    } catch (FileNotFoundException e) {
      fail();
    }
  }

  @Test
  public void decimationTest() {

    long interval40Hz = (TimeSeriesUtils.ONE_HZ_INTERVAL / 40);
    long interval = TimeSeriesUtils.ONE_HZ_INTERVAL;

    double[] timeSeries = new double[160];

    for (int i = 0; i < 160; ++i) {
      timeSeries[i] = i;
    }

    System.out.println(Arrays.toString(timeSeries));

    // set anti-aliasing frequency to HALF of sampling rate (i.e., nyq rate)
    double[] filtered = FFTResult.lowPassFilter(timeSeries, 40., 0.5);

    System.out.println(Arrays.toString(filtered));

    timeSeries = TimeSeriesUtils.decimate(timeSeries, interval40Hz, interval);

    System.out.println(Arrays.toString(timeSeries));

    assertEquals(timeSeries.length, 4);
    for (int i = 0; i < timeSeries.length; ++i) {
      assertEquals(timeSeries[i], filtered[40 * i], 0.1);
      // is the data nearly halfway between points?
      assertEquals(timeSeries[i], Math.max(0., 40. * i - 20.), 3.0);
    }
  }

  @Test
  public void demeaningTest() {

    // tests that demean does what it says it does and that
    // the results are applied in-place

    double[] numbers = {1,2,3,4,5};

    double[] numList = numbers.clone();
    double[] demeaned = numList.clone();

    TimeSeriesUtils.demeanInPlace(demeaned);

    for (int i = 0; i < numList.length; ++i) {
      assertEquals(demeaned[i], numList[i]-3, 1E-15);
    }

  }

  @Test
  public void detrendAtEndsTest() {
    double[] x = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,
        18, 19, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4,
        3, 2, 1 };

    x = TimeSeriesUtils.detrendEnds(x);

    assertEquals(x[0], 0, 1E-5);
    assertEquals(x[x.length-1], 0, 1E-5);
  }

  @Test
  public void detrendingCycleTest() {

    double[] x = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,
        18, 19, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4,
        3, 2, 1 };

    // List<Number> toDetrend = Arrays.asList(x);

    double[] answer = { -9d, -8d, -7d, -6d, -5d, -4d, -3d, -2d, -1d, 0d, 1d, 2d,
        3d, 4d, 5d, 6d, 7d, 8d, 9d, 10d, 9d, 8d, 7d, 6d, 5d, 4d, 3d, 2d, 1d, 0d,
        -1d, -2d, -3d, -4d, -5d, -6d, -7d, -8d, -9d };


    x = TimeSeriesUtils.detrend(x);

    for (int i = 0; i < x.length; i++) {
      assertEquals( x[i],  answer[i], 0.5);
    }

  }

  @Test
  public void detrendingLinearTest() {

    Number[] x = { 1, 2, 3, 4, 5, 6, 7, 8, 9};

    List<Number> toDetrend = Arrays.asList(x);
    TimeSeriesUtils.detrend(toDetrend);

    for (Number num : toDetrend) {
      assertEquals(num.doubleValue(), 0.0, 0.001);
    }

  }

  @Test
  public void dumbDivisionTest() {
    int div = 12;
    double num = 1.44;
    double res = num / div;
    assertEquals(0.12, res, 1E-10);
  }

  @Test
  public void firstSampleCorrect() {
    String fname = folder + "random_cal_lowfrq/BHZ.512.seed";
    try {
      String data = TimeSeriesUtils.getMplexNameList(fname).get(0);
      DataBlock db = TimeSeriesUtils.getTimeSeries(fname, data);
      long start = db.getStartTime();
      Map<Long, double[]> timeseries = db.getDataMap();
      double[] firstContiguous = timeseries.get(start);
      long sum = 0;
      for (Number n : firstContiguous) {
        sum += n.longValue();
      }
      assertEquals(2902991374L, sum);
      assertEquals(1652432, firstContiguous.length);
      // System.out.println(timeseries.get(start)[0]);

    } catch (FileNotFoundException | SeedFormatException | CodecException e) {
      fail();
      e.printStackTrace();
    }
  }

  @Test
  public void firstSampleCorrect2() {
    String fname = folder + "random_cal_lowfrq/BC0.512.seed";
    try {
      String data = TimeSeriesUtils.getMplexNameList(fname).get(0);
      DataBlock db = TimeSeriesUtils.getTimeSeries(fname, data);
      long start = db.getStartTime();
      Map<Long, double[]> timeseries = db.getDataMap();
      double[] firstContiguous = timeseries.get(start);
      assertEquals(561682, firstContiguous.length);
      int halfLen = firstContiguous.length / 2;
      firstContiguous = Arrays.copyOfRange(firstContiguous, 0, halfLen);
      long sum = 0;
      for (Number n : firstContiguous) {
        sum += n.longValue();
      }
      assertEquals(707752187L, sum);
      OffsetDateTime cCal = OffsetDateTime.ofInstant(db.getStartInstant(), ZoneOffset.UTC);
      String correctDate = "2017.08.02 | 00:00:00.019-Z";
      DateTimeFormatter dtf = DateTimeFormatter.ofPattern("YYYY.MM.dd | HH:mm:ss.SSS-X");
      String inputDate = cCal.format(dtf);
      assertEquals(inputDate, correctDate);
      assertEquals(20.0, db.getSampleRate(), 1E-20);
      // System.out.println(timeseries.get(start)[0]);

    } catch (FileNotFoundException | SeedFormatException | CodecException e) {
      fail();
      e.printStackTrace();
    }
  }

  @Before
  public void getReferencedData() {

    // place in sprockets folder under 'from-sensor-test/[test-name]'

    String refSubfolder = TestUtils.SUBPAGE + "cowi-multitests/";
    String filename = "C100823215422_COWI.LHx";
    String filename2 = "DT000110.LH1";
    try {
      TestUtils.downloadTestData(refSubfolder, filename, refSubfolder, filename);
      TestUtils.downloadTestData(refSubfolder, filename2, refSubfolder, filename2);
    } catch (IOException e) {
      e.printStackTrace();
    }

    refSubfolder = TestUtils.SUBPAGE + "blocktrim/";
    try {
      TestUtils.downloadTestData(refSubfolder, fileID, refSubfolder, fileID);
    } catch (IOException e) {
      e.printStackTrace();
    }

    refSubfolder = TestUtils.SUBPAGE + "kiev-step/";
    filename = "_BC0.512.seed";
    try {
      TestUtils.downloadTestData(refSubfolder, filename, refSubfolder, filename);
    } catch (IOException e) {
      e.printStackTrace();
    }

    refSubfolder = TestUtils.SUBPAGE + "multiplex/";
    filename = "cat.seed";
    try {
      TestUtils.downloadTestData(refSubfolder, filename, refSubfolder, filename);
    } catch (IOException e) {
      e.printStackTrace();
    }

    refSubfolder = TestUtils.SUBPAGE + "random_cal_lowfrq/";
    filename = "BHZ.512.seed";
    filename2 = "BC0.512.seed";
    try {
      TestUtils.downloadTestData(refSubfolder, filename, refSubfolder, filename);
      TestUtils.downloadTestData(refSubfolder, filename2, refSubfolder, filename2);
    } catch (IOException e) {
      e.printStackTrace();
    }

  }

  @Test
  public void inputFileReaderCreatesXYSeries() {
    DataInput dis;
    List<Number> data = new ArrayList<Number>();
    String filename1 = folder + "blocktrim/" + fileID;
    try {
      dis = new DataInputStream(
            new BufferedInputStream(
            new FileInputStream(filename1) ) );
      while ( true ) {

        try {
          // long interval = 0L;
          SeedRecord sr = SeedRecord.read(dis,4096);
          if(sr instanceof DataRecord) {
            DataRecord dr = (DataRecord)sr;
            // DataHeader dh = dr.getHeader();

            DecompressedData decomp = dr.decompress();

            // get the original datatype of the series (loads data faster)
            // otherwise the decompressed data gets converted (cloned) as
            // the other type instead
            int dataType = decomp.getType();

            // This is probably the best way to do this since
            // we have to add each point individually and type convert anyway

            switch (dataType) {
            case B1000Types.INTEGER:
              int[] decomArrayInt = decomp.getAsInt();
              for (int dataPoint : decomArrayInt ) {
                data.add(dataPoint);
              }
              break;
            case B1000Types.FLOAT:
              float[] decomArrayFlt = decomp.getAsFloat();
              for (float dataPoint : decomArrayFlt ) {
                data.add(dataPoint);
              }
              break;
            case B1000Types.SHORT:
              short[] decomArrayShr = decomp.getAsShort();
              for (short dataPoint : decomArrayShr ) {
                data.add(dataPoint);
              }
              break;
            default:
              double[] decomArrayDbl = decomp.getAsDouble();
              for (double dataPoint : decomArrayDbl ) {
                data.add(dataPoint);
              }
              break;
            }
          }
        } catch(EOFException e) {
          break;
        }

      }

      // quickly get the one name in the list
      Set<String> names = TimeSeriesUtils.getMplexNameSet(filename1);
      List<String> nameList = new ArrayList<String>(names);
      System.out.println("DATA BLOCK SIZE: " + data.size());

      DataBlock testAgainst =
          TimeSeriesUtils.getTimeSeries(filename1, nameList.get(0) );
      assertEquals( data.size(), testAgainst.getData().length );

    } catch (FileNotFoundException e) {
      assertNull(e);
    } catch (SeedFormatException e) {
      assertNull(e);
    } catch (IOException e) {
      assertNull(e);
    } catch (UnsupportedCompressionType e) {
      assertNull(e);
    } catch (CodecException e) {
      assertNull(e);
    }
  }

  @Test
  public void seisFileCanParseFile() {
    String filename1 = folder + "blocktrim/" + fileID;
    try {
      DataInput dis = new DataInputStream( new BufferedInputStream(
          new FileInputStream(filename1) ) );
      try{
        while(true) {
          SeedRecord sr = SeedRecord.read(dis,4096);
          if (sr instanceof DataRecord) {
            DataRecord dr = (DataRecord)sr;

            String loc = dr.getHeader().getLocationIdentifier();
            assertTrue( loc.equals(location) );
            String stat = dr.getHeader().getStationIdentifier().trim();
            assertTrue( stat.equals(station) );

            String chan = dr.getHeader().getChannelIdentifier();
            assertTrue( chan.equals(channel) );
          }
        }
      } catch (EOFException e) {
        assertNotNull(e); // I haaates it! I haaaaaaaaaates it!
      } catch (SeedFormatException e) {
        assertNull(e);
      } catch (IOException e) {
        assertNull(e);
      }

    } catch (FileNotFoundException e) {
      assertNull(e);
    }

  }

  @Test
  public void seisFileGivesCorrectSampleRateAndInterval() {
    DataInput dis;
    String filename1 = folder + "blocktrim/" + fileID;
    try {
      while (true) {
        dis = new DataInputStream( new BufferedInputStream(
            new FileInputStream(filename1) ) );
        SeedRecord sr = SeedRecord.read(dis,4096);
        if(sr instanceof DataRecord) {
          DataRecord dr = (DataRecord)sr;

          int fact = dr.getHeader().getSampleRateFactor();
          int mult = dr.getHeader().getSampleRateMultiplier();

          //System.out.println(fact+","+mult);

          double rate = dr.getHeader().getSampleRate();
          assertTrue((double)fact/mult == rate);

          // checking the correct values for the intervals

          double multOf1Hz = rate/TimeSeriesUtils.ONE_HZ;
          long inverse = TimeSeriesUtils.ONE_HZ_INTERVAL/(long)multOf1Hz;

          long interval = TimeSeriesUtils.ONE_HZ_INTERVAL*mult/fact;

          assertEquals( inverse, interval);
          // System.out.println(interval);

          break;

        }
      }
    } catch (FileNotFoundException e) {
      assertNull(e); // only reading one record;
    } catch (SeedFormatException e) {
      assertNull(e);
    } catch (IOException e) {
      assertNull(e);
    }
  }

  @Test
  public void testCOWIDemean() {
    String filename = folder + "cowi-multitests/C100823215422_COWI.LHx";
    String dataname = "US_COWI_  _LHN";
    DataBlock db;
    try {
      db = TimeSeriesUtils.getTimeSeries(filename, dataname);
      String startString = "2010-236T02:00:00.0";
      String endString = "2010-236T13:00:00.0";
      long st = TestUtils.timeStringToEpochMilli(startString);
      long ed = TestUtils.timeStringToEpochMilli(endString);
      db.trim(st, ed);
      double[] data = db.getData();
      data = TimeSeriesUtils.demean(data);
      double mean = 0.;
      for (double point : data) {
        mean += point;
      }
      mean /= data.length;
      assertEquals(0., mean, 1E-10);
    } catch (FileNotFoundException | SeedFormatException | CodecException e) {
      e.printStackTrace();
      fail();
    }

  }

  @Test
  public final void testDemean1to9() throws Exception {
    double[] x = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    double[] expected = { -4d, -3d, -2d, -1d, 0d, 1d, 2d, 3d, 4d };
    TimeSeriesUtils.demeanInPlace(x);
    for (int i = 0; i < x.length; i++) {
      assertEquals(x[i], expected[i], 1E-15);
    }
  }

  @Test
  public final void testDetrendLinear2() throws Exception {
    double[] x = { -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6,
        7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20 };

    x = TimeSeriesUtils.detrend(x);
    for (int i = 0; i < x.length; i++) {
      assertEquals(new Double(Math.round(x[i])), new Double(0));
    }
  }

  @Test
  public void testNormByMax() {
    String fname = folder + "kiev-step/_BC0.512.seed";
    try {
      DataBlock db = TimeSeriesUtils.getFirstTimeSeries(fname);
      String startString = "2018-038T15:25:00.0";
      String endString = "2018-038T16:00:00.0";
      long st = TestUtils.timeStringToEpochMilli(startString);
      long ed = TestUtils.timeStringToEpochMilli(endString);
      db.trim(st, ed);
      double[] data = db.getData();
      System.out.println("AVG OF RAW DATA: " + TimeSeriesUtils.getMean(data));
      data = TimeSeriesUtils.demean(data);
      System.out.println("AVG OF DEMEANED DATA: " + TimeSeriesUtils.getMean(data));
      double[] normed = TimeSeriesUtils.normalize(data);
      int minIdx = 0; int maxIdx = 0;
      double min = 0; double max = 0;
      for (int i = 0; i < data.length; ++i) {
        if (data[i] < min) {
          min = data[i];
          minIdx = i;
        } else if (data[i] > max) {
          max = data[i];
          maxIdx = i;
        }
      }
      System.out.println("(max, min): " + max + ", " + min);
      //assertTrue(Math.abs(min) < Math.abs(max));
      assertEquals(1.0, normed[maxIdx], 1E-3);
      assertTrue(Math.abs(normed[minIdx]) < 1.0);
      assertEquals(data[minIdx]/data[maxIdx], normed[minIdx], 1E-3);
      assertEquals(-0.79, normed[minIdx], 0.01);
    } catch (FileNotFoundException | SeedFormatException | CodecException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
      fail();
    }
  }

  @Test
  public void timeDataCorrect() {
    String fname = folder + "random_cal_lowfrq/BHZ.512.seed";
    try {
      String data = TimeSeriesUtils.getMplexNameList(fname).get(0);
      DataBlock db = TimeSeriesUtils.getTimeSeries(fname, data);

      Map<Long, double[]> dataMap = db.getDataMap();
      List<Long> regions = new ArrayList<Long>( dataMap.keySet() );
      Collections.sort(regions);
      for (int i = 0; i < regions.size(); ++i) {
        long time = regions.get(i);
        System.out.println(dataMap.get(time).length);
      }

      Instant timeInst = Instant.ofEpochMilli( db.getStartTime() );
      OffsetDateTime dt = OffsetDateTime.ofInstant(timeInst, ZoneOffset.UTC);
      DateTimeFormatter dtf = DateTimeFormatter.ofPattern("uuuu.MM.dd'T'HH:mm:ss.SSS");
      String inputDate = dtf.format(dt);
      // System.out.println(inputDate);
      String correctDate = "2017.08.02T00:00:00.019";
      assertEquals(inputDate, correctDate);
    } catch (FileNotFoundException | SeedFormatException | CodecException e) {
      e.printStackTrace();
      fail();
    }
  }

}
