package asl.sensor.utils;

import static asl.sensor.utils.TimeSeriesUtils.concatAll;
import static asl.sensor.utils.TimeSeriesUtils.downsample;
import static asl.sensor.utils.TimeSeriesUtils.euclidGCD;
import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import asl.sensor.input.DataBlock;
import asl.sensor.test.TestUtils;
import edu.iris.dmc.seedcodec.B1000Types;
import edu.iris.dmc.seedcodec.CodecException;
import edu.iris.dmc.seedcodec.DecompressedData;
import edu.sc.seis.seisFile.mseed.DataRecord;
import edu.sc.seis.seisFile.mseed.SeedFormatException;
import edu.sc.seis.seisFile.mseed.SeedRecord;
import java.io.BufferedInputStream;
import java.io.DataInput;
import java.io.DataInputStream;
import java.io.EOFException;
import java.io.FileInputStream;
import java.io.IOException;
import java.time.Instant;
import java.time.OffsetDateTime;
import java.time.ZoneOffset;
import java.time.format.DateTimeFormatter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Set;
import org.junit.Test;

public class TimeSeriesUtilsTest {

  private static final String folder = TestUtils.TEST_DATA_LOCATION + TestUtils.SUBPAGE;

  private final String station = "TST5";
  private final String location = "00";
  private final String channel = "BH0";
  private final String fileID = station + "_" + location + "_" + channel + ".512.seed";

  @Test
  public void canGetFile() {
    try {
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
      assertEquals(names.size(), 3);
    } catch (IOException | SeedFormatException e) {
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

    double[] numbers = {1, 2, 3, 4, 5};

    double[] numList = numbers.clone();
    double[] demeaned = numList.clone();

    TimeSeriesUtils.demeanInPlace(demeaned);

    for (int i = 0; i < numList.length; ++i) {
      assertEquals(demeaned[i], numList[i] - 3, 1E-15);
    }

  }

  @Test
  public void detrendAtEndsTest() {
    double[] x = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,
        18, 19, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4,
        3, 2, 1};

    x = TimeSeriesUtils.detrendEnds(x);

    assertEquals(x[0], 0, 1E-5);
    assertEquals(x[x.length - 1], 0, 1E-5);
  }

  @Test
  public final void testDetrendLinear() {
    double[] x = {-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,
        18, 19, 20};

    x = TimeSeriesUtils.detrend(x);
    for (double aX : x) {
      assertEquals(aX, 0.0, 1E-5);
    }
  }

  @Test
  public void detrendingCycleTest() {
    /*
     * i = 0 to 15
     * x = sin(i)+2
     * Best Fit: y=0.0044008753x + 2.0879739023
     */
    double[] x = {2.0, 2.8414709848, 2.9092974268, 2.1411200081, 1.2431975047, 1.0410757253,
        1.7205845018, 2.6569865987, 2.9893582466, 2.4121184852, 1.4559788891, 1.0000097934,
        1.463427082, 2.4201670368, 2.9906073557, 2.6502878402};
    double[] answer = {-0.0879739023, 0.7490962071, 0.8125217738, 0.0399434797, -0.862379899,
        -1.0689025537, -0.3937946526, 0.5382065689, 0.8661773415, 0.2845367048,
        -0.6760037667, -1.1363737377, -0.6773573245, 0.2749817549, 0.8410211985,
        0.4963008076};

    x = TimeSeriesUtils.detrend(x);

    assertArrayEquals(x, answer, 1E-5);
  }

  @Test
  public final void testDetrendNoSlope() {
    double[] x = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    x = TimeSeriesUtils.detrend(x);
    for (double aX : x) {
      assertEquals(aX, 0.0, 1E-5);
    }
  }

  @Test
  public final void downsample_thirdOfFrequency() {
    double[] in = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
    int factor = 3;
    double[] expected = {1, 4, 7, 10};
    double[] out = downsample(in, factor);
    assertArrayEquals(expected, out, 1E-5);
  }

  @Test
  public final void euclidGCD_basicTest() {
    assertEquals(50L, euclidGCD(100, 50));
    assertEquals(50L, euclidGCD(50, 100));
    assertEquals(3L, euclidGCD(999, 105));
    assertEquals(2L, euclidGCD(2147483612, 2137483646));
  }

  @Test
  public void concatAllTest() {
    double[] array1 = {1, 2, 3, 4.5, 6, 7};
    double[] array2 = {111, 112, 113, 114.5, 116, 117};
    double[] array3 = {11, 12, 13, 14.5, 16, 17};
    double[] arrayAnswer = {1, 2, 3, 4.5, 6, 7, 111, 112, 113, 114.5, 116, 117, 11, 12, 13, 14.5,
        16, 17};
    List<double[]> list = new ArrayList<>();
    list.add(array1);
    list.add(array2);
    list.add(array3);
    double[] arrayResult = concatAll(list);
    assertArrayEquals(arrayAnswer, arrayResult, 1E-10);
  }

  @Test
  public void concatAllTestEmpty() {
    double[] arrayAnswer = {};
    List<double[]> list = new ArrayList<>();
    double[] arrayResult = concatAll(list);
    assertArrayEquals(arrayAnswer, arrayResult, 1E-10);
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

    } catch (IOException | SeedFormatException | CodecException e) {
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

    } catch (IOException | SeedFormatException | CodecException e) {
      fail();
      e.printStackTrace();
    }
  }

  @Test
  public void inputFileReaderCreatesXYSeries() {
    DataInput dis;
    List<Number> data = new ArrayList<>();
    String filename1 = folder + "blocktrim/" + fileID;
    try {
      dis = new DataInputStream(
          new BufferedInputStream(
              new FileInputStream(filename1)));
      while (true) {

        try {
          // long interval = 0L;
          SeedRecord sr = SeedRecord.read(dis, 4096);
          if (sr instanceof DataRecord) {
            DataRecord dr = (DataRecord) sr;
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
                for (int dataPoint : decomArrayInt) {
                  data.add(dataPoint);
                }
                break;
              case B1000Types.FLOAT:
                float[] decomArrayFlt = decomp.getAsFloat();
                for (float dataPoint : decomArrayFlt) {
                  data.add(dataPoint);
                }
                break;
              case B1000Types.SHORT:
                short[] decomArrayShr = decomp.getAsShort();
                for (short dataPoint : decomArrayShr) {
                  data.add(dataPoint);
                }
                break;
              default:
                double[] decomArrayDbl = decomp.getAsDouble();
                for (double dataPoint : decomArrayDbl) {
                  data.add(dataPoint);
                }
                break;
            }
          }
        } catch (EOFException e) {
          break;
        }

      }

      // quickly get the one name in the list
      Set<String> names = TimeSeriesUtils.getMplexNameSet(filename1);
      List<String> nameList = new ArrayList<>(names);
      System.out.println("DATA BLOCK SIZE: " + data.size());

      DataBlock testAgainst =
          TimeSeriesUtils.getTimeSeries(filename1, nameList.get(0));
      assertEquals(data.size(), testAgainst.getData().length);

    } catch (IOException | SeedFormatException | CodecException e) {
      assertNull(e);
    }
  }

  @Test
  public void seisFileCanParseFile() {
    String filename1 = folder + "blocktrim/" + fileID;
    try {
      DataInput dis = new DataInputStream(new BufferedInputStream(
          new FileInputStream(filename1)));
      try {
        while (true) {
          SeedRecord sr = SeedRecord.read(dis, 4096);
          if (sr instanceof DataRecord) {
            DataRecord dr = (DataRecord) sr;

            String loc = dr.getHeader().getLocationIdentifier();
            assertTrue(loc.equals(location));
            String stat = dr.getHeader().getStationIdentifier().trim();
            assertTrue(stat.equals(station));

            String chan = dr.getHeader().getChannelIdentifier();
            assertTrue(chan.equals(channel));
          }
        }
      } catch (EOFException e) {
        assertNotNull(e); // I haaates it! I haaaaaaaaaates it!
      } catch (SeedFormatException | IOException e) {
        assertNull(e);
      }

    } catch (IOException e) {
      assertNull(e);
    }

  }

  @Test
  public void seisFileGivesCorrectSampleRateAndInterval() {
    DataInput dis;
    String filename1 = folder + "blocktrim/" + fileID;
    try {
      while (true) {
        dis = new DataInputStream(new BufferedInputStream(
            new FileInputStream(filename1)));
        SeedRecord sr = SeedRecord.read(dis, 4096);
        if (sr instanceof DataRecord) {
          DataRecord dr = (DataRecord) sr;

          int fact = dr.getHeader().getSampleRateFactor();
          int mult = dr.getHeader().getSampleRateMultiplier();

          double rate = dr.getSampleRate();
          assertTrue((double) fact / mult == rate);

          // checking the correct values for the intervals
          double multOf1Hz = rate / TimeSeriesUtils.ONE_HZ;
          long inverse = TimeSeriesUtils.ONE_HZ_INTERVAL / (long) multOf1Hz;

          long interval = TimeSeriesUtils.ONE_HZ_INTERVAL * mult / fact;

          assertEquals(inverse, interval);
          break;
        }
      }
    } catch (IOException | SeedFormatException e) {
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
    } catch (IOException | SeedFormatException | CodecException e) {
      e.printStackTrace();
      fail();
    }
  }

  @Test
  public final void testDemean1to9() {
    double[] x = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    double[] expected = {-4d, -3d, -2d, -1d, 0d, 1d, 2d, 3d, 4d};
    TimeSeriesUtils.demeanInPlace(x);
    for (int i = 0; i < x.length; i++) {
      assertEquals(x[i], expected[i], 1E-15);
    }
  }

  @Test
  public final void testDetrendLinear2() {
    double[] x = {-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6,
        7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20};

    x = TimeSeriesUtils.detrend(x);
    for (double aX : x) {
      assertEquals(Math.round(aX), 0, 1E-7);
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
      int minIdx = 0;
      int maxIdx = 0;
      double min = 0;
      double max = 0;
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
      assertEquals(data[minIdx] / data[maxIdx], normed[minIdx], 1E-3);
      assertEquals(-0.79, normed[minIdx], 0.01);
    } catch (IOException | SeedFormatException | CodecException e) {
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

      Instant timeInst = Instant.ofEpochMilli(db.getStartTime());
      OffsetDateTime dt = OffsetDateTime.ofInstant(timeInst, ZoneOffset.UTC);
      DateTimeFormatter dtf = DateTimeFormatter.ofPattern("uuuu.MM.dd'T'HH:mm:ss.SSS");
      String inputDate = dtf.format(dt);
      // System.out.println(inputDate);
      String correctDate = "2017.08.02T00:00:00.019";
      assertEquals(inputDate, correctDate);
    } catch (IOException | SeedFormatException | CodecException e) {
      e.printStackTrace();
      fail();
    }
  }

}
