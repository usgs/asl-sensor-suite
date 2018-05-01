package asl.sensor.input;

import static asl.sensor.test.TestUtils.RESP_LOCATION;
import static asl.sensor.test.TestUtils.getSeedFolder;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import asl.sensor.test.TestUtils;
import java.io.IOException;
import java.time.OffsetDateTime;
import java.util.ArrayList;
import java.util.List;
import org.junit.Test;
import asl.sensor.gui.InputPanel;
import asl.sensor.input.DataBlock;
import asl.sensor.input.DataStore;
import asl.sensor.utils.TimeSeriesUtils;
import edu.iris.dmc.seedcodec.CodecException;
import edu.sc.seis.seisFile.mseed.SeedFormatException;

public class DataStoreTest {

  public static String folder = TestUtils.TEST_DATA_LOCATION + TestUtils.SUBPAGE;

  public String station = "TST5";
  public String location = "00";
  public String channel = "BH0";
  public String fileID = station+"_"+location+"_"+channel+".512.seed";

  @Test
  public void trim_BCIPData_timeAndLengthMatch() {
    String respName = RESP_LOCATION + "RESP.CU.BCIP.00.BHZ_2017_268";
    String dataFolderName = getSeedFolder("CU", "BCIP", "2017", "268");
    String calName =  dataFolderName + "CB_BC0.512.seed";
    String sensOutName = dataFolderName + "00_EHZ.512.seed";

    DataStore ds = DataStoreUtils.createFromNames(respName, calName, sensOutName);
    OffsetDateTime cCal = TestUtils.getStartCalendar(ds);
    cCal = cCal.withHour(18).withMinute(49).withSecond(0).withNano(0);
    long start = cCal.toInstant().toEpochMilli();

    cCal = cCal.withHour(19).withMinute(4);
    long end = cCal.toInstant().toEpochMilli();

    ds.trim(start, end);

    assertEquals(start, ds.getBlock(0).getStartTime());
    assertEquals(end, ds.getBlock(0).getEndTime());

    assertEquals(start, ds.getBlock(1).getStartTime());
    assertEquals(end, ds.getBlock(1).getEndTime());

    assertEquals(180000, ds.getBlock(0).getData().length);
    assertEquals(180000, ds.getBlock(1).getData().length);
  }

  @Test
  public void commonTimeTrimMatchesLength() {
    String filename = folder + "blocktrim/" + fileID;
    DataStore ds = new DataStore();
    DataBlock db;
    try {
      db = TimeSeriesUtils.getFirstTimeSeries(filename);
      String filter = db.getName();
      ds.setBlock(0, db);

      int left = 250;
      int right = 750;

      int oldSize = db.size();

      // tested in DataPanelTest
      long loc1 = InputPanel.getMarkerLocation(db, left);
      long loc2 = InputPanel.getMarkerLocation(db, right);
      //  tested in DataBlockTest
      db.trim(loc1, loc2);

      ds.setBlock(1, filename, filter);
      ds.setBlock(2, filename, filter);

      // function under test
      ds.trimToCommonTime();

      assertEquals( ds.getBlock(1).getStartTime(), loc1);
      assertEquals( ds.getBlock(1).getEndTime(), loc2);
      assertEquals( db.size(), ds.getBlock(1).size() );
      assertNotEquals( db.size(), oldSize );
    } catch (IOException | SeedFormatException | CodecException e) {
      e.printStackTrace();
      fail();
    }
  }

  @Test
  public void decimationMatchesFrequency() {
    long interval40Hz = (TimeSeriesUtils.ONE_HZ_INTERVAL / 40);
    long interval25Hz = (TimeSeriesUtils.ONE_HZ_INTERVAL / 25);

    long start = 0;

    // range of 4 seconds
    double[] series25Hz = new double[100];
    double[] series40Hz = new double[160];

    for (int i = 0; i < 100; ++i) {
      series25Hz[i] = i * Math.sin(i);
    }

    for (int i = 0; i < 160; ++i) {
      series40Hz[i] = i * Math.sin(i);
    }

    DataBlock block25Hz = new DataBlock(series25Hz, interval25Hz, "25", start);
    DataBlock block40Hz = new DataBlock(series40Hz, interval40Hz, "40", start);

    DataStore ds = new DataStore();
    ds.setBlock(0, block25Hz);
    ds.setBlock(1, block40Hz);
    ds.matchIntervals();

    assertEquals( ds.getBlock(1).getInterval(), interval25Hz );
    assertEquals( ds.getBlock(0).size(), ds.getBlock(1).size() );
    // make sure that the data has been initialized (i.e., not all 0)
    // if data wasn't being set correctly, result would be all zeros
    boolean notAllZero = false;
    for ( Number val : ds.getBlock(1).getData() ) {
      if (val.doubleValue() != 0.) {
        notAllZero = true;
      }
    }
    assertTrue(notAllZero);
  }

}