package asl.sensor.test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;
import java.io.FileNotFoundException;
import java.io.IOException;
import org.junit.Before;
import org.junit.Test;
import asl.sensor.gui.InputPanel;
import asl.sensor.input.DataBlock;
import asl.sensor.input.DataStore;
import asl.sensor.utils.TimeSeriesUtils;
import edu.iris.dmc.seedcodec.CodecException;
import edu.sc.seis.seisFile.mseed.SeedFormatException;

public class DataStoreTest {

  public static String folder = TestUtils.DL_DEST_LOCATION + TestUtils.SUBPAGE;

  public String station = "TST5";
  public String location = "00";
  public String channel = "BH0";
  public String fileID = station+"_"+location+"_"+channel+".512.seed";

  @Before
  public void getReferencedData() {

    // place in sprockets folder under 'from-sensor-test/[test-name]'
    String refSubfolder = TestUtils.SUBPAGE + "blocktrim/";
    try {
      TestUtils.downloadTestData(refSubfolder, fileID, refSubfolder, fileID);
    } catch (IOException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }

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
    } catch (FileNotFoundException | SeedFormatException | CodecException e) {
      e.printStackTrace();
      fail();
    }


  }

  @Test
  public void decimationMatchesFrequency() {
    long interval40Hz = (TimeSeriesUtils.ONE_HZ_INTERVAL / 40);
    long interval25Hz = (TimeSeriesUtils.ONE_HZ_INTERVAL / 25);
    // long interval = TimeSeriesUtils.ONE_HZ_INTERVAL;

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