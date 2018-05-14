package asl.sensor.input;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;

import asl.sensor.gui.InputPanel;
import asl.sensor.test.TestUtils;
import asl.sensor.utils.TimeSeriesUtils;
import edu.iris.dmc.seedcodec.CodecException;
import edu.sc.seis.seisFile.mseed.SeedFormatException;
import java.io.FileNotFoundException;
import java.io.IOException;
import org.junit.Test;

public class DataBlockTest {

  public static String folder = TestUtils.TEST_DATA_LOCATION + TestUtils.SUBPAGE;

  public String station = "TST5";
  public String location = "00";
  public String channel = "BH0";
  public String fileID = station + "_" + location + "_" + channel + ".512.seed";

  @Test
  public void trimsCorrectly() {
    int left = InputPanel.SLIDER_MAX / 4;
    int right = 3 * InputPanel.SLIDER_MAX / 4;

    String filename = folder + "blocktrim/" + fileID;

    try {
      DataBlock db = TimeSeriesUtils.getFirstTimeSeries(filename);

      int sizeOld = db.size();

      // these get tested in DataPanelTest
      long loc1 = InputPanel.getMarkerLocation(db, left);
      long loc2 = InputPanel.getMarkerLocation(db, right);

      db.trim(loc1, loc2);

      assertEquals(loc1, db.getStartTime());
      assertEquals(sizeOld / 2, db.size());

    } catch (SeedFormatException | CodecException | IOException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }

  }

  @Test
  public void trimsCOWICorrectly() {
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
      assertEquals(39600, data.length);
    } catch (FileNotFoundException | SeedFormatException | CodecException e) {
      e.printStackTrace();
      fail();
    }
  }

  @Test
  public void appendsCorrectly() {
    String subfolder = folder + "test-appending/";
    String fullLength = "_BC0.512.seed";
    String firstPart = "044._BC0.512.seed";
    String secondPart = "045._BC0.512.seed";
    DataBlock dbFull, dbAppended;
    try {
      dbFull = TimeSeriesUtils.getFirstTimeSeries(subfolder + fullLength);
      dbAppended = TimeSeriesUtils.getFirstTimeSeries(subfolder + firstPart);
      dbAppended.appendTimeSeries(subfolder + secondPart);
      long st = dbFull.getStartTime();
      long ed = dbFull.getEndTime();
      assertEquals(st, dbAppended.getStartTime());
      assertEquals(ed, dbAppended.getEndTime());
    } catch (SeedFormatException | CodecException | IOException e) {
      e.printStackTrace();
      fail();
    }
  }

}