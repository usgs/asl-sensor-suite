package asl.sensor.gui;

import static org.junit.Assert.assertEquals;

import asl.sensor.test.TestUtils;
import asl.utils.timeseries.TimeSeriesUtils;
import asl.utils.timeseries.DataBlock;
import org.junit.Test;

public class DataPanelTest {

  private static final String folder = TestUtils.TEST_DATA_LOCATION + TestUtils.SUBPAGE;

  @Test
  public void getsCorrectTrimming() throws Exception {
    int left = InputPanel.SLIDER_MAX / 4;
    int right = 3 * InputPanel.SLIDER_MAX / 4;
    int farLeft = 0;

    String channel = "BH0";
    String location = "00";
    String station = "TST5";
    String fileID = station + "_" + location + "_" + channel + ".512.seed";
    String filename = folder + "blocktrim/" + fileID;

    DataBlock db = TimeSeriesUtils.getFirstTimeSeries(filename);
    long start = db.getStartTime();
    long end = db.getEndTime();
    long interval = db.getInterval();
    int size = db.size();

    long timeRange = interval * size;
    long timeRangeFromExtremes = end - start;
    assertEquals(timeRange, timeRangeFromExtremes);

    long loc1 = InputPanel.getMarkerLocation(db, left);
    long loc2 = InputPanel.getMarkerLocation(db, right);
    long loc3 = InputPanel.getMarkerLocation(db, farLeft);

    assertEquals(loc3, start);
    assertEquals(loc2 - loc1, timeRange / 2); // range is 3/4-1/4 = 1/2 of data
    assertEquals(loc1, start + (interval * size / 4)); // correct start pt?
  }

}
