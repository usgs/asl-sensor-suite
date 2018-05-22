package asl.sensor.input;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertTrue;

import asl.sensor.gui.InputPanel;
import asl.sensor.test.TestUtils;
import asl.sensor.utils.TimeSeriesUtils;
import org.junit.Test;

public class DataBlockTest {

  private static final String folder = TestUtils.TEST_DATA_LOCATION + TestUtils.SUBPAGE;

  private final String station = "TST5";
  private final String location = "00";
  private final String channel = "BH0";
  private final String fileID = station + "_" + location + "_" + channel + ".512.seed";

  @Test
  public void trimsCorrectly() throws Exception {
    int left = InputPanel.SLIDER_MAX / 4;
    int right = 3 * InputPanel.SLIDER_MAX / 4;

    String filename = folder + "blocktrim/" + fileID;

    DataBlock db = TimeSeriesUtils.getFirstTimeSeries(filename);

    int sizeOld = db.size();

    // these get tested in DataPanelTest
    long loc1 = InputPanel.getMarkerLocation(db, left);
    long loc2 = InputPanel.getMarkerLocation(db, right);

    db.trim(loc1, loc2);

    assertEquals(loc1, db.getStartTime());
    assertEquals(sizeOld / 2, db.size());
  }

  @Test
  public void trimsCOWICorrectly() throws Exception {
    String filename = folder + "cowi-multitests/C100823215422_COWI.LHx";
    String dataname = "US_COWI_  _LHN";
    DataBlock db = TimeSeriesUtils.getTimeSeries(filename, dataname);
    String startString = "2010-236T02:00:00.0";
    String endString = "2010-236T13:00:00.0";
    long st = TestUtils.timeStringToEpochMilli(startString);
    long ed = TestUtils.timeStringToEpochMilli(endString);
    db.trim(st, ed);
    double[] data = db.getData();
    assertEquals(39600, data.length);
  }

  @Test
  public void appendsCorrectly() throws Exception {
    String subfolder = folder + "test-appending/";
    String fullLength = "_BC0.512.seed";
    String firstPart = "044._BC0.512.seed";
    String secondPart = "045._BC0.512.seed";
    DataBlock dbFull, dbAppended;
    dbFull = TimeSeriesUtils.getFirstTimeSeries(subfolder + fullLength);
    dbAppended = TimeSeriesUtils.getFirstTimeSeries(subfolder + firstPart);
    dbAppended.appendTimeSeries(subfolder + secondPart);
    long st = dbFull.getStartTime();
    long ed = dbFull.getEndTime();
    assertEquals(st, dbAppended.getStartTime());
    assertEquals(ed, dbAppended.getEndTime());
  }

  @Test
  public void unTrim_bothStartEndTimesReset() throws Exception {
    String filename = folder + "blocktrim/" + fileID;
    DataBlock block = TimeSeriesUtils.getFirstTimeSeries(filename);
    double[] initialData = block.getData();
    int initialLength = initialData.length;

    long initialStart = block.getStartTime();
    long initialEnd = block.getEndTime();
    // trim a quarter of the length off either end
    long quarterTrimDifference = (initialEnd - initialStart) / 4;
    long trimmedStart = initialStart + quarterTrimDifference;
    long trimmedEnd = initialEnd - quarterTrimDifference;
    block.trim(trimmedStart, trimmedEnd);
    int trimmedLength = block.getData().length;

    block.untrim();
    long untrimmedStart = block.getStartTime();
    long untrimmedEnd = block.getEndTime();
    //Was rebuildList set correctly too?
    double[] untrimmedData = block.getData();
    int untrimmedLength = untrimmedData.length;

    assertEquals(initialStart, untrimmedStart);
    assertEquals(initialEnd, untrimmedEnd);
    assertEquals(initialLength, untrimmedLength);
    assertArrayEquals(initialData, untrimmedData, 1E-6);
    assertFalse(initialData == untrimmedData);
    assertTrue(trimmedLength < untrimmedLength);
  }

  @Test
  public void unTrim_startTimeReset() throws Exception {
    String filename = folder + "blocktrim/" + fileID;
    DataBlock block = TimeSeriesUtils.getFirstTimeSeries(filename);
    double[] initialData = block.getData();
    int initialLength = initialData.length;

    long initialStart = block.getStartTime();
    long initialEnd = block.getEndTime();
    // trim a quarter of the length off
    long quarterTrimDifference = (initialEnd - initialStart) / 4;
    long trimmedStart = initialStart + quarterTrimDifference;
    block.trim(trimmedStart, initialEnd);
    int trimmedLength = block.getData().length;

    block.untrim();
    long untrimmedStart = block.getStartTime();
    //Was rebuildList set correctly too?
    double[] untrimmedData = block.getData();
    int untrimmedLength = untrimmedData.length;

    assertEquals(initialStart, untrimmedStart);
    assertEquals(initialLength, untrimmedLength);
    assertArrayEquals(initialData, untrimmedData, 1E-6);
    assertFalse(initialData == untrimmedData);
    assertTrue(trimmedLength < untrimmedLength);
  }

  @Test
  public void unTrim_endTimeReset() throws Exception {
    String filename = folder + "blocktrim/" + fileID;
    DataBlock block = TimeSeriesUtils.getFirstTimeSeries(filename);
    double[] initialData = block.getData();
    int initialLength = initialData.length;

    long initialStart = block.getStartTime();
    long initialEnd = block.getEndTime();
    // trim a quarter of the length off
    long quarterTrimDifference = (initialEnd - initialStart) / 4;
    long trimmedEnd = initialEnd - quarterTrimDifference;
    block.trim(initialStart, trimmedEnd);
    int trimmedLength = block.getData().length;

    block.untrim();
    long untrimmedEnd = block.getEndTime();
    //Was rebuildList set correctly too?
    double[] untrimmedData = block.getData();
    int untrimmedLength = untrimmedData.length;

    assertEquals(initialEnd, untrimmedEnd);
    assertEquals(initialLength, untrimmedLength);
    assertArrayEquals(initialData, untrimmedData, 1E-6);
    assertTrue(trimmedLength < untrimmedLength);
  }

  @Test
  public void unTrim_wasNotTrimmed_rebuildWasTrue() throws Exception {
    String filename = folder + "blocktrim/" + fileID;
    DataBlock block = TimeSeriesUtils.getFirstTimeSeries(filename);
    double[] initialData = block.getData();
    int initialLength = initialData.length;

    long initialStart = block.getStartTime();
    long initialEnd = block.getEndTime();

    block.trim(initialStart, initialEnd);

    block.untrim();
    long untrimmedStart = block.getStartTime();
    long untrimmedEnd = block.getEndTime();
    //Was rebuildList set correctly too?
    double[] untrimmedData = block.getData();
    int untrimmedLength = untrimmedData.length;

    assertEquals(initialStart, untrimmedStart);
    assertEquals(initialEnd, untrimmedEnd);
    assertEquals(initialLength, untrimmedLength);
    assertTrue(initialData == untrimmedData);
    assertEquals(initialData, untrimmedData);
  }

  @Test
  public void resample_upsamplingIgnored() throws Exception {
    String filename = folder + "blocktrim/" + fileID;
    DataBlock block = TimeSeriesUtils.getFirstTimeSeries(filename);
    long initialInterval = block.getInitialInterval();
    double[] initialData = block.getData();
    block.resample(block.getInitialInterval() / 2);
    assertEquals(initialInterval, block.getInterval());
    assertArrayEquals(initialData, block.getData(), 1E-10);
    assertTrue(initialData == block.getData());
  }

  @Test
  public void resample_decimation() throws Exception {
    String filename = folder + "blocktrim/" + fileID;
    DataBlock block = TimeSeriesUtils.getFirstTimeSeries(filename);
    long initialInterval = block.getInitialInterval();
    double[] initialData = block.getData();
    long slowerInterval = initialInterval * 2;
    block.resample(slowerInterval);
    assertNotEquals(initialInterval, block.getInterval());
    assertEquals(initialData.length, block.getData().length * 2);
  }

  @Test
  public void getSampleRate_basicTest() throws Exception {
    String filename = folder + "blocktrim/" + fileID;
    DataBlock block = TimeSeriesUtils.getFirstTimeSeries(filename);
    assertEquals(40, block.getSampleRate(), 1E-10);
  }

  @Test
  public void getData_rebuildTrue_doesItRebuildCorrectly() throws Exception {
    String filename = folder + "blocktrim/" + fileID;
    DataBlock block = TimeSeriesUtils.getFirstTimeSeries(filename);
    double[] initialData = block.getData();
    int initialLength = initialData.length;

    long initialStart = block.getStartTime();
    long initialEnd = block.getEndTime();
    // trim a quarter of the length off either end
    long quarterTrimDifference = (initialEnd - initialStart) / 4;
    long trimmedStart = initialStart + quarterTrimDifference;
    long trimmedEnd = initialEnd - quarterTrimDifference;
    block.trim(trimmedStart, trimmedEnd);
    int trimmedLength = block.getData().length;

    block.untrim();
    //Was rebuildList set correctly too?
    double[] untrimmedData = block.getData();
    int untrimmedLength = untrimmedData.length;

    assertEquals(initialLength, untrimmedLength);
    assertArrayEquals(initialData, untrimmedData, 1E-6);
    assertFalse(initialData == untrimmedData);
    assertTrue(trimmedLength < untrimmedLength);
  }

  @Test
  public void dataBlock_dataBlockIn_didItCloneCorrectly() throws Exception {
    String filename = folder + "blocktrim/" + fileID;
    DataBlock block = TimeSeriesUtils.getFirstTimeSeries(filename);
    DataBlock clonedBlock = new DataBlock(block);
    assertEquals(block.getName(), clonedBlock.getName());
    assertFalse(block.getDataMap() == clonedBlock.getDataMap());
    assertFalse(block.getData() == clonedBlock.getData());
    assertArrayEquals(block.getData(), clonedBlock.getData(), 1E-10);
    assertEquals(block.getStartTime(), clonedBlock.getStartTime());
    assertEquals(block.getEndTime(), clonedBlock.getEndTime());
    assertEquals(block.getInterval(), clonedBlock.getInterval());
    // another assert to check that objects are in fact distinct -- different interval values
    clonedBlock.resample(clonedBlock.getInterval() * 2);
    assertNotEquals(block.getInterval(), clonedBlock.getInterval());
  }

}