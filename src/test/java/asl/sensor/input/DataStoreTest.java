package asl.sensor.input;

import static asl.sensor.test.TestUtils.RESP_LOCATION;
import static asl.sensor.test.TestUtils.getSeedFolder;
import static asl.utils.NumericUtils.euclidLCM;
import static asl.utils.TimeSeriesUtils.ONE_HZ_INTERVAL;
import static asl.utils.TimeSeriesUtils.getFirstTimeSeries;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertTrue;

import asl.sensor.gui.InputPanel;
import asl.sensor.test.TestUtils;
import asl.utils.input.DataBlock;
import java.time.OffsetDateTime;
import org.junit.Test;

public class DataStoreTest {

  private static final String folder = TestUtils.TEST_DATA_LOCATION + TestUtils.SUBPAGE;

  @Test
  public void trim_BCIPData_timeAndLengthMatch_firstValuesMatch() {
    String respName = RESP_LOCATION + "RESP.CU.BCIP.00.BHZ_2017_268";
    String dataFolderName = getSeedFolder("CU", "BCIP", "2017", "268");
    String calName = dataFolderName + "CB_BC0.512.seed";
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

    double[] dataIn = ds.getBlock(0).getData();
    double[] dataOut = ds.getBlock(1).getData();

    assertEquals(180000, dataIn.length);
    assertEquals(180000, dataOut.length);
    double[] dataOutFrontExpected = {3144, 3193, 3444, 5648, 13154, 23191, 30286, 34405, 34299,
        28745};
    double[] dataInFrontExpected = {19794, 28962, 259083, 561389, 581833, 488582, 333575, -57416,
        -276434, -142344};

    double[] dataInBackExpected = {1074917, 1407337, 1224311, 807605, 385362, -27351,
        -128151, 18800, -84357, -145258};
    double[] dataOutBackExpected = {256572, 235521, 217799, 205593, 202182, 204611, 204693,
        204504, 208574, 210550};

    for (int i = 0; i < 10; i++) {
      assertEquals(dataOutFrontExpected[i], dataOut[i], 1E-5);
      assertEquals(dataInFrontExpected[i], dataIn[i], 1E-5);
      assertEquals(dataOutBackExpected[i], dataOut[180000 - 1 - i], 1E-5);
      assertEquals(dataInBackExpected[i], dataIn[180000 - 1 - i], 1E-5);
    }
  }

  @Test
  public void trimToCommonTime_commonTimeTrimMatchesLength() throws Exception {
    String channel = "BH0";
    String location = "00";
    String station = "TST5";
    String fileID = station + "_" + location + "_" + channel + ".512.seed";
    String filename = folder + "blocktrim/" + fileID;
    DataStore ds = new DataStore();
    DataBlock db;
    db = getFirstTimeSeries(filename);
    String filter = db.getName();
    ds.setBlock(0, db);

    int left = 250;
    int right = 750;

    int oldSize = db.size();

    // tested in DataPanelTest
    long loc1 = InputPanel.getMarkerLocation(db, left);
    long loc2 = InputPanel.getMarkerLocation(db, right);
    db.trim(loc1, loc2);

    ds.setBlock(1, filename, filter);
    ds.setBlock(2, filename, filter);

    // function under test
    ds.trimToCommonTime();

    assertEquals(ds.getBlock(1).getStartTime(), loc1);
    assertEquals(ds.getBlock(1).getEndTime(), loc2);
    assertEquals(db.size(), ds.getBlock(1).size());
    assertNotEquals(db.size(), oldSize);
  }

  @Test
  public void decimationMatchesFrequency() {
    long interval40Hz = ONE_HZ_INTERVAL / 40;
    long interval25Hz = ONE_HZ_INTERVAL / 25;

    long target = euclidLCM(interval25Hz, interval40Hz);
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

    assertEquals(target, ds.getBlock(1).getInterval());
    assertEquals(ds.getBlock(0).size(), ds.getBlock(1).size());
    // make sure that the data has been initialized (i.e., not all 0)
    // if data wasn't being set correctly, result would be all zeros
    boolean notAllZero = false;
    for (Number val : ds.getBlock(1).getData()) {
      if (val.doubleValue() != 0.) {
        notAllZero = true;
      }
    }
    assertTrue(notAllZero);
  }

  @Test
  public void bothComponentsSet_bothSet_True() {
    String respName = RESP_LOCATION + "RESP.CU.BCIP.00.BHZ_2017_268";
    String dataFolderName = getSeedFolder("CU", "BCIP", "2017", "268");
    String sensOutName = dataFolderName + "00_EHZ.512.seed";
    DataStore dataStore = DataStoreUtils.createFromNames(respName, null, sensOutName);

    assertTrue(dataStore.bothComponentsSet(1));
  }

  @Test
  public void bothComponentsSet_dataSet_respMissing_False() {

    String dataFolderName = getSeedFolder("CU", "BCIP", "2017", "268");
    String calName = dataFolderName + "CB_BC0.512.seed";
    DataStore dataStore = DataStoreUtils.createFromNames(null, calName, null);

    assertFalse(dataStore.bothComponentsSet(0));
  }

  @Test
  public void bothComponentsSet_dataMissing_respSet_False() {
    String respName = RESP_LOCATION + "RESP.CU.BCIP.00.BHZ_2017_268";
    DataStore dataStore = DataStoreUtils.createFromNames(respName, null, null);
    assertFalse(dataStore.bothComponentsSet(1));
  }

  @Test
  public void bothComponentsSet_bothMissing_False() {
    String respName = RESP_LOCATION + "RESP.CU.BCIP.00.BHZ_2017_268";
    DataStore dataStore = DataStoreUtils.createFromNames(respName, null, null);
    assertFalse(dataStore.bothComponentsSet(0));
  }

  @Test
  public void removeData_removesSpecificIndex() throws Exception {
    String channel = "BH0";
    String location = "00";
    String station = "TST5";
    String fileID = station + "_" + location + "_" + channel + ".512.seed";
    String filename = folder + "blocktrim/" + fileID;
    DataStore ds = new DataStore();
    DataBlock db;
    db = getFirstTimeSeries(filename);
    ds.setBlock(0, db);
    String respName = RESP_LOCATION + "RESP.CU.BCIP.00.BHZ_2017_268";
    ds.setResponse(0, respName);
    assertTrue(ds.bothComponentsSet(0));
    ds.removeData(0);
    assertFalse(ds.responseIsSet(0));
    assertFalse(ds.blockIsSet(0));
    assertFalse(ds.bothComponentsSet(0));
  }

  @Test
  public void removeBlock_removesSpecificIndex() throws Exception {
    String channel = "BH0";
    String location = "00";
    String station = "TST5";
    String fileID = station + "_" + location + "_" + channel + ".512.seed";
    String filename = folder + "blocktrim/" + fileID;
    DataStore ds = new DataStore();
    DataBlock db;
    db = getFirstTimeSeries(filename);
    ds.setBlock(0, db);
    String respName = RESP_LOCATION + "RESP.CU.BCIP.00.BHZ_2017_268";
    ds.setResponse(0, respName);
    assertTrue(ds.bothComponentsSet(0));
    ds.removeBlock(0);
    assertTrue(ds.responseIsSet(0));
    assertFalse(ds.blockIsSet(0));
    assertFalse(ds.bothComponentsSet(0));
  }

  @Test
  public void isAnythingSet_RESPOnly() {
    String respName = RESP_LOCATION + "RESP.CU.BCIP.00.BHZ_2017_268";
    DataStore dataStore = DataStoreUtils.createFromNames(respName, null, null);
    assertTrue(dataStore.isAnythingSet());
  }

  @Test
  public void isAnythingSet_RESPAndData() {
    String respName = RESP_LOCATION + "RESP.CU.BCIP.00.BHZ_2017_268";
    String dataFolderName = getSeedFolder("CU", "BCIP", "2017", "268");
    String sensOutName = dataFolderName + "00_EHZ.512.seed";
    DataStore dataStore = DataStoreUtils.createFromNames(respName, null, sensOutName);
    assertTrue(dataStore.isAnythingSet());
  }

  @Test
  public void isAnythingSet_dataOnly() {
    String dataFolderName = getSeedFolder("CU", "BCIP", "2017", "268");
    String sensOutName = dataFolderName + "00_EHZ.512.seed";
    DataStore dataStore = DataStoreUtils.createFromNames(null, null, sensOutName);
    assertTrue(dataStore.isAnythingSet());
  }

  @Test
  public void isAnythingSet_nothingSet() {
    DataStore dataStore = DataStoreUtils.createFromNames(null, null, null);
    assertFalse(dataStore.isAnythingSet());
  }

}