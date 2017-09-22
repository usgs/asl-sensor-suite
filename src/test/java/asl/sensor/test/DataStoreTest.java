package asl.sensor.test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

import asl.sensor.gui.InputPanel;
import asl.sensor.input.DataBlock;
import asl.sensor.input.DataStore;
import asl.sensor.utils.TimeSeriesUtils;

public class DataStoreTest {

  public String network = "XX";
  public String station = "TST5";
  public String location = "00";
  public String channel = "BH0";
  
  public String fileID = station+"_"+location+"_"+channel;
  public String filter = network + "_" + fileID;
  
  public String filename1 = "./data/"+fileID+".512.seed";
  
  @Test
  public void commonTimeTrimMatchesLength() {
    
    DataStore ds = new DataStore();
    ds.setBlock(0, filename1, filter);
    
    int left = 250;
    int right = 750;
    
    DataBlock db = ds.getBlock(0);
    
    int oldSize = db.size();
    
    // tested in DataPanelTest
    long loc1 = InputPanel.getMarkerLocation(db, left);
    long loc2 = InputPanel.getMarkerLocation(db, right);
    //  tested in DataBlockTest
    db.trim(loc1, loc2);
    
    ds.setBlock(1, filename1, filter);
    ds.setBlock(2, filename1, filter);
    
    // function under test
    ds.trimToCommonTime();
    
    assertEquals( ds.getBlock(1).getStartTime(), loc1);
    assertEquals( ds.getBlock(1).getEndTime(), loc2);
    assertEquals( db.size(), ds.getBlock(1).size() );
    assertNotEquals( db.size(), oldSize );
    
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