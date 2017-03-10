package asl.sensor.test;

import static org.junit.Assert.*;

import java.io.FileNotFoundException;
import java.util.ArrayList;

import org.junit.Test;

import asl.sensor.gui.InputPanel;
import asl.sensor.input.DataBlock;
import asl.sensor.utils.TimeSeriesUtils;

public class DataBlockTest {

  public String station = "TST5";
  public String location = "00";
  public String channel = "BH0";
  
  public String fileID = station+"_"+location+"_"+channel;
  
  public String filename1 = "./data/"+fileID+".512.seed";
  
  @Test
  public void trimsCorrectly() {
    int left = InputPanel.SLIDER_MAX / 4;
    int right = 3 * InputPanel.SLIDER_MAX / 4;
    
    String name;
    try {
      name = new ArrayList<String>( 
            TimeSeriesUtils.getMplexNameSet(filename1)
          ).get(0);
      
      DataBlock db = TimeSeriesUtils.getTimeSeries(filename1, name);

      int sizeOld = db.size();
      
      
      // these get tested in DataPanelTest
      long loc1 = InputPanel.getMarkerLocation(db, left);
      long loc2 = InputPanel.getMarkerLocation(db, right);
      
      db.trim(loc1, loc2);
      
      assertEquals( loc1, db.getStartTime() );
      assertEquals( sizeOld/2, db.size() );
      
    } catch (FileNotFoundException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
      fail();
    }
    

  }
  
}