package asl.sensor.test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;
import java.io.FileNotFoundException;
import org.junit.Test;
import asl.sensor.gui.InputPanel;
import asl.sensor.input.DataBlock;
import asl.sensor.utils.TimeSeriesUtils;
import edu.iris.dmc.seedcodec.CodecException;
import edu.sc.seis.seisFile.mseed.SeedFormatException;

public class DataPanelTest {

  public String station = "TST5";
  public String location = "00";
  public String channel = "BH0";

  public String fileID = station+"_"+location+"_"+channel;

  public String filename1 = "./test-data/blocktrim/"+fileID+".512.seed";

  @Test
  public void getsCorrectTrimming() {
    int left = InputPanel.SLIDER_MAX / 4;
    int right = 3 * InputPanel.SLIDER_MAX / 4;
    int farLeft = 0;

    System.out.println(left+","+right);

    try {
      DataBlock db = TimeSeriesUtils.getFirstTimeSeries(filename1);
      long start = db.getStartTime();
      long end = db.getEndTime();
      long interval = db.getInterval();
      // long end = db.getEndTime();
      int size = db.size();

      long timeRange = interval*size;
      long timeRangeFromExtremes = end - start;
      assertEquals(timeRange, timeRangeFromExtremes);

      long loc1 = InputPanel.getMarkerLocation(db, left);
      long loc2 = InputPanel.getMarkerLocation(db, right);
      long loc3 = InputPanel.getMarkerLocation(db, farLeft);

      assertEquals(loc3, start);
      assertEquals(loc2 - loc1, timeRange/2); // range is 3/4-1/4 = 1/2 of data
      assertEquals(loc1, start + (interval * size / 4) ); // correct start pt?
    } catch (FileNotFoundException | SeedFormatException | CodecException e) {
      e.printStackTrace();
      fail();
    }



  }

}
