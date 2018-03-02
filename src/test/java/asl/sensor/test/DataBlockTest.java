package asl.sensor.test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.time.LocalDateTime;
import java.time.ZoneOffset;
import java.time.format.DateTimeFormatter;
import org.junit.Before;
import org.junit.Test;
import asl.sensor.gui.InputPanel;
import asl.sensor.input.DataBlock;
import asl.sensor.utils.TimeSeriesUtils;
import edu.iris.dmc.seedcodec.CodecException;
import edu.sc.seis.seisFile.mseed.SeedFormatException;

public class DataBlockTest {

  public static String folder = TestUtils.DL_DEST_LOCATION + TestUtils.SUBPAGE;

  public String station = "TST5";
  public String location = "00";
  public String channel = "BH0";
  public String fileID = station+"_"+location+"_"+channel+".512.seed";

  @Before
  public void getReferencedData() {

    // place in sprockets folder under 'from-sensor-test/[test-name]'

    String refSubfolder = TestUtils.SUBPAGE + "cowi-multitests/";
    String filename = "C100823215422_COWI.LHx";
    String filename2 = "DT000110.LH1";
    try {
      TestUtils.downloadTestData(refSubfolder, filename, refSubfolder, filename);
      TestUtils.downloadTestData(refSubfolder, filename2, refSubfolder, filename2);
    } catch (IOException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }

    refSubfolder = TestUtils.SUBPAGE + "blocktrim/";
    try {
      TestUtils.downloadTestData(refSubfolder, fileID, refSubfolder, fileID);
    } catch (IOException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }

  }

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

      assertEquals( loc1, db.getStartTime() );
      assertEquals( sizeOld/2, db.size() );

    } catch (FileNotFoundException e) {
      e.printStackTrace();
      fail();
    } catch (SeedFormatException | CodecException e) {
      e.printStackTrace();
      fail();
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
        DateTimeFormatter dtf = DateTimeFormatter.ofPattern("uuuu-DDD'T'HH:mm:ss.S");
        long st = LocalDateTime.parse(startString, dtf).toInstant(ZoneOffset.UTC).toEpochMilli();
        long ed = LocalDateTime.parse(endString, dtf).toInstant(ZoneOffset.UTC).toEpochMilli();
        db.trim(st, ed);
        double[] data = db.getData();
        assertEquals(39600, data.length);
      } catch (FileNotFoundException e) {
        e.printStackTrace();
        fail();
      } catch (SeedFormatException | CodecException e) {
        e.printStackTrace();
        fail();
      }

  }

}