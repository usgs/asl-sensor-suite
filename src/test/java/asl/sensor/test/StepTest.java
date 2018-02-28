package asl.sensor.test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;
import java.time.LocalDateTime;
import java.time.ZoneOffset;
import java.time.format.DateTimeFormatter;
import org.junit.Test;
import asl.sensor.experiment.StepExperiment;
import asl.sensor.input.DataStore;
import edu.iris.dmc.seedcodec.CodecException;
import edu.sc.seis.seisFile.mseed.SeedFormatException;

public class StepTest {

  @Test
  public void testKievGoodCorner() {
    DataStore ds = new DataStore();
    String folder = "test-data/kiev-step/";
    String fname1 = "_BC0.512.seed";
    String fname2 = "00_BHZ.512.seed";
    try {
      ds.setBlock(0, folder + fname1);
      ds.setBlock(1, folder + fname2);
      ds.setEmbedResponse(1, "STS1T5_Q330HR");
      String startString = "2018-038T15:25:00.0";
      String endString = "2018-038T16:00:00.0";
      DateTimeFormatter dtf = DateTimeFormatter.ofPattern("uuuu-DDD'T'HH:mm:ss.S");
      long st = LocalDateTime.parse(startString, dtf).toInstant(ZoneOffset.UTC).toEpochMilli();
      long ed = LocalDateTime.parse(endString, dtf).toInstant(ZoneOffset.UTC).toEpochMilli();
      ds.trim(st, ed);
      StepExperiment se = new StepExperiment();
      se.runExperimentOnData(ds);
      double[] fitParams = se.getFitParams();
      assertEquals(366.97, 1./fitParams[0], 0.5);
      assertEquals(0.7196, fitParams[1], 0.0005);
    } catch (SeedFormatException | CodecException e) {
      e.printStackTrace();
      fail();
    }
  }

  @Test
  public void testMDJGoodCorner() {
    DataStore ds = new DataStore();
    String folder = "test-data/mdj-step/";
    String fname1 = "_BC0.512.seed";
    String fname2 = "00_BHZ.512.seed";
    try {
      ds.setBlock(0, folder + fname1);
      ds.setBlock(1, folder + fname2);
      ds.setEmbedResponse(1, "STS1T5_Q330HR");
      String startString = "2017-248T05:02:00.0";
      String endString = "2017-248T05:30:00.0";
      DateTimeFormatter dtf = DateTimeFormatter.ofPattern("uuuu-DDD'T'HH:mm:ss.S");
      long st = LocalDateTime.parse(startString, dtf).toInstant(ZoneOffset.UTC).toEpochMilli();
      long ed = LocalDateTime.parse(endString, dtf).toInstant(ZoneOffset.UTC).toEpochMilli();
      ds.trim(st, ed);
      StepExperiment se = new StepExperiment();
      se.runExperimentOnData(ds);
      double[] fitParams = se.getFitParams();
      assertEquals(364.5, 1./fitParams[0], 0.5);
      assertEquals(0.72, fitParams[1], 0.005);
    } catch (SeedFormatException | CodecException e) {
      e.printStackTrace();
      fail();
    }
  }

}
