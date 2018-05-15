package asl.sensor.experiment;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import asl.sensor.input.DataStore;
import asl.sensor.test.TestUtils;
import edu.iris.dmc.seedcodec.CodecException;
import edu.sc.seis.seisFile.mseed.SeedFormatException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import org.apache.commons.math3.complex.Complex;
import org.junit.Test;

public class StepTest {

  public static String folder = TestUtils.TEST_DATA_LOCATION + TestUtils.SUBPAGE;

  @Test
  public void testKievGoodCorner() {
    DataStore ds = new DataStore();
    String testFolder = folder + "kiev-step/";
    String fname1 = "_BC0.512.seed";
    String fname2 = "00_BHZ.512.seed";
    try {
      ds.setBlock(0, testFolder + fname1);
      ds.setBlock(1, testFolder + fname2);
      ds.setEmbedResponse(1, "resps/STS1T5_Q330HR");
      String startString = "2018-038T15:25:00.0";
      String endString = "2018-038T16:00:00.0";
      long st = TestUtils.timeStringToEpochMilli(startString);
      long ed = TestUtils.timeStringToEpochMilli(endString);
      ds.trim(st, ed);
      StepExperiment se = new StepExperiment();
      se.runExperimentOnData(ds);
      double[] fitParams = se.getFitParams();
      assertEquals(366.97, 1. / fitParams[0], 0.5);
      assertEquals(0.7196, fitParams[1], 0.0005);
    } catch (IOException | SeedFormatException | CodecException e) {
      e.printStackTrace();
      fail();
    }
  }

  @Test
  public void testKonoGoodCorner() {
    DataStore ds = new DataStore();
    String testFolder = folder + "kono-step/";
    String fname1 = "_BC1.512.seed";
    String fname2 = "10_BHZ.512.seed";
    try {
      ds.setBlock(0, testFolder + fname1);
      ds.setBlock(1, testFolder + fname2);
      ds.setEmbedResponse(1, "resps/STS2gen3_Q330HR");
      String startString = "2018-037T20:02:00.0";
      String endString = "2018-037T20:17:00.0";
      long st = TestUtils.timeStringToEpochMilli(startString);
      long ed = TestUtils.timeStringToEpochMilli(endString);
      ds.trim(st, ed);
      StepExperiment se = new StepExperiment();
      se.runExperimentOnData(ds);
      double[] fitParams = se.getFitParams();
      assertEquals(120.00, 1. / fitParams[0], 0.5);
      assertEquals(0.7035, fitParams[1], 0.0005);
    } catch (IOException | SeedFormatException | CodecException e) {
      e.printStackTrace();
      fail();
    }
  }

  @Test
  public void testMDJGoodCorner() {
    DataStore ds = new DataStore();
    String testFolder = folder + "mdj-step/";
    String fname1 = "_BC0.512.seed";
    String fname2 = "00_BHZ.512.seed";
    try {
      ds.setBlock(0, testFolder + fname1);
      ds.setBlock(1, testFolder + fname2);
      ds.setEmbedResponse(1, "resps/STS1T5_Q330HR");
      String startString = "2017-248T05:00:00.0";
      String endString = "2017-248T05:30:00.0";
      long st = TestUtils.timeStringToEpochMilli(startString);
      long ed = TestUtils.timeStringToEpochMilli(endString);
      ds.trim(st, ed);
      StepExperiment se = new StepExperiment();
      se.runExperimentOnData(ds);
      double[] fitParams = se.getFitParams();
      assertEquals(367.5, 1. / fitParams[0], 0.5);
      assertEquals(0.715, fitParams[1], 0.005);
    } catch (IOException | SeedFormatException | CodecException e) {
      e.printStackTrace();
      fail();
    }
  }

  @Test
  public void testINCNAnomaly() {
    DataStore ds = new DataStore();
    String testFolder = folder + "incn-step/";
    String fname1 = "_BC0.512.seed";
    String fname2 = "00_BHZ.512.seed";
    try {
      ds.setBlock(0, testFolder + fname1);
      ds.setBlock(1, testFolder + fname2);
      ds.setEmbedResponse(1, "resps/TR360_Q330HR");
      String startString = "2018-102T22:22:00.0";
      String endString = "2018-102T23:00:00.0";
      long st = TestUtils.timeStringToEpochMilli(startString);
      long ed = TestUtils.timeStringToEpochMilli(endString);
      ds.trim(st, ed);
      StepExperiment se = new StepExperiment();
      se.runExperimentOnData(ds);
      double[] fitParams = se.getFitParams();
      assertEquals(371.3, 1. / fitParams[0], 0.5);
      assertEquals(0.719, fitParams[1], 0.005);
    } catch (IOException | SeedFormatException | CodecException e) {
      e.printStackTrace();
      fail();
    }
  }

  @Test
  public void testWaterLevelCalc() {
    // water level calc intended to invert value unless it is 0, then set it to 0 instead
    // note that the values are scaled WRT to the maximum value of the data arrays
    List<Complex> data = new ArrayList<>();
    data.add(new Complex(1, 2));
    data.add(Complex.ZERO);
    data.add(new Complex(2, 2));
    data.add(new Complex(1, -1));
    data.add(Complex.ZERO);
    data.add(new Complex(100, -3));
    data.add(new Complex(5, -10));
    data.add(new Complex(20));
    Complex[] arr = data.toArray(new Complex[]{});
    data = new ArrayList<>();
    data.add(new Complex(.2, -.4));
    data.add(Complex.ZERO);
    data.add(new Complex(.25, -.25));
    data.add(new Complex(.5, .5));
    data.add(Complex.ZERO);
    data.add(new Complex(0.00999101, 2.99730243E-4));
    data.add(new Complex(0.04, 0.08));
    data.add(new Complex(0.05));
    arr = StepExperiment.setWaterLevel(arr);
    for (int i = 0; i < arr.length; ++i) {
      assertTrue(Complex.equals(arr[i], data.get(i), 1E-7));
    }
  }

}
