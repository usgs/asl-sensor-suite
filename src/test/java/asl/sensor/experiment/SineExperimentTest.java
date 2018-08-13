package asl.sensor.experiment;

import static org.junit.Assert.assertEquals;

import asl.sensor.input.DataStore;
import asl.sensor.test.TestUtils;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

public class SineExperimentTest {

  private static final String folder = TestUtils.TEST_DATA_LOCATION + TestUtils.SUBPAGE;

  private static DataStore dataStore;

  @Before
  public void setUpData() throws Exception{
    dataStore = new DataStore();
    String calFname = folder + "sine-test/" + "_BC0.512.seed";
    String outFname = folder + "sine-test/" + "00_BHZ.512.seed";
    dataStore.setBlock(0, calFname);
    dataStore.setBlock(1, outFname);
    String startTimeString = "2015-166T20:23:08.5";
    String endTimeString = "2015-166T21:02:46.7";
    long start = TestUtils.timeStringToEpochMilli(startTimeString);
    long end = TestUtils.timeStringToEpochMilli(endTimeString);
    dataStore.trim(start, end);
  }

  @After
  public void cleanUpData() {
    dataStore = null;
  }


  @Test
  public void testCalFrequency() {
    SineExperiment sexp = new SineExperiment();
    sexp.runExperimentOnData(dataStore);
    double ratio = sexp.getCalAmplitude() / sexp.getOutAmplitude();
    assertEquals(0.02768, ratio, 1E-3);
    double freq = sexp.getEstSineFreq();
    assertEquals(250, freq, 2.);
  }

}
