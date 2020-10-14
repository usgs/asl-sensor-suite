package asl.sensor.experiment;

import static asl.sensor.experiment.OrientedSineExperiment.AMP_CUTOFF;
import static asl.sensor.experiment.OrientedSineExperiment.FREQ_CUTOFF;
import static java.lang.Math.abs;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import asl.sensor.input.DataStore;
import asl.sensor.test.TestUtils;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

public class OrientedSineExperimentTest {

  private static final String folder = TestUtils.TEST_DATA_LOCATION + TestUtils.SUBPAGE;
  private static DataStore dataStore;

  @Before
  public void setUpData() throws Exception{
    dataStore = new DataStore();
    String northName = folder + "oriented-sine-test/" + "00_HH1.512.seed";
    String eastName = folder + "oriented-sine-test/" + "00_HH2.512.seed";
    String vertName = folder + "oriented-sine-test/" + "00_HHZ.512.seed";
    dataStore.setBlock(0, northName);
    dataStore.setBlock(1, eastName);
    dataStore.setBlock(2, vertName);
    String startTimeString = "2020-281T20:52:00.0";
    String endTimeString = "2020-281T20:54:59.0";
    // these will be better for amplitude
    //String startTimeString = "2015-166T20:26:00.0";
    //String endTimeString = "2015-166T21:02:30.0";
    long start = TestUtils.timeStringToEpochMilli(startTimeString);
    long end = TestUtils.timeStringToEpochMilli(endTimeString);
    dataStore.trim(start, end);
  }

  @After
  public void cleanUpData() {
    dataStore = null;
  }


  @Test
  public void testFrequency() {
    OrientedSineExperiment sexp = new OrientedSineExperiment();
    sexp.setDoRotation(true);
    sexp.setRotationTrillium(false);
    sexp.runExperimentOnData(dataStore);
    double[] frequencies = sexp.getFrequencies();
    double meanFreq = sexp.getMeanFrequency();
    double[] freqErrors = sexp.getFrequencyErrors();
    for (int i = 0; i < frequencies.length; ++i) {
      double error = abs(frequencies[i] - meanFreq) / Math.abs(meanFreq);
      assertEquals(freqErrors[i], error, 1E-4);
      assertTrue(freqErrors[i] < FREQ_CUTOFF);
    }
  }

  @Test
  public void testAmplitude() {
    OrientedSineExperiment sexp = new OrientedSineExperiment();
    sexp.setDoRotation(true);
    sexp.setRotationTrillium(false);
    sexp.runExperimentOnData(dataStore);
    double[] amplitudes = sexp.getAmplitudes();
    double meanAmp = sexp.getMeanAmplitude();
    double[] ampErrors = sexp.getAmplitudeErrors();

    for (int i = 0; i < amplitudes.length; ++i) {
      double error = abs(amplitudes[i] - meanAmp) / Math.abs(meanAmp);
      assertEquals(ampErrors[i], error, 1E-4);
      assertTrue("Error value at index " + i + " greater than cutoff: " + ampErrors[i],
          ampErrors[i] < AMP_CUTOFF);
    }
  }

}
