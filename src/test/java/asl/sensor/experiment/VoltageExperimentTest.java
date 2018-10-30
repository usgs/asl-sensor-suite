package asl.sensor.experiment;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import asl.sensor.gui.ExperimentPanel;
import asl.sensor.input.DataStore;
import asl.sensor.input.InstrumentResponse;
import asl.sensor.test.TestUtils;
import edu.iris.dmc.seedcodec.CodecException;
import edu.sc.seis.seisFile.mseed.SeedFormatException;
import java.io.IOException;
import java.util.Calendar;
import org.junit.Test;

public class VoltageExperimentTest {

  private static final String folder = TestUtils.TEST_DATA_LOCATION + TestUtils.SUBPAGE;

  @Test
  public void passesGoodDataCorrectly() throws IOException {

    DataStore ds = new DataStore();

    String dataFolder = folder + "voltage-passes-AZI0/";
    String prefix = "AZI0_00_BH";
    String[] identifier = new String[3];
    identifier[0] = "1";
    identifier[1] = "2";
    identifier[2] = "Z";
    String extension = ".512.seed";

    InstrumentResponse ir = InstrumentResponse.loadEmbeddedResponse("STS2SGgen3_Q330HR");

    for (int i = 0; i < identifier.length; ++i) {
      String fName = dataFolder + prefix + identifier[i] + extension;
      try {
        ds.setBlock(i, fName);
        ds.setResponse(i, ir);
      } catch (SeedFormatException | CodecException | IOException e) {
        e.printStackTrace();
        fail();
      }
    }

    Calendar cCal = Calendar.getInstance(ExperimentPanel.DATE_TIME_FORMAT.get().getTimeZone());
    cCal.setTimeInMillis(ds.getBlock(0).getEndTime());
    System.out.println(cCal.get(Calendar.HOUR_OF_DAY));
    cCal.set(Calendar.HOUR_OF_DAY, 1);
    cCal.set(Calendar.MINUTE, 17);
    cCal.set(Calendar.SECOND, 30);
    long start = cCal.getTime().getTime();
    cCal.set(Calendar.MINUTE, 18);
    cCal.set(Calendar.SECOND, 33);
    long end = cCal.getTime().getTime();

    ds.trim(start, end);

    VoltageExperiment ve = new VoltageExperiment();
    assertTrue(ve.hasEnoughData(ds));
    ve.runExperimentOnData(ds);

    double[] gainValues = ve.getAllGainValues();
    double[] sensitivities = ve.getAllSensitivities();
    double[] percentDifferences = ve.getPercentDifferences();

    double gain = ir.getGain()[2];
    double[] expectedGains = {gain, gain, gain};
    double[] expectedSensitivities = {1677406.04, 1677602.62, 1677748.32};
    assertArrayEquals(expectedGains, gainValues,1E-15);
    assertArrayEquals(expectedSensitivities, sensitivities, 1E-2);
    for (double percentDifference : percentDifferences) {
      assertTrue(percentDifference < 0.05);
    }
  }

}
