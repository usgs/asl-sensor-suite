package asl.sensor.output;

import static asl.sensor.test.TestUtils.RESP_LOCATION;
import static asl.sensor.test.TestUtils.getSeedFolder;
import static org.junit.Assert.assertEquals;

import asl.sensor.CalProcessingServer;
import asl.sensor.test.TestUtils;
import asl.utils.input.InstrumentResponse;
import edu.iris.dmc.seedcodec.CodecException;
import edu.sc.seis.seisFile.mseed.SeedFormatException;
import java.io.IOException;
import org.apache.commons.math3.complex.Complex;
import org.junit.Test;

public class CalServerTest {

  private static final String folder = TestUtils.TEST_DATA_LOCATION + TestUtils.SUBPAGE;

  @Test
  public void testRandomNoOOBException() throws SeedFormatException, CodecException, IOException {
    String respName = "STS2gen3_Q330HR";
    String dataFolderName = getSeedFolder("IU", "FUNA", "2019", "073");
    String calName = dataFolderName + "_BC0.512.seed";
    String sensOutName = dataFolderName + "00_BHZ.512.seed";
    String startTime = "2019-03-14T18:51:00Z";
    String endTime = "2019-03-14T22:51:00Z";

    CalProcessingServer server = new CalProcessingServer();
    server.runRand(calName, sensOutName, respName, true, startTime, endTime, true);

    // now, do we get a null pointer exception?
  }

  @Test
  public void testRandomResults() throws SeedFormatException, CodecException, IOException {
    String respName = RESP_LOCATION + "RESP.IU.KIEV.00.BH1";
    String dataFolderName = getSeedFolder("IU", "KIEV", "2018", "044");
    String calName = dataFolderName + "_BC0.512.seed";
    String sensOutName = dataFolderName + "00_BH1.512.seed";

    String startDateTime = "2018-02-13T23:37:00Z";
    String endDateTime = "2018-02-14T07:37:00Z";

    dataFolderName = getSeedFolder("IU", "KIEV", "2018", "045");
    String calName2 = dataFolderName + "_BC0.512.seed";
    String sensOutName2 = dataFolderName + "00_BH1.512.seed";

    CalProcessingServer server = new CalProcessingServer();
    CalResult result = server.runRand(calName, calName2, sensOutName, sensOutName2, respName, false,
        startDateTime, endDateTime, true);

    double[] initPolesDoubles = result.getNumerMap().get("Initial_poles");
    double[] fitPolesDoubles = result.getNumerMap().get("Best_fit_poles");

    Complex[] initPoles = new Complex[initPolesDoubles.length/2];
    Complex[] fitPoles = new Complex[initPoles.length];
    for (int i = 0; i < fitPoles.length; ++i) {
      int polesRawIndex = i * 2;
      double initRealPart = initPolesDoubles[polesRawIndex];
      double initImagPart = initPolesDoubles[polesRawIndex + 1];
      double fitRealPart = fitPolesDoubles[polesRawIndex];
      double fitImagPart = fitPolesDoubles[polesRawIndex + 1];
      initPoles[i] = new Complex(initRealPart, initImagPart);
      fitPoles[i] = new Complex(fitRealPart, fitImagPart);
    }

    Complex[] expectedInitPoles =
        new InstrumentResponse(respName).getPoles().toArray(new Complex[]{});

    Complex[] expectedFitPoles = {
        new Complex(-0.012540, -0.011495),
        new Complex(-0.012540,  0.011495),
        expectedInitPoles[2],
        expectedInitPoles[3],
    };

    assertEquals(expectedInitPoles.length, initPoles.length);
    assertEquals(expectedFitPoles.length, fitPoles.length);

    for (int i = 0; i < fitPoles.length; i++) {
      assertEquals(expectedInitPoles[i].getReal(), initPoles[i].getReal(), 1.5E-4);
      assertEquals(expectedInitPoles[i].getImaginary(), initPoles[i].getImaginary(), 1.5E-4);

      assertEquals(expectedFitPoles[i].getReal(), fitPoles[i].getReal(), 1.5E-4);
      assertEquals(expectedFitPoles[i].getImaginary(), fitPoles[i].getImaginary(), 1.5E-4);
    }

  }

  @Test
  public void testStepResults() throws IOException, CodecException, SeedFormatException {
    String testFolder = folder + "kiev-step/";
    String calInputFile = testFolder + "_BC0.512.seed";
    String calOutputFile = testFolder +  "00_BHZ.512.seed";
    String resp = "STS1T5_Q330HR";
    final boolean EMBED = true;

    String startString = "2018-02-07T15:20:00+00:00";
    String endString = "2018-02-07T15:59:00+00:00";

    CalProcessingServer server = new CalProcessingServer();
    CalResult result = server.runStep(calInputFile, calOutputFile,
            resp, EMBED, startString, endString);

    double fitCorner = result.numerMap.get("Fit_corner")[0];
    double fitDamping = result.numerMap.get("Fit_damping")[0];

    assertEquals(366.97, 1./fitCorner, 0.5);
    assertEquals(0.7196, fitDamping, 0.0005);
  }

  @Test
  public void testSineResults() throws IOException, CodecException, SeedFormatException {
    String calInputFile = folder + "sine-test/" + "_BC0.512.seed";
    String calOutputFile = folder + "sine-test/" + "00_BHZ.512.seed";
    String startTimeString = "2015-06-15T20:23:08+00:00";
    String endTimeString = "2015-06-15T21:02:46+00:00";

    CalProcessingServer server = new CalProcessingServer();
    CalResult result = server.runSine(calInputFile, calOutputFile, startTimeString, endTimeString);

    double ratio = result.numerMap.get("Calibration_to_output_ratio")[0];
    double freq = result.numerMap.get("Estimated_signal_frequency")[0];

    assertEquals(0.02768, ratio, 1E-3);
    assertEquals(250, freq, 2.);
  }
}
