package asl.sensor.experiment;

import static asl.sensor.experiment.RandomizedExperiment.DEFAULT_NYQUIST_PERCENT_LIMIT;
import static asl.sensor.experiment.RandomizedExperiment.getResponseCorrection;
import static asl.sensor.test.TestUtils.RESP_LOCATION;
import static asl.sensor.test.TestUtils.getSeedFolder;
import static asl.utils.NumericUtils.TAU;
import static asl.utils.NumericUtils.atanc;
import static asl.utils.response.ResponseBuilders.deepCopyResponse;
import static asl.utils.response.ResponseParser.loadEmbeddedResponse;
import static asl.utils.response.ResponseParser.parseResponse;
import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import asl.sensor.CalProcessingServer;
import asl.sensor.ExperimentFactory;
import asl.sensor.gui.RandomizedPanel;
import asl.sensor.input.DataStore;
import asl.sensor.input.DataStoreUtils;
import asl.sensor.output.CalResult;
import asl.sensor.test.TestUtils;
import asl.utils.ReportingUtils;
import asl.utils.response.ChannelMetadata;
import asl.utils.response.ChannelMetadata.ResponseStageException;
import asl.utils.response.PolesZeros.Pole;
import asl.utils.response.ResponseUnits;
import asl.utils.response.ResponseUnits.ResolutionType;
import asl.utils.response.ResponseUnits.SensorType;
import edu.iris.dmc.seedcodec.CodecException;
import edu.sc.seis.seisFile.mseed.SeedFormatException;
import java.awt.Font;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.text.DecimalFormat;
import java.time.OffsetDateTime;
import java.time.ZoneOffset;
import java.time.ZonedDateTime;
import java.time.format.DateTimeFormatter;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.complex.ComplexFormat;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;
import org.apache.pdfbox.pdmodel.PDDocument;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.LogarithmicAxis;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.junit.Test;

public class RandomizedExperimentTest {

  private static final String folder = TestUtils.TEST_DATA_LOCATION + TestUtils.SUBPAGE;
  private static final String testRespName =
      folder + "random-high-32+70i/RESP.XX.NS088..BHZ.STS1.360.2400";


  @Test
  public void testSampleRateMismatchError() {
    String respName = "STS2gen3_Q330HR";
    String dataFolderName = getSeedFolder("IU", "FUNA", "2019", "073");
    String calName = dataFolderName + "_BC0.512.seed";
    String sensOutName = dataFolderName + "00_BHZ.512.seed";

    DataStore ds = DataStoreUtils.createFromNamesEmbedResp(respName, calName, sensOutName);

    String startTime = "2019,073,18:51";
    DateTimeFormatter dateTimeFormatter =
        DateTimeFormatter.ofPattern("uuuu,DDD,HH:mm").withZone(ZoneOffset.UTC);
    long startCal = ZonedDateTime.parse(startTime, dateTimeFormatter).toInstant().toEpochMilli();
    String endTime = "2019,073,22:51";
    long endCal = ZonedDateTime.parse(endTime, dateTimeFormatter).toInstant().toEpochMilli();

    ds.trim(startCal, endCal);

    RandomizedExperiment randomExp = new RandomizedExperiment();
    randomExp.setLowFrequencyCalibration(true);
    randomExp.setCapactiveCalibration(false);
    randomExp.runExperimentOnData(ds);

    // now, do we get a null pointer exception?

  }

  @Test
  public void testCorrectionResponseLoadsCorrectly() throws IOException {
    RandomizedExperiment rExp = new RandomizedExperiment();
    rExp.setCorrectionResponse(SensorType.TR120);
  }

  @Test
  public void testCorrectionResponse_WontLoadNonTrillium() throws IOException {
    RandomizedExperiment rExp = new RandomizedExperiment();
    rExp.setCorrectionResponse(SensorType.STS2gen3);
    assertNull(rExp.getCorrectionToApply());
    rExp.setCorrectionResponse(null);
    assertNull(rExp.getCorrectionToApply());
  }

  @Test
  public void testEvaluationOfJacobian() throws IOException, ResponseStageException {
    String fname = folder + "resp-parse/TST5_response.txt";
    ChannelMetadata ir;
    ir = parseResponse(fname);
    double[] freqs = new double[80];
    for (int i = 0; i < freqs.length; ++i) {
      freqs[i] = i;
    }
    final boolean lowFreq = false;
    RealVector initialGuess, initialPoleGuess, initialZeroGuess;
    initialPoleGuess = ir.getPoleZeroStage().polesToVector(lowFreq, Double.MAX_VALUE);
    initialZeroGuess = ir.getPoleZeroStage().zerosToVector(lowFreq, Double.MAX_VALUE);
    int numZeros = initialZeroGuess.getDimension();
    initialGuess = initialZeroGuess.append(initialPoleGuess);
    Complex[] corrections = getResponseCorrection(freqs, null);
    Pair<RealVector, RealMatrix> jacobianResult =
        RandomizedExperiment.jacobian(initialGuess, freqs, numZeros, ir, lowFreq, corrections);

    // evaluate the data for reference, limiting to the PZ stage
    Complex[] result = ir.applyResponseToInputUnscaled(freqs, 1, 2);
    double[] testData = new double[2 * result.length];
    for (int i = 0; i < result.length; ++i) {
      int argIdx = i + result.length;
      Complex c = result[i];
      testData[i] = c.abs();
      testData[argIdx] = atanc(c);
    }
    RandomizedExperiment.scaleValues(testData, freqs, lowFreq);

    // test that the Jacobian first result is the actual evaluation
    double[] functionEvaluation = jacobianResult.getFirst().toArray();
    
    assertArrayEquals(testData, functionEvaluation, 1E-2);
    // now compare the forward difference to the first difference in the array
    RealVector testAgainst = initialGuess.copy();
    double changingVar = testAgainst.getEntry(0);
    double jacobianDiff = 100 * Math.ulp(changingVar);
    testAgainst.setEntry(0, changingVar + jacobianDiff);
    ChannelMetadata testDiffResponse =
        ir.buildResponseFromFitVector(testAgainst.toArray(), lowFreq, numZeros);
    // again, we restrict to PZ stage (1)
    Complex[] diffResult = testDiffResponse.applyResponseToInputUnscaled(freqs, 1, 2);
    double[] testDiffData = new double[2 * diffResult.length];
    for (int i = 0; i < diffResult.length; ++i) {
      int argIdx = i + result.length;
      Complex c = diffResult[i];
      testDiffData[i] = c.abs();
      testDiffData[argIdx] = atanc(c);
    }
    RandomizedExperiment.scaleValues(testDiffData, freqs, lowFreq);
    double[] firstJacobian = new double[testDiffData.length];
    for (int i = 0; i < testDiffData.length; ++i) {
      firstJacobian[i] = testDiffData[i] - testData[i];
      firstJacobian[i] /= jacobianDiff;
    }
    double[] testFirstJacobianAgainst = jacobianResult.getSecond().getColumnVector(0).toArray();
    assertArrayEquals(testFirstJacobianAgainst, firstJacobian, 1E-3);
  }

  private DataStore setUpTest1() {

    String dataFolderName = folder + "random-high-32+70i/";
    String calName = dataFolderName + "_EC0.512.seed";
    String sensOutName = dataFolderName + "00_EHZ.512.seed";

    DataStore ds = DataStoreUtils.createFromNames(testRespName, calName, sensOutName);
    OffsetDateTime cCal = TestUtils.getStartCalendar(ds);

    cCal = cCal.withMinute(36);
    cCal = cCal.withSecond(0);
    long start = cCal.toInstant().toEpochMilli();

    cCal = cCal.withMinute(49); // was 41 -- changed to 49 to allow sensor output to settle
    long end = cCal.toInstant().toEpochMilli();

    ds.trim(start, end);

    return ds;
  }

  @Test
  public void testPlotScalingCorrect() {
    String keyMustContain = "Calc. resp.";

    DataStore ds = setUpTest1();
    RandomizedExperiment rCal = (RandomizedExperiment)
        ExperimentFactory.RANDOMCAL.createExperiment();
    rCal.setLowFrequencyCalibration(false);
    assertTrue(rCal.hasEnoughData(ds));
    rCal.runExperimentOnData(ds);

    List<XYSeriesCollection> xysc = rCal.getData();
    int indexOfCalcCurve = 0;
    for (int i = 0; i < xysc.get(0).getSeriesCount(); ++i) {
      String key = (String) xysc.get(0).getSeriesKey(i);
      if (key.startsWith(keyMustContain)) {
        indexOfCalcCurve = i;
        break;
      }
    }

    // to check that scaling is correct, assert data near the normalization is close to 0
    XYSeries calcCurve = xysc.get(0).getSeries(indexOfCalcCurve);

    for (int i = 0; i < calcCurve.getItemCount(); ++i) {
      double x = (double) calcCurve.getX(i);
      if (x > 1 && x < 2.5) {
        double y = (double) calcCurve.getY(i);
        assertTrue("Curve value above expected normalization at index " + i  + " (freq. "
            + x + "): Got amplitude value of " + y, Math.abs(y) < 0.15);
      }
    }


  }

  // @Test
  public void testCalculationResult1() throws IOException, ResponseStageException {

    String currentDir = System.getProperty("user.dir");
    DataStore ds = setUpTest1();
    ChannelMetadata ir = ds.getResponse(1);

    RandomizedExperiment rCal = (RandomizedExperiment)
        ExperimentFactory.RANDOMCAL.createExperiment();

    rCal.setLowFrequencyCalibration(false);

    assertTrue(rCal.hasEnoughData(ds));
    rCal.runExperimentOnData(ds);

    double bestResid = rCal.getFitResidual();

    int width = 1280;
    int height = 960;

    List<XYSeriesCollection> xysc = rCal.getData();


    String[] yAxisTitles = new String[]{"Resp(f), dB", "Angle / TAU"};
    JFreeChart[] jfcl = new JFreeChart[yAxisTitles.length];

    String xAxisTitle = "Frequency (Hz)";
    NumberAxis xAxis = new LogarithmicAxis(xAxisTitle);
    Font bold = xAxis.getLabelFont().deriveFont(Font.BOLD);
    xAxis.setLabelFont(bold);

    StringBuilder sb = new StringBuilder();

    String[] resultString = rCal.getInsetStrings();
    for (String resultPart : resultString) {
      sb.append(resultPart);
      sb.append('\n');
    }
    sb.append(rCal.getFormattedDateRange());
    sb.append('\n');
    sb.append("Input files:\n");
    sb.append(ds.getBlock(0).getName());
    sb.append(" (calibration)\n");
    sb.append(ds.getBlock(1).getName());
    sb.append(" (sensor output)\n");
    sb.append("Response file used:\n");
    sb.append(ds.getResponse(1).getName());
    sb.append("\n \n");

    String page1 = sb.toString();

    String[] addtlPages = (RandomizedPanel.getAdditionalReportPages(rCal));
    // technically 'page 2' but really second part of first dataset report
    // and I'm too lazy to rename everything to reflect that
    String page1Part2 = addtlPages[0];

    sb = new StringBuilder();

    // expected best fit params, for debugging
    sb.append("BELOW RESULTS FOR EXPECTED BEST FIT (YELLOW CURVE)\n");
    double[] expectedParams = new double[]{-3.580104E+1, +7.122400E+1};
    ir = ir.buildResponseFromFitVector(expectedParams, false, 0);
    ir.setName("Best-fit params");
    ds.setResponse(1, ir);
    rCal.runExperimentOnData(ds);

    // residual from other code's best-fit parameters
    // compare to best-fit residual and assume difference is < 5%
    double expectedResid = rCal.getInitResidual();

    // add initial curve from expected fit params to report
    XYSeries expectedInitialCurve = rCal.getData().get(0).getSeries(0);
    xysc.get(0).addSeries(expectedInitialCurve);
    XYSeries expectedInitialAngle = rCal.getData().get(1).getSeries(0);
    xysc.get(1).addSeries(expectedInitialAngle);

    resultString = rCal.getInsetStrings();
    for (String resultPart : resultString) {
      sb.append(resultPart);
      sb.append('\n');
    }

    for (int i = 0; i < jfcl.length; ++i) {

      jfcl[i] = ChartFactory.createXYLineChart(
          ExperimentFactory.RANDOMCAL.getName(),
          xAxisTitle,
          yAxisTitles[i],
          xysc.get(i),
          PlotOrientation.VERTICAL,
          true,
          false,
          false);

      XYPlot xyp = jfcl[i].getXYPlot();

      //xyp.clearAnnotations();
      //xyp.addAnnotation(xyt);

      xyp.setDomainAxis(xAxis);
    }

    String page2 = sb.toString();

    PDDocument pdf = new PDDocument();
    ReportingUtils.chartsToPDFPage(width, height, pdf, jfcl);
    ReportingUtils.textListToPDFPages(pdf, page1, page1Part2, page2);

    String testResultFolder = currentDir + "/testResultImages/";
    File dir = new File(testResultFolder);
    if (!dir.exists()) {
      dir.mkdir();
    }

    String testResult =
        testResultFolder + "Random-Calib-Test-1.pdf";
    pdf.save(new File(testResult));
    pdf.close();
  }

  @Test
  public void testCalServerEntryMethods() {
    String testFolder = folder + "test-crashed-on-cal/";
    String calInFile = testFolder + "_BC0.512.seed";
    String sensorOutFile = testFolder + "00_BHZ.512.seed";
    String respFile = testFolder + "RESP.US.MVCO.00.BHZ";
    String start = "2018-01-30T07:55:00+00:00";
    String end = "2018-01-30T11:55:00+00:00";
    try {
      CalProcessingServer cps = new CalProcessingServer();
      // empty string because this calibration does not require a trillium correction
      CalResult rCalResult = cps.runRand(calInFile, sensorOutFile, respFile,
          false, start, end, true,
          DEFAULT_NYQUIST_PERCENT_LIMIT, "");
    } catch (IOException | SeedFormatException | CodecException | ResponseStageException e) {
      e.printStackTrace();
      fail();
    }
  }

  @Test
  public void runExperiment_BCIP_HFCalibration() {
    // using an embedded resp because the actual resp file used for this test just isn't a very
    // good fit to begin with. Ideally we switch this to data that a resp fits to well before
    // solving, and verify that the RMS is low
    String respName =
        ResponseUnits.getFilenameFromComponents(SensorType.STS2_5, ResolutionType.HIGH);
    String dataFolderName = getSeedFolder("CU", "BCIP", "2017", "268");
    String calName = dataFolderName + "CB_BC0.512.seed";
    String sensOutName = dataFolderName + "00_EHZ.512.seed";

    DataStore ds = DataStoreUtils.createFromNamesEmbedResp(respName, calName, sensOutName);

    OffsetDateTime cCal = TestUtils.getStartCalendar(ds);
    cCal = cCal.withHour(18).withMinute(49).withSecond(0).withNano(0);
    long start = cCal.toInstant().toEpochMilli();
    cCal = cCal.withHour(19).withMinute(4);
    long end = cCal.toInstant().toEpochMilli();
    ds.trim(start, end);

    RandomizedExperiment rCal = (RandomizedExperiment)
        ExperimentFactory.RANDOMCAL.createExperiment();

    rCal.setLowFrequencyCalibration(false);
    rCal.setNyquistMultiplier(.8);

    assertTrue(rCal.hasEnoughData(ds));
    rCal.runExperimentOnData(ds);
    double initResidual = rCal.getInitResidual();
    double fitResidual = rCal.getFitResidual();
    ChannelMetadata fitResponse = rCal.getFitResponse();

    double[][] calculatedDataSeries = rCal.getData().get(0).getSeries(1).toArray();
    double[] frequencyTest = calculatedDataSeries[0];
    double[] calculatedResponseCurveAmp = calculatedDataSeries[1];

    Complex[] fitResponseCurve = fitResponse.applyResponseToInput(frequencyTest);
    double[] fitAmpAndPhase = new double[2 * fitResponseCurve.length];
    for (int i = 0; i < fitResponseCurve.length; ++i) {
      int phaseIndex = i + fitResponseCurve.length;
      fitAmpAndPhase[i] = fitResponseCurve[i].abs();
      fitAmpAndPhase[phaseIndex] = atanc(fitResponseCurve[i]);
    }
    RandomizedExperiment.scaleValues(fitAmpAndPhase, frequencyTest, false);
    assertEquals(calculatedResponseCurveAmp.length, fitAmpAndPhase.length / 2);
    assertEquals(calculatedResponseCurveAmp.length, fitResponseCurve.length);

    double[] fitAmp = Arrays.copyOfRange(fitAmpAndPhase, 0, fitResponseCurve.length);
    double ampRMS = 0.;
    for (int i = 0; i < fitAmp.length; ++i) {
      ampRMS += Math.pow(fitAmp[i] - calculatedResponseCurveAmp[i], 2);
    }

    ampRMS = Math.sqrt(ampRMS / fitResponseCurve.length);
    assertTrue("ampRMS larger than expected (1.5), is " + ampRMS, ampRMS < 1.5);

    assertTrue(fitResidual <= initResidual);
    // assertTrue("Percent error not less than 30: " + pctError, pctError < 30);
  }

  @Test
  public void runExperiment_KIEV_LFCalibration() {
    String respName = RESP_LOCATION + "RESP.IU.KIEV.00.BH1";
    String dataFolderName = getSeedFolder("IU", "KIEV", "2018", "044");
    String calName = dataFolderName + "_BC0.512.seed";
    String sensOutName = dataFolderName + "00_BH1.512.seed";

    DataStore ds = DataStoreUtils.createFromNames(respName, calName, sensOutName);

    dataFolderName = getSeedFolder("IU", "KIEV", "2018", "045");
    calName = dataFolderName + "_BC0.512.seed";
    sensOutName = dataFolderName + "00_BH1.512.seed";

    ds = DataStoreUtils.appendFromNames(ds, calName, sensOutName);

    OffsetDateTime cCal = TestUtils.getStartCalendar(ds);
    cCal = cCal.withHour(23).withMinute(37).withSecond(0).withNano(0);
    long start = cCal.toInstant().toEpochMilli();

    cCal = TestUtils.getEndCalendar(ds);
    cCal = cCal.withHour(7).withMinute(37);
    long end = cCal.toInstant().toEpochMilli();

    ds.trim(start, end);

    Complex[] expectedPoles = {
        new Complex(-0.01247, -0.01174),
        new Complex(-0.01247,  0.01174)
    };

    RandomizedExperiment rCal = (RandomizedExperiment)
        ExperimentFactory.RANDOMCAL.createExperiment();

    rCal.setLowFrequencyCalibration(true);

    assertTrue(rCal.hasEnoughData(ds));
    rCal.runExperimentOnData(ds);
    List<Complex> fitPoles = rCal.getFitPoles();

    for (int i = 0; i < fitPoles.size(); i++) {
      assertEquals(expectedPoles[i].getReal(), fitPoles.get(i).getReal(), 5E-4);
      assertEquals(expectedPoles[i].getImaginary(), fitPoles.get(i).getImaginary(), 5E-4);
    }

    assertTrue(rCal.getFitResidual() < rCal.getInitResidual());
  }

  @Test
  public void rssdCheckAgainstNaNErrorEstimates() {

    String resp = ResponseUnits.getFilenameFromComponents(SensorType.STS6, ResolutionType.HIGH);

    String dataFolderName1 = getSeedFolder("IU", "RSSD", "2019", "319");
    String dataFolderName2 = getSeedFolder("IU", "RSSD", "2019", "320");
    String calFilename = "_BC0.512.seed";
    String calFirstDay = dataFolderName1 + calFilename;
    String calSecondDay = dataFolderName2 + calFilename;
    String sensorFilename = "00_BHZ.512.seed";
    String sensorFirstDay = dataFolderName1 + sensorFilename;
    String sensorSecondDay = dataFolderName2 + sensorFilename;

    DataStore ds = DataStoreUtils.createFromNamesEmbedResp(resp, calFirstDay, sensorFirstDay);
    DataStoreUtils.appendFromNames(ds, calSecondDay, sensorSecondDay);

    long start = TestUtils.timeStringToEpochMilli("2019-319T22:00:48.2")+91;
    long end = TestUtils.timeStringToEpochMilli("2019-320T06:08:50.2")+32;
    ds.trim(start, end);

    RandomizedExperiment rExp = new RandomizedExperiment();
    rExp.setLowFrequencyCalibration(true);
    rExp.runExperimentOnData(ds);

    for (Complex error : rExp.getPoleErrors().values()) {
      assertNotEquals(error.getReal(), Double.NaN);
    }
    for (Complex error : rExp.getZeroErrors().values()) {
      assertNotEquals(error.getReal(), Double.NaN);
    }
  }

  @Test
  public void hrvHasReasonableRMSMeasure() throws FileNotFoundException, ResponseStageException {
    String respName = RESP_LOCATION + "RESP.IU.HRV.00.BHZ";
    String dataFolderName = getSeedFolder("IU", "HRV", "2018", "192");
    String calName = dataFolderName + "CB_BC0.512.seed";
    String sensOutName = dataFolderName + "00_EHZ.512.seed";

    DataStore ds = DataStoreUtils.createFromNames(respName, calName, sensOutName);
    Complex expectedFitPole = new Complex(-37.39683437804999, -82.25199086188465);
    ChannelMetadata ir = deepCopyResponse(ds.getResponse(1));
    List<Pole> poles = ir.getPoleZeroStage().getPoleDoubleList();
    // expectedFitPole = poles.get(4);
    poles.set(4, new Pole(expectedFitPole));
    poles.set(5, new Pole(expectedFitPole.conjugate()));
    ir.getPoleZeroStage().setPoles(poles);
    ds.setResponse(1, ir);

    RandomizedExperiment rCal = new RandomizedExperiment();
    rCal.setLowFrequencyCalibration(false);
    rCal.setNyquistMultiplier(0.3);
    assertTrue(rCal.hasEnoughData(ds));
    rCal.runExperimentOnData(ds);

    double initialResidual = rCal.getInitResidual();
    double fitResidual = rCal.getFitResidual();

    assertTrue("Residual value over 25: " + fitResidual,fitResidual < 25.);
  }

  @Test // TODO: evaluate something in all of this
  public void verifyRespCurveBehavior() throws IOException {
    String respName1 = RESP_LOCATION + "RESP.IU.SNZO.00.BHZ";
    String dataFolderName = getSeedFolder("IU", "SNZO", "2019", "086");
    String calName = dataFolderName + "CB_BC0.512.seed";
    String sensOutName = dataFolderName + "00_EHZ.512.seed";
    String startTime = "2019,086,15:31:00";
    DateTimeFormatter dateTimeFormatter =
        DateTimeFormatter.ofPattern("uuuu,DDD,HH:mm:ss").withZone(ZoneOffset.UTC);

    DataStore ds = DataStoreUtils.createFromNames(respName1, calName, sensOutName);
    ds.setResponse(1, loadEmbeddedResponse(SensorType.TR360, ResolutionType.HIGH));
    assertFalse(ds.responseIsSet(0));
    long startCal = ZonedDateTime.parse(startTime, dateTimeFormatter).toInstant().toEpochMilli();
    String endTime = "2019,086,15:56:00";
    long endCal = ZonedDateTime.parse(endTime, dateTimeFormatter).toInstant().toEpochMilli();
    ds.trim(startCal, endCal);

    RandomizedExperiment rCal = new RandomizedExperiment();
    rCal.setLowFrequencyCalibration(false);
    rCal.setNyquistMultiplier(0.6);
    assertTrue(rCal.hasEnoughData(ds));
    rCal.runExperimentOnData(ds);

    // second entry in data list is phase series for plots, first one is initial repsonse curve
    // second index in array is the series of y-values
    double[] results = rCal.getData().get(1).getSeries(0).toArray()[1];
    double[] allZeros = new double[results.length]; // array values initialized to zero
    // if this assert fails, there is an issue with the calculation of response curves in
    // the random cal solver
    assertFalse(Arrays.equals(allZeros, results));
  }

  @Test
  public void hrvHasReasonablePoleFit() throws FileNotFoundException {
    String respName = RESP_LOCATION + "RESP.IU.HRV.00.BHZ";
    String dataFolderName = getSeedFolder("IU", "HRV", "2018", "192");
    String calName = dataFolderName + "CB_BC0.512.seed";
    String sensOutName = dataFolderName + "00_EHZ.512.seed";

    DataStore ds = DataStoreUtils.createFromNames(respName, calName, sensOutName);
    Complex expectedFitPole = new Complex(-34.92592022858792, -82.25199086188465);

    RandomizedExperiment rCal = new RandomizedExperiment();
    rCal.setLowFrequencyCalibration(false);
    rCal.setNyquistMultiplier(0.3);
    assertTrue(rCal.hasEnoughData(ds));
    rCal.runExperimentOnData(ds);

    List<Complex> initialPoles = rCal.getFitPoles();

    assertEquals(2, initialPoles.size());
    assertEquals(expectedFitPole.getReal(), initialPoles.get(0).getReal(), 0.5);
    assertEquals(expectedFitPole.getImaginary(), initialPoles.get(0).getImaginary(), 0.5);
  }

  @Test
  public void kievHasCorrectError() {
    String respName = RESP_LOCATION + "RESP.IU.KIEV.00.BH1";
    String dataFolderName = getSeedFolder("IU", "KIEV", "2018", "044");
    String calName = dataFolderName + "_BC0.512.seed";
    String sensOutName = dataFolderName + "00_BH1.512.seed";

    DataStore ds = DataStoreUtils.createFromNames(respName, calName, sensOutName);

    dataFolderName = getSeedFolder("IU", "KIEV", "2018", "045");
    calName = dataFolderName + "_BC0.512.seed";
    sensOutName = dataFolderName + "00_BH1.512.seed";

    ds = DataStoreUtils.appendFromNames(ds, calName, sensOutName);

    OffsetDateTime cCal = TestUtils.getStartCalendar(ds);
    cCal = cCal.withHour(23).withMinute(37).withSecond(0).withNano(0);
    long start = cCal.toInstant().toEpochMilli();

    cCal = TestUtils.getEndCalendar(ds);
    cCal = cCal.withHour(7).withMinute(37);
    long end = cCal.toInstant().toEpochMilli();

    ds.trim(start, end);

    RandomizedExperiment rCal = (RandomizedExperiment)
        ExperimentFactory.RANDOMCAL.createExperiment();

    rCal.setLowFrequencyCalibration(true);

    assertTrue(rCal.hasEnoughData(ds));
    rCal.runExperimentOnData(ds);

    Map<Complex, Complex> poleErrors = rCal.getPoleErrors();
    Map<Complex, Complex> zeroErrors = rCal.getZeroErrors();

    ComplexFormat cf = new ComplexFormat(new DecimalFormat("#.##########"));

    assertEquals(0, zeroErrors.size());
    assertEquals(2, poleErrors.size());

    Complex[] evaluatedPoles = poleErrors.keySet().toArray(new Complex[]{});

    for (Complex pole : evaluatedPoles) {
      // these values are clearly dependent on weighting scheme for calculated calibration curve
      Complex expectedPoleError = new Complex(0.0000520969, 0.0000590451);
      Complex evaluatedPoleError = poleErrors.get(pole);
      String message = "Difference between expected "
          + "and evaluated error values outside of precision bound:\n\t"
          + cf.format(expectedPoleError) + " , " + cf.format(evaluatedPoleError);
      double error = 1E-5;
      assertTrue(message,
          expectedPoleError.getReal() + error >= evaluatedPoleError.getReal());
      assertTrue(message,
          expectedPoleError.getImaginary() + error >= evaluatedPoleError.getImaginary());
    }

  }

  @Test
  public void hasEnoughData_missingInputData() {
    String respName = RESP_LOCATION + "RESP.CU.BCIP.00.BHZ_2017_268";
    String dataFolderName = getSeedFolder("CU", "BCIP", "2017", "268");
    String sensOutName = dataFolderName + "00_EHZ.512.seed";

    DataStore ds = DataStoreUtils.createFromNames(respName, null, sensOutName);
    RandomizedExperiment rCal = (RandomizedExperiment)
        ExperimentFactory.RANDOMCAL.createExperiment();
    assertFalse(rCal.hasEnoughData(ds));
  }

  @Test
  public void hasEnoughData_missingOutputData() {
    String respName = RESP_LOCATION + "RESP.CU.BCIP.00.BHZ_2017_268";
    String dataFolderName = getSeedFolder("CU", "BCIP", "2017", "268");
    String calName = dataFolderName + "CB_BC0.512.seed";

    DataStore ds = DataStoreUtils.createFromNames(respName, calName, null);
    RandomizedExperiment rCal = (RandomizedExperiment)
        ExperimentFactory.RANDOMCAL.createExperiment();
    assertFalse(rCal.hasEnoughData(ds));
  }

  @Test
  public void hasEnoughData_missingBothData() {
    String respName = RESP_LOCATION + "RESP.CU.BCIP.00.BHZ_2017_268";

    DataStore ds = DataStoreUtils.createFromNames(respName, null, null);
    RandomizedExperiment rCal = (RandomizedExperiment)
        ExperimentFactory.RANDOMCAL.createExperiment();
    assertFalse(rCal.hasEnoughData(ds));
  }

  @Test
  public void hasEnoughData_hasEnoughData() {
    String respName = RESP_LOCATION + "RESP.CU.BCIP.00.BHZ_2017_268";
    String dataFolderName = getSeedFolder("CU", "BCIP", "2017", "268");
    String calName = dataFolderName + "CB_BC0.512.seed";
    String sensOutName = dataFolderName + "00_EHZ.512.seed";

    DataStore ds = DataStoreUtils.createFromNames(respName, calName, sensOutName);
    RandomizedExperiment rCal = (RandomizedExperiment)
        ExperimentFactory.RANDOMCAL.createExperiment();
    assertTrue(rCal.hasEnoughData(ds));
  }

  @Test
  public void verifyThatCurveTrimmingIsCorrect() throws IOException {
    double[] freqs = new double[]{1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    ChannelMetadata resp = loadEmbeddedResponse(SensorType.STS2gen3, ResolutionType.HIGH);
    Complex[] curve = resp.applyResponseToInput(freqs);
    double[] ampAndPhase = new double[2 * curve.length];
    for (int i = 0; i < curve.length; ++i) {
      ampAndPhase[i] = curve[i].abs();
      ampAndPhase[curve.length+i] = atanc(curve[i]);
    }
    for (int j = 0; j < freqs.length; ++j) {
      final double[] freqsExceptOne = new double[freqs.length - 1];
      System.arraycopy(freqs, 0, freqsExceptOne, 0, j);
      if (j + 1 < freqs.length) {
        System.arraycopy(freqs, j + 1, freqsExceptOne, j, freqsExceptOne.length - j);
      }
      Complex[] trimmedCurve = resp.applyResponseToInput(freqsExceptOne);
      double[] expectedTrimmedAmpAndPhase = new double[trimmedCurve.length * 2];
      for (int i = 0; i < trimmedCurve.length; ++i) {
        expectedTrimmedAmpAndPhase[i] = trimmedCurve[i].abs();
        expectedTrimmedAmpAndPhase[i+trimmedCurve.length] = atanc(trimmedCurve[i]);
      }
      // a copy and paste from randomizedexperiment with the variable names changed
      final double[] observedCurve = new double[ampAndPhase.length - 2];
      int trimOffset = ampAndPhase.length / 2; // where the phase component starts
      // copy the first j points to this array (from 0 to j-1)
      System.arraycopy(ampAndPhase, 0, observedCurve, 0, j);
      // now copy from from j+1 to where the phase component starts
      System.arraycopy(ampAndPhase, j + 1, observedCurve, j,
          trimOffset - j - 1);
      // now the first j phase components (note offset for destination due to missing phase term)
      System.arraycopy(ampAndPhase, trimOffset, observedCurve,
          trimOffset - 1, j);
      System.arraycopy(ampAndPhase, trimOffset + j + 1, observedCurve,
          trimOffset - 1 + j, trimOffset - j - 1);

      assertArrayEquals(expectedTrimmedAmpAndPhase, observedCurve, 0);
    }
  }

  @Test
  public void testSyntheticLowFreq333() throws IOException {
    String testLocation = folder + "synthetic-lowfreq-cals/";
    String calSignalFile = testLocation + "_BC0.FAKE.333.mseed";
    String calOutputFile = testLocation + "00_BHZ.FAKE.333.mseed";
    String testRespName = "STS6_Q330HR";
    DataStore ds =
        DataStoreUtils.createFromNamesEmbedResp(testRespName, calSignalFile, calOutputFile);

    RandomizedExperiment rExp = new RandomizedExperiment();
    rExp.setLowFrequencyCalibration(true);
    rExp.runExperimentOnData(ds);
    List<Complex> poles = rExp.getFitPoles();
    // assert that there are two poles, complex conjugates
    assertEquals(2, poles.size());
    assertNotEquals(0., poles.get(0).getImaginary());
    assertEquals(poles.get(0).getReal(), poles.get(1).getReal(), 0.);
    assertEquals(poles.get(0).getImaginary(), -1 * poles.get(1).getImaginary(), 0.);
    assertEquals(333,  TAU / poles.get(0).abs(), 1);
  }

  @Test
  public void testSyntheticLowFreq125() throws IOException {
    String testLocation = folder + "synthetic-lowfreq-cals/";
    String calSignalFile = testLocation + "_BC0.FAKE.125.mseed";
    String calOutputFile = testLocation + "00_BHZ.FAKE.125.mseed";
    String testRespName = "STS2gen3_Q330HR";
    DataStore ds =
        DataStoreUtils.createFromNamesEmbedResp(testRespName, calSignalFile, calOutputFile);

    RandomizedExperiment rExp = new RandomizedExperiment();
    rExp.setLowFrequencyCalibration(true);
    rExp.runExperimentOnData(ds);
    List<Complex> poles = rExp.getFitPoles();
    // assert that there are two poles, complex conjugates
    assertEquals(2, poles.size());
    assertNotEquals(0., poles.get(0).getImaginary());
    assertEquals(poles.get(0).getReal(), poles.get(1).getReal(), 0.);
    assertEquals(poles.get(0).getImaginary(), -1 * poles.get(1).getImaginary(), 0.);
    assertEquals(125,  TAU / poles.get(0).abs(), 1);
  }

  @Test
  public void testSyntheticHighFreq333() throws IOException {
    String testLocation = folder + "synthetic-lowfreq-cals/";
    String calSignalFile = testLocation + "_BC0.FAKE.333.mseed";
    String calOutputFile = testLocation + "00_BHZ.FAKE.333.mseed";
    String testRespName = testLocation + "RESP.XX.FAKE.333.BHZ";
    DataStore ds =
        DataStoreUtils.createFromNames(testRespName, calSignalFile, calOutputFile);

    RandomizedExperiment rExp = new RandomizedExperiment();
    rExp.setLowFrequencyCalibration(false);
    rExp.runExperimentOnData(ds);

    // the included response is intended to be as close to the optimum fit for this function as
    // possible. so to test a valid result we'll make sure the results agree to a minimum level of
    // precision (the results will not be exactly identical under most likely circumstances)
    // which we will set to be the third significant figure
    List<Complex> poles = rExp.getFitPoles();
    // also, as some poles are close to optimum, we allow for them not all being changed
    assertTrue(poles.size() <= 5);
    for (int i = 0; i < poles.size(); ++i) {
      double orderOfMagnitude = Math.floor(Math.log10(poles.get(i).abs()));
      double delta = 2 * Math.pow(10, orderOfMagnitude - 3);
      assertEquals(rExp.getInitialPoles().get(i).abs(), poles.get(i).abs(), delta);
    }

    List<Complex> zeros = rExp.getFitZeros();
    assertTrue(zeros.size() <= 6);
    for (int i = 0; i < zeros.size(); ++i) {
      double orderOfMagnitude = Math.ceil(Math.log10(zeros.get(i).abs()));
      double delta = 2 * Math.pow(10, orderOfMagnitude - 3);
      assertEquals(rExp.getInitialZeros().get(i).abs(), zeros.get(i).abs(), delta);
    }
  }

  @Test
  public void testSyntheticHighFreq125() throws IOException {
    String testLocation = folder + "synthetic-lowfreq-cals/";
    String calSignalFile = testLocation + "_BC0.FAKE.125.mseed";
    String calOutputFile = testLocation + "00_BHZ.FAKE.125.mseed";
    String testRespName = testLocation + "RESP.XX.FAKE.125.BHZ";
    DataStore ds =
        DataStoreUtils.createFromNames(testRespName, calSignalFile, calOutputFile);

    RandomizedExperiment rExp = new RandomizedExperiment();
    rExp.setLowFrequencyCalibration(false);
    rExp.setNyquistMultiplier(0.5);
    rExp.runExperimentOnData(ds);

    // the included response is intended to be as close to the optimum fit for this function as
    // possible. so to test a valid result we'll make sure the results agree to a minimum level of
    // precision (the results will not be exactly identical under most likely circumstances)
    // which we will set to be the third significant figure

    List<Complex> poles = rExp.getFitPoles();
    // 9 poles have a frequency above 1 hz, but the last 3 are above 2000 Hz
    // that's far too high to include in the fit, so we don't adjust them at all
    // leaving us with only 6 poles that the solver is going to adjust
    assertEquals(6, poles.size());
    for (int i = 0; i < poles.size(); ++i) {
      double orderOfMagnitude = Math.floor(Math.log10(poles.get(i).abs()));
      double delta = Math.pow(10, orderOfMagnitude - 3);
      assertEquals(rExp.getInitialPoles().get(i).abs(), poles.get(i).abs(), delta);
    }

    List<Complex> zeros = rExp.getFitZeros();
    assertEquals(4, zeros.size());
    for (int i = 0; i < zeros.size(); ++i) {
      double orderOfMagnitude = Math.ceil(Math.log10(zeros.get(i).abs()));
      double delta = Math.pow(10, orderOfMagnitude - 3);
      assertEquals(rExp.getInitialZeros().get(i).abs(), zeros.get(i).abs(), delta);
    }
  }

  @Test
  public void assureThatResponseCorrectionIsGood() throws IOException {
    String respName = RESP_LOCATION + "RESP.IU.SJG.00.BHZ";
    String dataFolderName = getSeedFolder("IU", "SJG", "2019", "121");
    String sensOutName = dataFolderName + "00_EHZ.512.seed";
    String calInName = dataFolderName + "CB_BC0.512.seed";
    DataStore ds =
        DataStoreUtils.createFromNames(respName, calInName, sensOutName);

    String startTime = "2019,121,17:22:00";
    DateTimeFormatter dateTimeFormatter =
        DateTimeFormatter.ofPattern("uuuu,DDD,HH:mm:ss").withZone(ZoneOffset.UTC);
    long startCal = ZonedDateTime.parse(startTime, dateTimeFormatter).toInstant().toEpochMilli();
    String endTime = "2019,121,17:47:30";
    long endCal = ZonedDateTime.parse(endTime, dateTimeFormatter).toInstant().toEpochMilli();
    ds.trim(startCal, endCal);

    RandomizedExperiment randomExp = new RandomizedExperiment();
    randomExp.setCorrectionResponse(null);
    randomExp.setLowFrequencyCalibration(false);
    randomExp.setCapactiveCalibration(false);
    randomExp.runExperimentOnData(ds);

    double uncorrectedResidual = randomExp.getInitResidual();
    double uncorrectedFitResidual = randomExp.getFitResidual();
    randomExp.setCorrectionResponse(SensorType.TR360);
    randomExp.runExperimentOnData(ds);
    double correctedResidual = randomExp.getInitResidual();
    double correctedFitResidual = randomExp.getFitResidual();
    assertTrue(correctedResidual < uncorrectedResidual);
    // the correction CANNOT be compensated for by the solver on its own
    // the uncorrected fit is WORSE than the corrected initial response
    assertTrue(correctedResidual < uncorrectedFitResidual);
    assertTrue(correctedFitResidual < uncorrectedFitResidual);
  }
}