package asl.sensor.experiment;

import static asl.sensor.test.TestUtils.RESP_LOCATION;
import static asl.sensor.test.TestUtils.getSeedFolder;
import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import asl.sensor.CalProcessingServer;
import asl.sensor.ExperimentFactory;
import asl.sensor.gui.RandomizedPanel;
import asl.sensor.input.DataStore;
import asl.sensor.input.DataStoreUtils;
import asl.sensor.output.CalResult;
import asl.sensor.test.TestUtils;
import asl.utils.NumericUtils;
import asl.utils.ReportingUtils;
import asl.utils.input.InstrumentResponse;
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
  public void testEvaluationOfJacobian() throws IOException {
    String fname = folder + "resp-parse/TST5_response.txt";
    InstrumentResponse ir;
    ir = new InstrumentResponse(fname);
    double[] freqs = new double[80];
    for (int i = 0; i < freqs.length; ++i) {
      freqs[i] = i;
    }
    final boolean lowFreq = false;
    RealVector initialGuess, initialPoleGuess, initialZeroGuess;
    initialPoleGuess = ir.polesToVector(lowFreq, Double.MAX_VALUE);
    initialZeroGuess = ir.zerosToVector(lowFreq, Double.MAX_VALUE);
    int numZeros = initialZeroGuess.getDimension();
    initialGuess = initialZeroGuess.append(initialPoleGuess);
    Pair<RealVector, RealMatrix> jacobianResult =
        RandomizedExperiment.jacobian(initialGuess, freqs, numZeros, ir, lowFreq);

    // evaluate the data for reference
    Complex[] result = ir.applyResponseToInputUnscaled(freqs);
    double[] testData = new double[2 * result.length];
    for (int i = 0; i < result.length; ++i) {
      int argIdx = i + result.length;
      Complex c = result[i];
      testData[i] = c.abs();
      testData[argIdx] = NumericUtils.atanc(c);
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
    InstrumentResponse testDiffResponse =
        ir.buildResponseFromFitVector(testAgainst.toArray(), lowFreq, numZeros);
    Complex[] diffResult = testDiffResponse.applyResponseToInputUnscaled(freqs);
    double[] testDiffData = new double[2 * diffResult.length];
    for (int i = 0; i < diffResult.length; ++i) {
      int argIdx = i + result.length;
      Complex c = diffResult[i];
      testDiffData[i] = c.abs();
      testDiffData[argIdx] = NumericUtils.atanc(c);
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

  @Test
  public void testCalculationResult1() {

    String currentDir = System.getProperty("user.dir");
    try {

      DataStore ds = setUpTest1();
      InstrumentResponse ir = ds.getResponse(1);

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

      double pctDiff =
          Math.abs(100 * (bestResid - expectedResid) / bestResid);

      if (pctDiff > 15) {
        System.out.println(rCal.getFitPoles());
        System.out.println(rCal.getInitialPoles());
        System.out.println(bestResid + ", " + expectedResid);
      }

      //assertTrue("PCT DIFF EXPECTED <15%, GOT " + pctDiff, pctDiff < 15);

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
      System.out.println("Output result has been written");

    } catch (IOException e) {
      e.printStackTrace();
      fail();
    }
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
      CalResult rCalResult = cps.runRand(calInFile, sensorOutFile, respFile,
          false, start, end, true);
      System.out.println(Arrays.toString(rCalResult.getNumerMap().get("Best_fit_poles")));
      System.out.println(Arrays.toString(rCalResult.getNumerMap().get("Best_fit_zeros")));
    } catch (IOException | SeedFormatException | CodecException e) {
      e.printStackTrace();
      fail();
    }
  }

  @Test
  public void runExperiment_BCIP_HFCalibration() throws IOException {
    String respName = RESP_LOCATION + "RESP.CU.BCIP.00.BHZ_2017_268";
    String dataFolderName = getSeedFolder("CU", "BCIP", "2017", "268");
    String calName = dataFolderName + "CB_BC0.512.seed";
    String sensOutName = dataFolderName + "00_EHZ.512.seed";

    DataStore ds = DataStoreUtils.createFromNames(respName, calName, sensOutName);
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
    InstrumentResponse fitResponse = rCal.getFitResponse();


    double[][] calculatedDataSeries = rCal.getData().get(0).getSeries(1).toArray();
    double[] frequencyTest = calculatedDataSeries[0];
    double[] calculatedResponseCurveAmp = calculatedDataSeries[1];

    Complex[] fitResponseCurve = fitResponse.applyResponseToInput(frequencyTest);
    double[] fitAmpAndPhase = new double[2 * fitResponseCurve.length];
    double[] expectedAmpAndPhase = new double[2 * fitResponseCurve.length];
    for (int i = 0; i < fitResponseCurve.length; ++i) {
      int phaseIndex = i + fitResponseCurve.length;
      fitAmpAndPhase[i] = fitResponseCurve[i].abs();
      fitAmpAndPhase[phaseIndex] = NumericUtils.atanc(fitResponseCurve[i]);
    }
    RandomizedExperiment.scaleValues(fitAmpAndPhase, frequencyTest, false);
    RandomizedExperiment.scaleValues(expectedAmpAndPhase, frequencyTest, false);

    double[] fitAmp = Arrays.copyOfRange(fitAmpAndPhase, 0, fitResponseCurve.length);
    double ampRMS = 0.;
    for (int i = 0; i < fitResponseCurve.length; ++i) {
      ampRMS += Math.pow(fitAmp[i] - calculatedResponseCurveAmp[i], 2);
    }

    ampRMS = Math.sqrt(ampRMS / fitResponseCurve.length);
    assertTrue(ampRMS < 1.5);

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

    RandomizedExperiment rCal = (RandomizedExperiment)
        ExperimentFactory.RANDOMCAL.createExperiment();

    rCal.setLowFrequencyCalibration(true);

    assertTrue(rCal.hasEnoughData(ds));
    rCal.runExperimentOnData(ds);

    List<Complex> fitPoles = rCal.getFitPoles();
    Complex[] expectedPoles = {
        new Complex(-0.01247, -0.01179),
        new Complex(-0.01247,  0.01179)
    };
    for (int i = 0; i < fitPoles.size(); i++) {
      assertEquals(expectedPoles[i].getReal(), fitPoles.get(i).getReal(), 5E-4);
      assertEquals(expectedPoles[i].getImaginary(), fitPoles.get(i).getImaginary(), 5E-4);
    }

    assertEquals(3.96975, rCal.getFitResidual(), 1E-3);
    assertEquals(4.30781, rCal.getInitResidual(), 1E-3);
  }

  @Test
  public void hrvHasReasonableRMSMeasure() throws FileNotFoundException {
    String respName = RESP_LOCATION + "RESP.IU.HRV.00.BHZ";
    String dataFolderName = getSeedFolder("IU", "HRV", "2018", "192");
    String calName = dataFolderName + "CB_BC0.512.seed";
    String sensOutName = dataFolderName + "00_EHZ.512.seed";

    DataStore ds = DataStoreUtils.createFromNames(respName, calName, sensOutName);
    Complex expectedFitPole = new Complex(-37.39683437804999, -82.25199086188465);
    InstrumentResponse ir = new InstrumentResponse(ds.getResponse(1));
    List<Complex> poles = ir.getPoles();
    // expectedFitPole = poles.get(4);
    poles.set(4, expectedFitPole);
    poles.set(5, expectedFitPole.conjugate());
    ir.setPoles(poles);
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

  // @Test TODO: evaluate something in all of this
  public void verifyRespCurveExpectation() {
    String respName1 = RESP_LOCATION + "RESP.IU.SNZO.00.BHZ";
    String respName2 = RESP_LOCATION + "NEW.RESP.IU.SNZO.00.BHZ";
    String dataFolderName = getSeedFolder("IU", "SNZO", "2019", "086");
    String calName = dataFolderName + "CB_BC0.512.seed";
    String sensOutName = dataFolderName + "00_EHZ.512.seed";
    String startTime = "2019,086,15:31:00";
    DateTimeFormatter dateTimeFormatter =
        DateTimeFormatter.ofPattern("uuuu,DDD,HH:mm:ss").withZone(ZoneOffset.UTC);

    DataStore ds = DataStoreUtils.createFromNames(respName1, calName, sensOutName);
    long startCal = ZonedDateTime.parse(startTime, dateTimeFormatter).toInstant().toEpochMilli();
    String endTime = "2019,086,15:56:00";
    long endCal = ZonedDateTime.parse(endTime, dateTimeFormatter).toInstant().toEpochMilli();
    ds.trim(startCal, endCal);

    RandomizedExperiment rCal = new RandomizedExperiment();
    rCal.setLowFrequencyCalibration(false);
    rCal.setNyquistMultiplier(0.6);
    assertTrue(rCal.hasEnoughData(ds));
    rCal.runExperimentOnData(ds);

    InstrumentResponse ir = rCal.getFitResponse();
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
      Complex expectedPoleError = new Complex(0.014546479, 0.0085197934);
      Complex evaluatedPoleError = poleErrors.get(pole);
      String message = "Difference between expected "
          + "and evaluated poles outside of error bound:\n\t"
          + cf.format(expectedPoleError) + " , " + cf.format(evaluatedPoleError);
      assertTrue(message, Complex.equals(poleErrors.get(pole), expectedPoleError, 2E-5));
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
}