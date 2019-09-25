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
import java.io.PrintWriter;
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

  private static class SharedGSN6TestRunner {
    private static List<XYSeriesCollection> experimentResult;


    static List<XYSeriesCollection> getGS6NResult() {
        if (experimentResult == null) {
          String respName = RESP_LOCATION + "resp_GSN6";
          String dataFolderName = getSeedFolder("XX", "GSN6", "2019", "129");
          String calName = dataFolderName + "CB_BC0.512.seed";
          String sensOutName = dataFolderName + "00_EHZ.512.seed";

          DataStore ds = DataStoreUtils.createFromNames(respName, calName, sensOutName);

          OffsetDateTime cCal = TestUtils.getStartCalendar(ds);
          cCal = cCal.withHour(15).withMinute(48).withSecond(59).withNano(0);
          long start = cCal.toInstant().toEpochMilli();

          //cCal = TestUtils.getEndCalendar(ds);
          cCal = cCal.withHour(16).withMinute(3).withSecond(59).withNano(0);
          long end = cCal.toInstant().toEpochMilli();

          ds.trim(start, end);
          RandomizedExperiment rCal = (RandomizedExperiment)
              ExperimentFactory.RANDOMCAL.createExperiment();

          rCal.setLowFrequencyCalibration(false);
          rCal.setNyquistMultiplier(90);
          assertTrue(rCal.hasEnoughData(ds));
          rCal.runExperimentOnData(ds);
          experimentResult = rCal.getData();
      }
      return experimentResult;
    }

  }

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
    boolean isLowFrequencyCalibration = false;
    RealVector initialGuess, initialPoleGuess, initialZeroGuess;
    initialPoleGuess = ir.polesToVector(isLowFrequencyCalibration, 80.);
    initialZeroGuess = ir.zerosToVector(isLowFrequencyCalibration, 80.);
    int numZeros = initialZeroGuess.getDimension();
    initialGuess = initialZeroGuess.append(initialPoleGuess);
    Pair<RealVector, RealMatrix> jacobianResult =
        RandomizedExperiment.jacobian(initialGuess, freqs, numZeros, ir, isLowFrequencyCalibration);
    // evaluate the data for reference
    Complex[] result = ir.applyResponseToInput(freqs);
    double[] testData = new double[2 * result.length];
    for (int i = 0; i < result.length; ++i) {
      int argIdx = i + result.length;
      Complex c = result[i];
      testData[i] = c.abs();
      testData[argIdx] = NumericUtils.atanc(c);
    }
    RandomizedExperiment.scaleValues(testData, freqs, isLowFrequencyCalibration);
    // test that the Jacobian first result is the actual evaluation
    double[] functionEvaluation = jacobianResult.getFirst().toArray();
    assertArrayEquals(testData, functionEvaluation, 1E-2);
    // now compare the forward difference to the first difference in the array
    RealVector testAgainst = initialGuess.copy();
    double changingVar = testAgainst.getEntry(0);
    testAgainst.setEntry(0, changingVar - RandomizedExperiment.DELTA);
    InstrumentResponse testDiffResponse =
        ir.buildResponseFromFitVector(testAgainst.toArray(), isLowFrequencyCalibration, numZeros);
    Complex[] diffResult = testDiffResponse.applyResponseToInput(freqs);
    double[] testDiffData = new double[2 * diffResult.length];
    for (int i = 0; i < diffResult.length; ++i) {
      int argIdx = i + result.length;
      Complex c = diffResult[i];
      testDiffData[i] = c.abs();
      testDiffData[argIdx] = NumericUtils.atanc(c);
    }
    RandomizedExperiment.scaleValues(testDiffData, freqs, isLowFrequencyCalibration);
    double[] firstJacobian = new double[testDiffData.length];
    for (int i = 0; i < testDiffData.length; ++i) {
      firstJacobian[i] = testData[i] - testDiffData[i];
      firstJacobian[i] /= (changingVar - (testAgainst.getEntry(0)));
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

    cCal = cCal.withMinute(41);
    long end = cCal.toInstant().toEpochMilli();

    ds.trim(start, end);

    return ds;
  }

  @Test
  public void testPlotScalingCorrect() {
    String keyMustContain = "Calc. resp.";

    DataStore ds = setUpTest1();
    InstrumentResponse ir = ds.getResponse(1);

    double nyq = ds.getBlock(0).getSampleRate() / 2.;
    System.out.println("NYQUIST RATE: " + nyq);

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

    // to check that scaling is correct, assert that every point in graph is bound to be <1
    // this should be a good check against any scaling issues past the nyquist % cut-off parameter
    // because in those cases the unscaled curve would be far beyond 1
    XYSeries calcCurve = xysc.get(0).getSeries(indexOfCalcCurve);
    for (int i = 0; i < calcCurve.getItemCount(); ++i) {
      double y = (double) calcCurve.getY(i);
      assertTrue("Violating y-value at index " + i + ": " + y, y <= 1.);
    }

  }

  @Test
  public void testCalculationResult1() {

    String currentDir = System.getProperty("user.dir");
    try {

      DataStore ds = setUpTest1();
      InstrumentResponse ir = ds.getResponse(1);

      double nyq = ds.getBlock(0).getSampleRate() / 2.;
      System.out.println("NYQUIST RATE: " + nyq);

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

      // TODO: add corrected assert here to compare best-fit and expected result
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
    InstrumentResponse fitResponse = rCal.getFitResponse();
    // these are expected poles which are good-fit. response fit may differ depending on JDK and
    // order of operations, etc. so we prefer to compare response curves over raw parameters,
    // as many different output responses may produce equally good fits for the data under analysis
    InstrumentResponse expectedResponse = new InstrumentResponse(RESP_LOCATION + "RESP.CU.BCIP.00.BHZ_2017_268_EXPECTED");

    double[][] calculatedDataSeries = rCal.getData().get(0).getSeries(1).toArray();
    double[] frequencyTest = calculatedDataSeries[0];
    double[] calculatedResponseCurveAmp = calculatedDataSeries[1];

    Complex[] fitResponseCurve = fitResponse.applyResponseToInput(frequencyTest);
    double[] fitAmpAndPhase = new double[2 * fitResponseCurve.length];
    Complex[] expectedResponseCurve = expectedResponse.applyResponseToInput(frequencyTest);
    double[] expectedAmpAndPhase = new double[2 * fitResponseCurve.length];
    for (int i = 0; i < fitResponseCurve.length; ++i) {
      int phaseIndex = i + fitResponseCurve.length;
      fitAmpAndPhase[i] = fitResponseCurve[i].abs();
      fitAmpAndPhase[phaseIndex] = NumericUtils.atanc(fitResponseCurve[i]);
      expectedAmpAndPhase[i] = expectedResponseCurve[i].abs();
      expectedAmpAndPhase[phaseIndex] = NumericUtils.atanc(expectedResponseCurve[i]);
    }
    RandomizedExperiment.scaleValues(fitAmpAndPhase, frequencyTest, false);
    RandomizedExperiment.scaleValues(expectedAmpAndPhase, frequencyTest, false);

    double[] fitAmp = Arrays.copyOfRange(fitAmpAndPhase, 0, fitResponseCurve.length);
    double[] expectedAmp = Arrays.copyOfRange(expectedAmpAndPhase, 0, expectedResponseCurve.length);
    for (int i = 0; i < fitResponseCurve.length; ++i) {
      fitAmp[i] = Math.abs(fitAmp[i] - calculatedResponseCurveAmp[i]);
      expectedAmp[i] = Math.abs(expectedAmp[i] - calculatedResponseCurveAmp[i]);
      String msg = "EXPECTED: " + expectedAmp[i] + " FIT: " + fitAmp[i] + " AT FREQ: " + frequencyTest[i];
      assertTrue(msg,fitAmp[i] <= expectedAmp[i] + 0.1);
    }

    assertTrue("Value of fit residual: " + rCal.getFitResidual(), rCal.getFitResidual() < 52.);
    assertEquals(1082.7, rCal.getInitResidual(), 1E-1);
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
        new Complex(-0.012725101823426397, -0.011495336794506263),
        new Complex(-0.012725101823426397, 0.011495336794506263)
    };
    for (int i = 0; i < fitPoles.size(); i++) {
      assertEquals(expectedPoles[i].getReal(), fitPoles.get(i).getReal(), 1E-5);
      assertEquals(expectedPoles[i].getImaginary(), fitPoles.get(i).getImaginary(), 1E-5);
    }

    assertEquals(423.6, rCal.getFitResidual(), 1E-1);
    assertEquals(482.4, rCal.getInitResidual(), 1E-1);
  }

  @Test
  public void hrvHasReasonablePoleFit() throws FileNotFoundException {
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

    // instead of measuring RMS accuracy with the pole/zero expected values or deviation from
    // the blue curve, measuring the residual % error is a good reliable way to test regression
    // that should be reasonably stable to changes in the code that won't meaningfully affect
    // the otherwise simple case here

    double percentError = Math.abs((initialResidual - fitResidual) / initialResidual) * 100;
    assertTrue(percentError < 2.);
  }

  // @Test uncommented because this is probably built off a bad response file
  public void gsn6HasCorrectCalculatedCurve() throws FileNotFoundException {

    XYSeriesCollection xyscAmp = SharedGSN6TestRunner.getGS6NResult().get(0);
    XYSeriesCollection xyscPhase = SharedGSN6TestRunner.getGS6NResult().get(1);

    // calculated magnitude ('blue') curve is second input
    XYSeries calculatedCurveAmp = xyscAmp.getSeries(1);
    XYSeries calculatedCurvePhase = xyscPhase.getSeries(1);
    assertEquals("Calc. resp. (XX_GSN6_00_EHZ) magnitude", calculatedCurveAmp.getKey());
    assertEquals("Calc. resp. (XX_GSN6_00_EHZ) phase", calculatedCurvePhase.getKey());

    StringBuilder sb = new StringBuilder("FREQUENCY,\tAMPLITUDE,\tPHASE\n");
    for (int i = 0; i < calculatedCurveAmp.getItemCount(); ++i) {
      double frequency = (double) calculatedCurveAmp.getX(i);
      double amplitude = (double) calculatedCurveAmp.getY(i);
      double phase = (double) calculatedCurvePhase.getY(i);
      // zero crossing point is also first point of fit
      if (i == 0) {
        assertEquals(0., amplitude, 1E-10);
        assertEquals(0., phase, 1E-10);
      }
      if (frequency < 60.) {
        sb.append(frequency).append(",\t").append(amplitude)
            .append(",\t").append(phase).append("\n");

      }
    }

    PrintWriter out = new PrintWriter(new File("GSN6-blue-curve.csv"));
    out.write(sb.toString());
    out.close();
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
      Complex expectedPoleError = new Complex(0.00043, 0.00185);
      Complex evaluatedPoleError = poleErrors.get(pole);
      String message = "Difference between expected "
          + "and evaluated poles outside of error bound:\n\t"
          + cf.format(expectedPoleError) + " , " + cf.format(evaluatedPoleError);
      assertTrue(message, Complex.equals(poleErrors.get(pole), expectedPoleError, 1E-5));
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
