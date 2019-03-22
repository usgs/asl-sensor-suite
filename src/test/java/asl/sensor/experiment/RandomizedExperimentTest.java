package asl.sensor.experiment;

import static asl.sensor.test.TestUtils.RESP_LOCATION;
import static asl.sensor.test.TestUtils.getSeedFolder;
import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import asl.sensor.output.CalResult;
import java.awt.Font;
import java.io.File;
import java.io.IOException;
import java.time.Instant;
import java.time.OffsetDateTime;
import java.time.ZoneOffset;
import java.time.ZonedDateTime;
import java.time.format.DateTimeFormatter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Calendar;
import java.util.List;
import java.util.Random;
import org.apache.commons.math3.complex.Complex;
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
import asl.sensor.CalProcessingServer;
import asl.sensor.ExperimentFactory;
import asl.sensor.gui.RandomizedPanel;
import asl.sensor.input.DataStore;
import asl.sensor.input.DataStoreUtils;
import asl.sensor.input.InstrumentResponse;
import asl.sensor.test.TestUtils;
import asl.sensor.utils.NumericUtils;
import asl.sensor.utils.ReportingUtils;
import edu.iris.dmc.seedcodec.CodecException;
import edu.sc.seis.seisFile.mseed.SeedFormatException;

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

  @Test
  public void responseCorrectConvertedToVectorHighFreq() throws Exception {
    String fname = folder + "resp-parse/TST5_response.txt";
    InstrumentResponse ir;
    ir = new InstrumentResponse(fname);
    List<Complex> poles = new ArrayList<>(ir.getPoles());
    // using an unnecessarily high nyquist rate here
    RealVector high = ir.polesToVector(false, 1E8);

    int complexIndex = 2; // start at second pole
    int vectorIndex = 0;

    while (vectorIndex < high.getDimension()) {
      // return current index
      double real = high.getEntry(vectorIndex++);
      double imag = high.getEntry(vectorIndex++);

      double poleImag = poles.get(complexIndex).getImaginary();

      assertEquals(real, poles.get(complexIndex).getReal(), 0.0);
      assertEquals(imag, poleImag, 0.0);

      if (poleImag != 0) {
        // complex conjugate case
        ++complexIndex;
        assertEquals(real, poles.get(complexIndex).getReal(), 0.0);
        assertEquals(imag, -poles.get(complexIndex).getImaginary(), 0.0);
      }
      ++complexIndex;
    }
  }

  @Test
  public void responseCorrectlyConvertedToVectorLowFreq() {
    String fname = folder + "resp-parse/TST5_response.txt";
    InstrumentResponse ir;
    try {

      ir = new InstrumentResponse(fname);
      List<Complex> poles = new ArrayList<>(ir.getPoles());
      // again, use a very high nyquist rate
      RealVector low = ir.polesToVector(true, 1E8);

      // only test lower two poless
      assertEquals(low.getEntry(0), poles.get(0).getReal(), 0.0);
      assertEquals(low.getEntry(1), poles.get(0).getImaginary(), 0.0);

      assertEquals(low.getEntry(0), poles.get(1).getReal(), 0.0);
      assertEquals(low.getEntry(1), -poles.get(1).getImaginary(), 0.0);

    } catch (IOException e) {
      fail();
      e.printStackTrace();
    }
  }

  @Test
  public void responseSetCorrectlyHighFreq() {
    String fname = folder + "resp-parse/TST5_response.txt";
    InstrumentResponse ir;

    try {
      ir = new InstrumentResponse(fname);

      List<Complex> poles = new ArrayList<>(ir.getPoles());
      List<Complex> replacements = new ArrayList<>();

      int start = 2;
      if (poles.get(0).getImaginary() == 0) {
        start = 1;
      }

      for (int i = start; i < poles.size(); ++i) {
        if (poles.get(i).getImaginary() == 0) {
          Complex c = poles.get(i);
          replacements.add(c.subtract(1));
          int next = i + 1;
          while (next < poles.size() && poles.get(next).equals(c)) {
            ++next; // skip duplicates
          }
        } else {
          Complex c = poles.get(i);
          c = c.subtract(new Complex(1, 1));
          replacements.add(c);
          ++i;
        }
      }

      //System.out.println(poles);
      //System.out.println(replacements);

      double[] newPoles = new double[replacements.size() * 2];
      for (int i = 0; i < newPoles.length; i += 2) {
        int poleIdx = i / 2;
        Complex c = replacements.get(poleIdx);
        newPoles[i] = c.getReal();
        newPoles[i + 1] = c.getImaginary();
      }

      InstrumentResponse ir2 =
          ir.buildResponseFromFitVector(newPoles, false, 0);

      List<Complex> testList = ir2.getPoles();

      System.out.println(poles);
      System.out.println(testList);
      // System.out.println(poles2);

      int offsetIdx = 0;
      for (int i = 0; i < poles.size(); ++i) {
        if (i < start) {
          assertEquals(poles.get(i), testList.get(i));
        } else {
          Complex c = replacements.get(offsetIdx);
          assertEquals(testList.get(i), c);
          if (poles.get(i).getImaginary() != 0) {
            Complex c1 = new Complex(1, 1);
            assertEquals(poles.get(i), c.add(c1));
            ++i;
            Complex c2 = new Complex(1, -1);
            assertEquals(testList.get(i), c.conjugate());
            assertEquals(poles.get(i), c.conjugate().add(c2));
          } else {
            assertEquals(poles.get(i), c.add(1));
          }
          ++offsetIdx;
        }
      }

    } catch (IOException e) {
      e.printStackTrace();
    }

  }

  @Test
  public void responseSetCorrectlyLowFreq() {
    String fname = folder + "resp-parse/TST5_response.txt";
    InstrumentResponse ir;
    try {
      ir = new InstrumentResponse(fname);
      List<Complex> poles = new ArrayList<>(ir.getPoles());

      double[] newPoles = new double[2];
      newPoles[0] = 0.;
      newPoles[1] = 1.;

      Complex c = new Complex(newPoles[0], newPoles[1]);

      InstrumentResponse ir2 =
          ir.buildResponseFromFitVector(newPoles, true, 0);
      List<Complex> poles2 = ir2.getPoles();

      List<Complex> testList = new ArrayList<>(poles);
      testList.set(0, c);
      testList.set(1, c.conjugate());

      for (int i = 0; i < poles.size(); ++i) {
        if (i < 2) {
          assertNotEquals(poles.get(i), poles2.get(i));
          assertEquals(poles2.get(i), testList.get(i));
        }
      }


    } catch (IOException e) {
      fail();
      e.printStackTrace();
    }

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
    assertEquals(1082.7313334829698, rCal.getInitResidual(), 1E-7);
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

    assertEquals(423.7415521942539, rCal.getFitResidual(), 1E-6);
    assertEquals(482.45559437599235, rCal.getInitResidual(), 1E-7);
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
