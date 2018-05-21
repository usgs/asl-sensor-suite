package asl.sensor.experiment;

import static asl.sensor.test.TestUtils.RESP_LOCATION;
import static asl.sensor.test.TestUtils.getSeedFolder;
import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;
import java.awt.Font;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.time.OffsetDateTime;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
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
import asl.sensor.CalProcessingServer;
import asl.sensor.CalProcessingServer.RandData;
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
      //System.out.println(testList);
      int offsetIdx = 0;
      for (int i = 0; i < poles.size(); ++i) {
        if (i < start) {
          assertTrue(poles.get(i).equals(testList.get(i)));
        } else {
          Complex c = replacements.get(offsetIdx);
          assertTrue(testList.get(i).equals(c));
          if (poles.get(i).getImaginary() != 0) {
            Complex c1 = new Complex(1, 1);
            assertTrue(poles.get(i).equals(c.add(c1)));
            ++i;
            Complex c2 = new Complex(1, -1);
            assertTrue(testList.get(i).equals(c.conjugate()));
            assertTrue(poles.get(i).equals(c.conjugate().add(c2)));
          } else {
            assertTrue(poles.get(i).equals(c.add(1)));
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

      // System.out.println(testList);
      // System.out.println(poles);
      // System.out.println(poles2);

      for (int i = 0; i < poles.size(); ++i) {
        if (i < 2) {
          assertFalse(poles.get(i).equals(poles2.get(i)));
          assertTrue(poles2.get(i).equals(testList.get(i)));
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

      String[] resultString = RandomizedPanel.getInsetString(rCal);
      for (String resultPart : resultString) {
        sb.append(resultPart);
        sb.append('\n');
      }
      sb.append(RandomizedPanel.getTimeStampString(rCal));
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

      resultString = RandomizedPanel.getInsetString(rCal);
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
      RandData rd = cps.populateDataAndRun(calInFile, sensorOutFile, respFile,
          false, start, end, true);
      System.out.println(Arrays.toString(rd.getFitPoles()));
      System.out.println(Arrays.toString(rd.getFitZeros()));
    } catch (IOException | SeedFormatException | CodecException e) {
      e.printStackTrace();
      fail();
    }
  }

  @Test
  public void runExperiment_BCIP_HFCalibration() {
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

    assertTrue(rCal.hasEnoughData(ds));
    rCal.runExperimentOnData(ds);
    rCal.setLowFrequencyCalibration(false);
    List<Complex> fitPoles = rCal.getFitPoles();
    Complex[] expectedPoles = {
        new Complex(-306.7741224387797, 0),
        new Complex(-3.4804079210157486, 0),
        new Complex(-101.27715855875556, -387.9300826976112),
        new Complex(-101.27715855875556, 387.9300826976112)
    };
    ComplexFormat cf = new ComplexFormat(new DecimalFormat("#.#####"));
    StringBuilder report = new StringBuilder("The best-fit values from the solver: ");
    for (int i = 0; i < fitPoles.size(); i++) {
      report.append(cf.format(expectedPoles[i]) + ", ");
    }
    for (int i = 0; i < fitPoles.size(); i++) {
      assertEquals(report.toString(), expectedPoles[i].getReal(), fitPoles.get(i).getReal(), 1E-5);
      assertEquals(report.toString(), expectedPoles[i].getImaginary(), fitPoles.get(i).getImaginary(), 1E-5);
    }

    assertEquals(50.11489080838925, rCal.getFitResidual(), 1E-7);
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
    rCal.setLowFrequencyCalibration(false);
    List<Complex> fitPoles = rCal.getFitPoles();
    Complex[] expectedPoles = {
        new Complex(-0.012781625484629284, -0.012442058263140014),
        new Complex(-0.012781625484629284, 0.012442058263140014)
    };
    for (int i = 0; i < fitPoles.size(); i++) {
      assertEquals(expectedPoles[i].getReal(), fitPoles.get(i).getReal(), 1E-5);
      assertEquals(expectedPoles[i].getImaginary(), fitPoles.get(i).getImaginary(), 1E-5);
    }

    assertEquals(197.1889105489712, rCal.getFitResidual(), 1E-7);
    assertEquals(414.3706105547109, rCal.getInitResidual(), 1E-7);
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
