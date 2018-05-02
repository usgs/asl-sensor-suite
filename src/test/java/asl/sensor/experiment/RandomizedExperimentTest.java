package asl.sensor.experiment;

import static asl.sensor.test.TestUtils.RESP_LOCATION;
import static asl.sensor.test.TestUtils.getSeedFolder;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import asl.sensor.input.DataStoreUtils;
import asl.sensor.test.TestUtils;
import java.awt.Font;
import java.io.File;
import java.io.IOException;
import java.time.OffsetDateTime;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.linear.RealVector;
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
import asl.sensor.gui.RandomizedPanel;
import asl.sensor.input.DataStore;
import asl.sensor.input.InstrumentResponse;
import asl.sensor.utils.ReportingUtils;
import edu.iris.dmc.seedcodec.CodecException;
import edu.sc.seis.seisFile.mseed.SeedFormatException;

public class RandomizedExperimentTest {

  public static String folder = TestUtils.TEST_DATA_LOCATION + TestUtils.SUBPAGE;
  String testRespName = folder + "random-high-32+70i/RESP.XX.NS088..BHZ.STS1.360.2400";



  /*
  @Test
  public void TestRandomCalCurves() throws Exception {
    String fname = folder + "kiev-random-lowfrq/";
    String cal = "_BC0.512.seed";
    String out = "00_BH1.512.seed";
    InstrumentResponse ir = InstrumentResponse.loadEmbeddedResponse("STS25_Q330HR");
    DataBlock calB = TimeSeriesUtils.getFirstTimeSeries(fname + cal);
    DataBlock outB = TimeSeriesUtils.getFirstTimeSeries(fname + out);
    DataStore ds = new DataStore();
    ds.setBlock(0, calB);
    ds.setBlock(1, outB);
    ds.setResponse(1, ir);

    String startString = "2018-044T23:37:00.0";
    long st = TestUtils.timeStringToEpochMilli(startString);
    long ed = st + (8 * 60 * 60 * 1000);
    ds.trim(st, ed);

    RandomizedExperiment re = new RandomizedExperiment();
    re.setLowFreq(true);
    re.runExperimentOnData(ds);

    assertEquals(262144/2 + 1, re.getUntrimmedPSDLength());

    Complex ref = new Complex(-0.01243, -0.01176);
    Complex got = re.getFitPoles().get(0);

    String msg = "Expected " + ref + " and got " + got;
    assertTrue(msg, Complex.equals(ref, got, 5E-4));
  }
  */

  @Test
  public void ResponseCorrectConvertedToVectorHighFreq() throws Exception{
    String fname = folder + "resp-parse/TST5_response.txt";
    InstrumentResponse ir;
    ir = new InstrumentResponse(fname);
    List<Complex> poles = new ArrayList<>(ir.getPoles());
    // using an unnecessarily high nyquist rate here
    RealVector high = ir.polesToVector(false, 1E8);

    int complexIndex = 2; // start at second pole
    int vectorIndex = 0;

    while ( vectorIndex < high.getDimension() ) {
      // return current index
      double real = high.getEntry(vectorIndex++);
      double imag = high.getEntry(vectorIndex++);

      double poleImag = poles.get(complexIndex).getImaginary();

      assertEquals( real, poles.get(complexIndex).getReal(), 0.0 );
      assertEquals( imag, poleImag, 0.0 );

      if (poleImag != 0) {
        // complex conjugate case
        ++complexIndex;
        assertEquals( real, poles.get(complexIndex).getReal(), 0.0 );
        assertEquals( imag, -poles.get(complexIndex).getImaginary(), 0.0 );
      }
      ++complexIndex;
    }
  }

  @Test
  public void ResponseCorrectlyConvertedToVectorLowFreq() {
    String fname = folder + "resp-parse/TST5_response.txt";
    InstrumentResponse ir;
    try {

      ir = new InstrumentResponse(fname);
      List<Complex> poles = new ArrayList<>(ir.getPoles());
      // again, use a very high nyquist rate
      RealVector low = ir.polesToVector(true, 1E8);

      // only test lower two poless
      assertEquals( low.getEntry(0), poles.get(0).getReal(), 0.0 );
      assertEquals( low.getEntry(1), poles.get(0).getImaginary(), 0.0 );

      assertEquals( low.getEntry(0), poles.get(1).getReal(), 0.0 );
      assertEquals( low.getEntry(1), -poles.get(1).getImaginary(), 0.0 );

    } catch (IOException e) {
      fail();
      e.printStackTrace();
    }
  }

  @Test
  public void ResponseSetCorrectlyHighFreq() {
    String fname = folder + "resp-parse/TST5_response.txt";
    InstrumentResponse ir;

    try {
      ir = new InstrumentResponse(fname);

      List<Complex> poles = new ArrayList<>(ir.getPoles());
      List<Complex> replacements = new ArrayList<>();

      int start = 2;
      if ( poles.get(0).getImaginary() == 0 ) {
        start = 1;
      }

      for (int i = start; i < poles.size(); ++i) {
        if ( poles.get(i).getImaginary() == 0 ) {
          Complex c = poles.get(i);
          replacements.add(c.subtract(1));
          int next = i+1;
          while (next < poles.size() && poles.get(next).equals(c)) {
            ++next; // skip duplicates
          }
        } else {
          Complex c = poles.get(i);
          c = c.subtract( new Complex(1, 1) );
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
          assertTrue( poles.get(i).equals( testList.get(i) ) );
        } else {
          Complex c = replacements.get(offsetIdx);
          assertTrue( testList.get(i).equals(c) );
          if ( poles.get(i).getImaginary() != 0 ) {
            Complex c1 = new Complex(1, 1);
            assertTrue(poles.get(i).equals( c.add(c1) ));
            ++i;
            Complex c2 = new Complex(1, -1);
            assertTrue( testList.get(i).equals( c.conjugate() ) );
            assertTrue( poles.get(i).equals( c.conjugate().add(c2) ) );
          } else {
            assertTrue( poles.get(i).equals(c.add(1)) );
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

      Complex c = new Complex( newPoles[0], newPoles[1] );

      InstrumentResponse ir2 =
          ir.buildResponseFromFitVector(newPoles, true, 0);
      List<Complex> poles2 = ir2.getPoles();

      List<Complex> testList = new ArrayList<>(poles);
      testList.set(0, c);
      testList.set( 1, c.conjugate() );

      // System.out.println(testList);
      // System.out.println(poles);
      // System.out.println(poles2);

      for (int i = 0; i < poles.size(); ++i) {
        if (i < 2) {
          assertFalse( poles.get(i).equals( poles2.get(i) ) );
          assertTrue( poles2.get(i).equals( testList.get(i) ) );
        }
      }


    } catch (IOException e) {
      fail();
      e.printStackTrace();
    }

  }

  public DataStore setUpTest1() throws IOException {

    String respName = testRespName;
    String dataFolderName = folder + "random-high-32+70i/";
    String calName =  dataFolderName + "_EC0.512.seed";
    String sensOutName = dataFolderName + "00_EHZ.512.seed";

    DataStore ds = DataStoreUtils.createFromNames(respName, calName, sensOutName);
    OffsetDateTime cCal = TestUtils.getStartCalendar(ds);

    cCal = cCal.withMinute(36);
    cCal = cCal.withSecond(0);
    long start = cCal.toInstant().toEpochMilli();

    cCal = cCal.withMinute(41);
    // System.out.println( "end: " + sdf.format( cCal.getTime() ) );
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
          ExperimentFactory.createExperiment(ExperimentEnum.RANDM);

      rCal.setLowFreq(false);

      assertTrue( rCal.hasEnoughData(ds) );
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
        sb.append( resultPart );
        sb.append('\n');
      }
      sb.append( RandomizedPanel.getTimeStampString(rCal) );
      sb.append('\n');
      sb.append("Input files:\n");
      sb.append( ds.getBlock(0).getName() );
      sb.append(" (calibration)\n");
      sb.append( ds.getBlock(1).getName() );
      sb.append(" (sensor output)\n");
      sb.append("Response file used:\n");
      sb.append( ds.getResponse(1).getName() );
      sb.append("\n \n");

      String page1 = sb.toString();

      String[] addtlPages = ( RandomizedPanel.getAdditionalReportPages(rCal) );
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
          Math.abs( 100 * (bestResid - expectedResid) / bestResid );

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
        sb.append( resultPart );
        sb.append('\n');
      }

      for (int i = 0; i < jfcl.length; ++i) {

        jfcl[i] = ChartFactory.createXYLineChart(
            ExperimentEnum.RANDM.getName(),
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
      if ( !dir.exists() ) {
        dir.mkdir();
      }

      String testResult =
          testResultFolder + "Random-Calib-Test-1.pdf";
      pdf.save( new File(testResult) );
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
      System.out.println( Arrays.toString(rd.getFitPoles()) );
      System.out.println( Arrays.toString(rd.getFitZeros()) );
    } catch (IOException | SeedFormatException | CodecException e) {
      e.printStackTrace();
      fail();
    }
  }

  @Test
  public void runExperiment_BCIP_HFCalibration() {
    String respName = RESP_LOCATION + "RESP.CU.BCIP.00.BHZ_2017_268";
    String dataFolderName = getSeedFolder("CU", "BCIP", "2017", "268");
    String calName =  dataFolderName + "CB_BC0.512.seed";
    String sensOutName = dataFolderName + "00_EHZ.512.seed";

    DataStore ds = DataStoreUtils.createFromNames(respName, calName, sensOutName);
    OffsetDateTime cCal = TestUtils.getStartCalendar(ds);
    cCal = cCal.withHour(18).withMinute(49).withSecond(0).withNano(0);
    long start = cCal.toInstant().toEpochMilli();

    cCal = cCal.withHour(19).withMinute(4);
    long end = cCal.toInstant().toEpochMilli();

    ds.trim(start, end);

    RandomizedExperiment rCal = (RandomizedExperiment)
        ExperimentFactory.createExperiment(ExperimentEnum.RANDM);

    rCal.setLowFreq(false);

    assertTrue( rCal.hasEnoughData(ds) );
    rCal.runExperimentOnData(ds);
    List<Complex> fitPoles = rCal.getFitPoles();
    Complex[] expectedPoles = {
        //new Complex(-0.037, 0.037),
        //new Complex(-0.037, -0.037),
        new Complex(-374.8, 0),
        //new Complex(-520.3,0),
        //new Complex(-1053,-1005),
        //new Complex(-1053,1005),
        //new Complex(-13300.0,0),
        new Complex(34.36739870776367,0),
        new Complex(-101.17033498402674,388.4469528350759),
        new Complex(-101.17033498402674,-388.4469528350759),
        new Complex(-241.33932705943866,0)};
    for(int i = 0; i < fitPoles.size(); i++){
      assertEquals(expectedPoles[i].getReal(), fitPoles.get(i).getReal(), 1E-5);
      assertEquals(expectedPoles[i].getImaginary(), fitPoles.get(i).getImaginary(), 1E-5);
    }


  }

  @Test
  public void hasEnoughData_missingInputData() {
    String respName = RESP_LOCATION + "RESP.CU.BCIP.00.BHZ_2017_268";
    String dataFolderName = getSeedFolder("CU", "BCIP", "2017", "268");
    String sensOutName = dataFolderName + "00_EHZ.512.seed";

    DataStore ds = DataStoreUtils.createFromNames(respName, null, sensOutName);
    RandomizedExperiment rCal = (RandomizedExperiment)
        ExperimentFactory.createExperiment(ExperimentEnum.RANDM);
    assertFalse( rCal.hasEnoughData(ds) );
  }

  @Test
  public void hasEnoughData_missingOutputData() {
    String respName = RESP_LOCATION + "RESP.CU.BCIP.00.BHZ_2017_268";
    String dataFolderName = getSeedFolder("CU", "BCIP", "2017", "268");
    String calName =  dataFolderName + "CB_BC0.512.seed";

    DataStore ds = DataStoreUtils.createFromNames(respName, calName, null);
    RandomizedExperiment rCal = (RandomizedExperiment)
        ExperimentFactory.createExperiment(ExperimentEnum.RANDM);
    assertFalse( rCal.hasEnoughData(ds) );
  }

  @Test
  public void hasEnoughData_missingBothData() {
    String respName = RESP_LOCATION + "RESP.CU.BCIP.00.BHZ_2017_268";

    DataStore ds = DataStoreUtils.createFromNames(respName, null, null);
    RandomizedExperiment rCal = (RandomizedExperiment)
        ExperimentFactory.createExperiment(ExperimentEnum.RANDM);
    assertFalse( rCal.hasEnoughData(ds) );
  }

  @Test
  public void hasEnoughData_hasEnoughData() {
    String respName = RESP_LOCATION + "RESP.CU.BCIP.00.BHZ_2017_268";
    String dataFolderName = getSeedFolder("CU", "BCIP", "2017", "268");
    String calName =  dataFolderName + "CB_BC0.512.seed";
    String sensOutName = dataFolderName + "00_EHZ.512.seed";

    DataStore ds = DataStoreUtils.createFromNames(respName, calName, sensOutName);
    RandomizedExperiment rCal = (RandomizedExperiment)
        ExperimentFactory.createExperiment(ExperimentEnum.RANDM);
    assertTrue( rCal.hasEnoughData(ds) );
  }
}