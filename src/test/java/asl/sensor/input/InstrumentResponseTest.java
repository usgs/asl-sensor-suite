package asl.sensor.input;

import static junit.framework.TestCase.assertNotNull;
import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.URL;
import java.nio.file.Paths;
import java.time.Instant;
import java.time.LocalDate;
import java.time.LocalDateTime;
import java.time.ZoneOffset;
import java.time.format.DateTimeFormatter;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;
import org.apache.pdfbox.pdmodel.PDDocument;
import org.junit.Test;
import asl.sensor.test.TestUtils;
import asl.sensor.utils.ReportingUtils;

public class InstrumentResponseTest {

  private static final String folder = TestUtils.TEST_DATA_LOCATION + TestUtils.SUBPAGE;
  private static final DateTimeFormatter DATE_TIME_FORMAT = InstrumentResponse.RESP_DT_FORMAT;

  @Test
  public void testGetClosestRespEpoch_beforeFirstEpoch() throws IOException {
    String filepath = TestUtils.RESP_LOCATION + "RESP.CU.BCIP.00.BHZ_2017_268";
    long start = 0;
    long end = 100;

    Instant expectedClosest =
        LocalDateTime.parse("2010,041,18:35:00", DATE_TIME_FORMAT).toInstant(ZoneOffset.UTC);
    Instant returnedClosest = InstrumentResponse.getRespFileClosestEpoch(filepath, start, end);
    assertEquals(expectedClosest, returnedClosest);
  }

  @Test
  public void testGetClosestRespEpoch_afterLastEpoch() throws IOException {
    String filepath = TestUtils.RESP_LOCATION + "RESP.CU.BCIP.00.BHZ_2017_268";
    long start =
        LocalDateTime.parse("2800,041,18:35:00", DATE_TIME_FORMAT).toInstant(ZoneOffset.UTC).toEpochMilli();
    long end =
        LocalDateTime.parse("2820,041,18:35:00", DATE_TIME_FORMAT).toInstant(ZoneOffset.UTC).toEpochMilli();

    Instant expectedClosest =
        LocalDateTime.parse("2015,055,00:00:00", DATE_TIME_FORMAT).toInstant(ZoneOffset.UTC);
    Instant returnedClosest = InstrumentResponse.getRespFileClosestEpoch(filepath, start, end);
    assertEquals(expectedClosest, returnedClosest);
  }

  @Test
  public void testGetClosestRespEpoch_betweenEpochBoundaries() throws IOException {
    String filepath = TestUtils.RESP_LOCATION + "RESP.CU.BCIP.00.BHZ_2017_268";
    long start =
        LocalDateTime.parse("2015,054,00:00:00", DATE_TIME_FORMAT).toInstant(ZoneOffset.UTC).toEpochMilli();
    long end =
        LocalDateTime.parse("2015,056,00:00:00", DATE_TIME_FORMAT).toInstant(ZoneOffset.UTC).toEpochMilli();

    // though data starts in an earlier epoch, it is closer to the start of this one
    Instant expectedClosest =
        LocalDateTime.parse("2015,055,00:00:00", DATE_TIME_FORMAT).toInstant(ZoneOffset.UTC);
    Instant returnedClosest = InstrumentResponse.getRespFileClosestEpoch(filepath, start, end);
    assertEquals(expectedClosest, returnedClosest);
  }

  @Test
  public void testGetClosestRespEpoch_betweenEpochBoundariesFirstCloser() throws IOException {
    String filepath = TestUtils.RESP_LOCATION + "RESP.CU.BCIP.00.BHZ_2017_268";
    long start =
        LocalDateTime.parse("2010,050,00:00:00", DATE_TIME_FORMAT).toInstant(ZoneOffset.UTC).toEpochMilli();
    long end =
        LocalDateTime.parse("2015,056,00:00:00", DATE_TIME_FORMAT).toInstant(ZoneOffset.UTC).toEpochMilli();

    Instant expectedClosest =
        LocalDateTime.parse("2010,041,18:35:00", DATE_TIME_FORMAT).toInstant(ZoneOffset.UTC);
    Instant returnedClosest = InstrumentResponse.getRespFileClosestEpoch(filepath, start, end);
    assertEquals(expectedClosest, returnedClosest);
  }

  @Test
  public void testGetClosestRespEpoch_insideEpochBoundaries() throws IOException {
    String filepath = TestUtils.RESP_LOCATION + "RESP.CU.BCIP.00.BHZ_2017_268";
    long start =
        LocalDateTime.parse("2012,050,00:00:00", DATE_TIME_FORMAT).toInstant(ZoneOffset.UTC).toEpochMilli();
    long end =
        LocalDateTime.parse("2012,052,00:00:00", DATE_TIME_FORMAT).toInstant(ZoneOffset.UTC).toEpochMilli();

    Instant expectedClosest =
        LocalDateTime.parse("2010,041,18:35:00", DATE_TIME_FORMAT).toInstant(ZoneOffset.UTC);
    Instant returnedClosest = InstrumentResponse.getRespFileClosestEpoch(filepath, start, end);
    assertEquals(expectedClosest, returnedClosest);
  }

  @Test
  public void testFileParse() {
    String filename = folder + "resp-parse/RESP.XX.NS087..BHZ.STS1.20.2400";

    try {
      InstrumentResponse ir = new InstrumentResponse(filename);

      assertEquals(TransferFunction.LAPLACIAN, ir.getTransferFunction());

      double nml = Double.parseDouble("3.948580E+03");
      assertEquals(nml, ir.getNormalization(), 0.0001);

      double nmf = Double.parseDouble("3.000000E-01");
      assertEquals(nmf, ir.getNormalizationFrequency(), 0.0001);

      double[] gn = {2.400000e+03, 2.400000e+03, 1.000000e+00};
      int maxStage = ir.getNumStages();
      assertEquals(maxStage, gn.length);
      for (int i = 0; i < maxStage; ++i) {
        assertEquals(gn[i], ir.getGain()[i], 1.);
      }
      //assertTrue( gnL.equals(ir.getGain() ) );

      assertEquals(Unit.VELOCITY, ir.getUnits());

      List<Complex> zrs = new ArrayList<>();
      zrs.add(new Complex(0.000000e+00, 0.000000e+00));
      zrs.add(new Complex(0.000000e+00, 0.000000e+00));
      assertEquals(zrs, ir.getZeros());

      List<Complex> pls = new ArrayList<>();
      pls.add(new Complex(-2.221000e-01, -2.221000e-01));
      pls.add(new Complex(-2.221000e-01, 2.221000e-01));
      pls.add(new Complex(-3.918000e+01, -4.912000e+01));
      pls.add(new Complex(-3.918000e+01, 4.912000e+01));
      assertEquals(pls, ir.getPoles());

    } catch (IOException e) {
      fail("Unexpected error trying to read response file");
    }

  }

  @Test
  public void testMultiEpoch() {

    String filename = folder + "resp-parse/multiepoch.txt";
    try {
      List<Pair<Instant, Instant>> eps = InstrumentResponse.getRespFileEpochs(filename);
      assertTrue(eps.size() > 1);
    } catch (IOException e) {
      fail();
      e.printStackTrace();
    }

  }

  @Test
  public void listsAllEpochs() {

    // epoch 1 ends same time epoch 2 begins
    String[] times = {"2016,193,00:00:00", "2016,196,00:00:00", "2016,224,00:00:00"};
    Instant[] insts = new Instant[times.length];
    for (int i = 0; i < times.length; ++i) {
      insts[i] = LocalDateTime.parse(times[i], DATE_TIME_FORMAT).toInstant(ZoneOffset.UTC);
    }
    List<Pair<Instant, Instant>> compareTo = new ArrayList<>();
    compareTo.add(new Pair<>(insts[0], insts[1]));
    compareTo.add(new Pair<>(insts[1], insts[2]));

    String filename = folder + "resp-parse/multiepoch.txt";
    try {
      List<Pair<Instant, Instant>> eps = InstrumentResponse.getRespFileEpochs(filename);
      for (int i = 0; i < eps.size(); ++i) {
        Pair<Instant, Instant> inst = eps.get(i);
        Pair<Instant, Instant> base = compareTo.get(i);
        assertEquals(inst.getFirst(), base.getFirst());
        assertEquals(inst.getSecond(), base.getSecond());
      }
    } catch (IOException e) {
      fail();
      e.printStackTrace();
    }
  }

  @Test
  public void testStringOutput() {

    String filename = folder + "resp-parse/RESP.XX.NS087..BHZ.STS1.20.2400";

    try {
      InstrumentResponse ir = new InstrumentResponse(filename);

      System.out.println(ir);

      PDDocument pdf = new PDDocument();

      ReportingUtils.textToPDFPage(ir.toString(), pdf);

      String currentDir = System.getProperty("user.dir");
      String testResultFolder = currentDir + "/testResultImages/";
      File dir = new File(testResultFolder);
      if (!dir.exists()) {
        dir.mkdir();
      }

      String testResult = testResultFolder + "response-report.pdf";
      pdf.save(new File(testResult));
      pdf.close();

    } catch (IOException e) {
      fail("Unexpected error trying to read response file");
    }
  }

  @Test
  public void vectorCreationRespectsDuplicatePoles() {

    InstrumentResponse ir;
    try {
      ir = new InstrumentResponse(TestUtils.RESP_LOCATION + "STS-5A_Q330HR_BH_40");
    } catch (IOException e) {
      fail();
      e.printStackTrace();
      return;
    }

    List<Complex> initPoles = new ArrayList<>(ir.getPoles());
    RealVector rv = ir.polesToVector(false, 100.);
    Complex c = new Complex(-20, 0);
    // poles at indices 2 and 3 are duplicated, have zero imaginary component
    // set them to a new value to test array resetting with diff. values
    initPoles.set(2, c);
    initPoles.set(3, c);
    // build new vector, is it the same?
    List<Complex> endPoles =
        ir.buildResponseFromFitVector(rv.toArray(), false, 0).getPoles();

    assertEquals(initPoles.size(), endPoles.size());

  }

  @Test
  public void parseTermAsDate_date_standard_input_field22() {
    Instant actual = InstrumentResponse
        .parseTermAsDate("B052F22     Start date:  2001,001,00:00:00");
    Instant expected = LocalDate.parse("2001-01-01").atStartOfDay().toInstant(ZoneOffset.UTC);
    assertEquals(expected, actual);
  }

  @Test
  public void parseTermAsDate_datetime_standard_input_field23() {
    Instant actual = InstrumentResponse.parseTermAsDate("B052F23     End date:  2007,337,10:15:30");
    Instant expected = LocalDateTime.parse("2007-12-03T10:15:30").toInstant(ZoneOffset.UTC);
    assertEquals(expected, actual);
  }

  @Test
  public void parseTermAsDate_noEndingTime_input_field23() {
    Instant actual = InstrumentResponse.parseTermAsDate("B052F23     End date:    No Ending Time");
    assertNull(actual);
  }

  @Test
  public void parseTermAsDate_starttime_no_seconds() {
    Instant actual = InstrumentResponse.parseTermAsDate("B052F23     End date:  2007,337,10:15");
    Instant expected = LocalDateTime.parse("2007-12-03T10:15:00").toInstant(ZoneOffset.UTC);
    assertEquals(expected, actual);
  }

  @Test
  public void parseTermAsDate_starttime_no_minutes() {
    Instant actual = InstrumentResponse.parseTermAsDate("B052F23     End date:  2007,337,01");
    Instant expected = LocalDateTime.parse("2007-12-03T01:00:00").toInstant(ZoneOffset.UTC);
    assertEquals(expected, actual);
  }

  @Test
  public void parseTermAsDate_starttime_no_hours() {
    Instant actual = InstrumentResponse.parseTermAsDate("B052F23     End date:  2007,337");
    Instant expected = LocalDateTime.parse("2007-12-03T00:00:00").toInstant(ZoneOffset.UTC);
    assertEquals(expected, actual);
  }


  @Test
  public void parseTermAsComplex_checkThatComplexArrayWasOrderedCorrectly() {
    Complex[] arr = new Complex[5];
    String[] lines = {"B053F15-18    0 -5.943130E+01  0.000000E+00  0.000000E+00  0.000000E+00",
        "B053F15-18    1 -2.271210E+01  2.710650E+01  0.000000E+00  0.000000E+00",
        "B053F15-18    2 -2.271210E+01 -2.710650E+01  0.000000E+00  0.000000E+00",
        "B053F15-18    3 -4.800400E-03  0.000000E+00  0.000000E+00  0.000000E+00",
        "B053F15-18    4 -7.394060E-02  0.000000E+00  7.594130E-04  0.000000E+00"};
    Complex[] testAgainst = {
        new Complex(-5.943130E1),
        new Complex(-2.271210E1, 2.710650E1),
        new Complex(-2.271210E1, -2.710650E1),
        new Complex(-4.800400E-3),
        new Complex(-7.394060E-2)
    };
    for (String line : lines) {
      InstrumentResponse.parseTermAsComplex(line, arr);
    }
    for (int i = 0; i < testAgainst.length; ++i) {
      assertTrue(Complex.equals(testAgainst[i], arr[i], 1E-10));
    }
  }

  @Test
  public void parseTermAsComplex_checkThatComplexParsedCorrectly() {
    Complex[] arr = new Complex[5];
    String line = "B053F15-18    2 -2.271210E+01 -2.710650E+01  0.000000E+00  0.000000E+00";
    InstrumentResponse.parseTermAsComplex(line, arr);
    Complex c = new Complex(-2.271210E1, -2.710650E1);
    assertTrue(Complex.equals(c, arr[2], 1E-10));
  }

  @Test
  public void parseTransferType_A_laplacianTransfer() {
    assertEquals(TransferFunction.LAPLACIAN, InstrumentResponse.parseTransferType("A"));
  }

  @Test
  public void parseTransferType_B_linearTransfer() {
    assertEquals(TransferFunction.LINEAR, InstrumentResponse.parseTransferType("B"));
  }

  /**
   * Test Composite Transfer from SEED specifications
   */
  @Test
  public void parseTransferType_C_defaultsToLaplacian() {
    assertEquals(TransferFunction.LAPLACIAN, InstrumentResponse.parseTransferType("C"));
  }

  /**
   * Test Digital (Z) Transfer from SEED specifications
   */
  @Test
  public void parseTransferType_D_defaultsToLaplacian() {
    assertEquals(TransferFunction.LAPLACIAN, InstrumentResponse.parseTransferType("D"));
  }

  @Test
  public void parseUnitType_velocity_checkCaseSensitivity() throws Exception {
    assertEquals(Unit.VELOCITY, InstrumentResponse.parseUnitType("M/s"));
  }

  @Test
  public void parseUnitType_velocity() throws Exception {
    assertEquals(Unit.VELOCITY, InstrumentResponse.parseUnitType("m/s"));
  }

  @Test
  public void parseUnitType_acceleration() throws Exception {
    assertEquals(Unit.ACCELERATION, InstrumentResponse.parseUnitType("m/s**2"));
  }

  @Test(expected = IOException.class)
  public void parseUnitType_defaultThrowsException() throws Exception {
    InstrumentResponse.parseUnitType("Amperage");
  }

  @Test
  public void skipToSelectedEpoch_methodReturnsEvenIfEpochMissing() throws Exception {
    BufferedReader reader = new BufferedReader(new InputStreamReader(InstrumentResponseTest.class
        .getResourceAsStream("/seismic-test-data/RESPs/RESP_with_empty_lines")));
    InstrumentResponse.skipToSelectedEpoch(reader, Instant.MAX);
    //Should be empty since it was fully read
    assertFalse(reader.ready());
  }

  @Test
  public void skipToSelectedEpoch_doesItStopOnCorrectEpochAfterSkippingOverEpoch_doesItSkipEmptyLines()
      throws Exception {
    LocalDate date = LocalDate.parse("2012-01-01");
    Instant epoch = date.atStartOfDay(ZoneOffset.UTC).toInstant();

    BufferedReader reader = new BufferedReader(new InputStreamReader(InstrumentResponseTest.class
        .getResourceAsStream("/seismic-test-data/RESPs/RESP_with_empty_lines")));
    InstrumentResponse.skipToSelectedEpoch(reader, epoch);
    assertTrue(reader.ready());

    //Verify that we are in the correct place
    assertEquals("B052F23     End date:    No Ending Time", reader.readLine());
  }

  @Test
  public void skipToSelectedEpoch_doesItStopOnFirstEpoch() throws Exception {
    LocalDate date = LocalDate.parse("2006-01-01");
    Instant epoch = date.atStartOfDay(ZoneOffset.UTC).toInstant();

    BufferedReader reader = new BufferedReader(new InputStreamReader(InstrumentResponseTest.class
        .getResourceAsStream("/seismic-test-data/RESPs/RESP_with_empty_lines")));
    InstrumentResponse.skipToSelectedEpoch(reader, epoch);
    assertTrue(reader.ready());

    //Verify that we are in the correct place
    assertEquals("B052F23     End date:    2012,001,00:00:00.0000", reader.readLine());
  }

  /**
   * Verify that responses.txt parses without error.
   * Verify that each file in responses.txt can be loaded.
   */
  @Test
  public void parseInstrumentList() throws IOException {
    Set<String> respFiles = InstrumentResponse.parseInstrumentList();
    assertEquals(12, respFiles.size());
    for (String respFile : respFiles) {
      InstrumentResponse response = InstrumentResponse.loadEmbeddedResponse(respFile);
      assertNotNull(response);
    }
  }

  @Test
  public void poleZerosMagnitudeIncreasing() throws IOException {
    InstrumentResponse resp =
        new InstrumentResponse(TestUtils.RESP_LOCATION + "STS6_Q330HR");
    Complex[] poles = resp.getPoles().toArray(new Complex[]{});
    Complex[] expected = {
        new Complex(-1.23000e-02, -1.23000e-02),
        new Complex(-1.23000e-02, +1.23000e-02),
        new Complex(-1.36833e+02),
        new Complex(-1.36833e+02),
        new Complex(-1.21500e+02, -6.47000e+02),
        new Complex(-1.21500e+02, +6.47000e+02),
        new Complex(-7.64103e+02),
        new Complex(-7.64103e+02),
        new Complex(-7.64103e+02),
        new Complex(-7.64103e+02),
        new Complex(-7.64103e+02),
        new Complex(-7.64103e+02),
        new Complex(-1.04409e+10)
    };
    for (int i = 0; i < poles.length; ++i) {
      assertTrue(Complex.equals(expected[i], poles[i], 1E-10));
    }
  }

  /**
   * This test tests that each field was parsed correctly and the correct epoch was parsed.
   */
  @Test
  public void parserDriver_notNullEpochParsesCorrectEpoch_firstEpoch() throws Exception {
    LocalDate date = LocalDate.parse("2006-01-01");
    Instant epochStart = date.atStartOfDay(ZoneOffset.UTC).toInstant();
    URL file = InstrumentResponseTest.class
        .getResource("/seismic-test-data/RESPs/RESP_with_empty_lines");

    InstrumentResponse response = new InstrumentResponse(Paths.get(file.toURI()).toString(),
        epochStart);

    assertEquals(epochStart, response.getEpochStart());

    LocalDate endDate = LocalDate.parse("2012-01-01");
    Instant epochEnd = endDate.atStartOfDay(ZoneOffset.UTC).toInstant();
    assertEquals(epochEnd, response.getEpochEnd());

    assertEquals(TransferFunction.LAPLACIAN, response.getTransferFunction());
    assertEquals(Unit.VELOCITY, response.getUnits());
    assertEquals(3.948580E3, response.getNormalization(), 1E-6);
    assertEquals(2E-2, response.getNormalizationFrequency(), 1E-6);
    assertEquals(3, response.getNumStages());

    double[] gain = {4.026530e+09, 2.400000e+03, 1.677721e+06};
    assertArrayEquals(gain, response.getGain(), 1E-6);

    List<Complex> zeros = new ArrayList<>();
    zeros.add(new Complex(0.000000e+00, 0));
    zeros.add(new Complex(0.000000e+00, 0));
    assertEquals(zeros, response.getZeros());

    List<Complex> poles = new ArrayList<>();
    poles.add(new Complex(-11.234000e-02, -12.234000e-02));
    poles.add(new Complex(-11.234000e-02, 15.234000e-02));
    poles.add(new Complex(-31.918000e+01, -5.912000e+01));
    poles.add(new Complex(-31.918000e+01, 5.912000e+01));
    assertEquals(poles, response.getPoles());

  }

  /**
   * This test tests that each field was parsed correctly and the correct epoch was parsed.
   */
  @Test
  public void parserDriver_notNullEpochParsesCorrectEpoch_lastEpoch() throws Exception {
    LocalDate date = LocalDate.parse("2012-01-01");
    Instant epochStart = date.atStartOfDay(ZoneOffset.UTC).toInstant();
    URL file = InstrumentResponseTest.class
        .getResource("/seismic-test-data/RESPs/RESP_with_empty_lines");

    InstrumentResponse response = new InstrumentResponse(Paths.get(file.toURI()).toString(),
        epochStart);

    assertEquals(epochStart, response.getEpochStart());

    assertNull(response.getEpochEnd());

    assertEquals(TransferFunction.LINEAR, response.getTransferFunction());
    assertEquals(Unit.ACCELERATION, response.getUnits());
    assertEquals(2.948580E3, response.getNormalization(), 1E-6);
    assertEquals(1E-2, response.getNormalizationFrequency(), 1E-6);
    assertEquals(6, response.getNumStages());

    double[] gain = {1.0, 1.400000e+02, 1.677721e+06, 1.0, 1.0, 1.026530e+09};
    assertArrayEquals(gain, response.getGain(), 1E-6);

    List<Complex> zeros = new ArrayList<>();
    zeros.add(new Complex(3.000000e+00, 0));
    assertEquals(zeros, response.getZeros());

    List<Complex> poles = new ArrayList<>();
    poles.add(new Complex(-2.234000e-02, -1.234000e-02));
    poles.add(new Complex(-2.234000e-02, 1.234000e-02));
    poles.add(new Complex(-1.918000e+01, 4.912000e+01));
    assertEquals(poles, response.getPoles());
  }

  /**
   * This test tests that each field was parsed correctly and the correct epoch was parsed.
   */
  @Test
  public void parserDriver_nullEpochParsesLastEpoch() throws Exception {
    LocalDate date = LocalDate.parse("2012-01-01");
    Instant epochStart = date.atStartOfDay(ZoneOffset.UTC).toInstant();
    URL file = InstrumentResponseTest.class
        .getResource("/seismic-test-data/RESPs/RESP_with_empty_lines");

    InstrumentResponse response = new InstrumentResponse(Paths.get(file.toURI()).toString());

    assertEquals(epochStart, response.getEpochStart());

    assertNull(response.getEpochEnd());

    assertEquals(TransferFunction.LINEAR, response.getTransferFunction());
    assertEquals(Unit.ACCELERATION, response.getUnits());
    assertEquals(2.948580E3, response.getNormalization(), 1E-6);
    assertEquals(1E-2, response.getNormalizationFrequency(), 1E-6);
    assertEquals(6, response.getNumStages());

    double[] gain = {1.0, 1.400000e+02, 1.677721e+06, 1.0, 1.0, 1.026530e+09};
    assertArrayEquals(gain, response.getGain(), 1E-6);

    List<Complex> zeros = new ArrayList<>();
    zeros.add(new Complex(3.000000e+00, 0));
    assertEquals(zeros, response.getZeros());

    List<Complex> poles = new ArrayList<>();
    poles.add(new Complex(-2.234000e-02, -1.234000e-02));
    poles.add(new Complex(-2.234000e-02, 1.234000e-02));
    poles.add(new Complex(-1.918000e+01, 4.912000e+01));
    assertEquals(poles, response.getPoles());
  }

  /**
   * This test tests that a resp with 10 gain stages (11 including stage 0) is parsed correctly.
   */
  @Test
  public void parserDriver_manyGainStages() throws Exception {
    LocalDate date = LocalDate.parse("2012-01-01");
    URL file = InstrumentResponseTest.class
        .getResource("/seismic-test-data/RESPs/BN.EKG.HHZ.resp");

    InstrumentResponse response = new InstrumentResponse(Paths.get(file.toURI()).toString());

    assertEquals(11, response.getNumStages());

    double[] gain = {1.246063e+09, 3.987400e+03, 3.125000e+05, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    assertArrayEquals(gain, response.getGain(), 1E-6);
  }

}
