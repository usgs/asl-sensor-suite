package asl.sensor.input;

import static junit.framework.TestCase.assertNotNull;
import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import asl.sensor.test.TestUtils;
import asl.sensor.utils.ReportingUtils;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.URISyntaxException;
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

public class InstrumentResponseTest {

  public static String folder = TestUtils.TEST_DATA_LOCATION + TestUtils.SUBPAGE;
  public static final DateTimeFormatter DATE_TIME_FORMAT = InstrumentResponse.RESP_DT_FORMAT;

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
      pls.add(new Complex(-2.221000e-01, 2.221000e-01));
      pls.add(new Complex(-2.221000e-01, -2.221000e-01));
      pls.add(new Complex(-3.918000e+01, 4.912000e+01));
      pls.add(new Complex(-3.918000e+01, -4.912000e+01));
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

    assertTrue(initPoles.size() == endPoles.size());

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
    assertEquals(null, actual);
  }

  @Test
  public void parseTermAsDate_starttime_no_seconds() {
    Instant actual = InstrumentResponse.parseTermAsDate("B052F23     End date:  2007,337,10:15");
    Instant expected = LocalDateTime.parse("2007-12-03T10:15:00").toInstant(ZoneOffset.UTC);
    assertEquals(expected, actual);
  }

  @Test
  public void parseTermAsComplex_checkThatComplexArrayWasOrderedCorrectly() {
    fail();
  }

  @Test
  public void parseTermAsComplex_checkThatComplexParsedCorrectly() {
    fail();
  }

  @Test
  public void parseTransferType_A_laplacianTransfer(){
    assertEquals(TransferFunction.LAPLACIAN, InstrumentResponse.parseTransferType("A"));
  }

  @Test
  public void parseTransferType_B_linearTransfer(){
    assertEquals(TransferFunction.LINEAR, InstrumentResponse.parseTransferType("B"));
  }

  /**
   * Test Composite Transfer from SEED specifications
   */
  @Test
  public void parseTransferType_C_defaultsToLaplacian(){
    assertEquals(TransferFunction.LAPLACIAN, InstrumentResponse.parseTransferType("C"));
  }

  /**
   * Test Digital (Z) Transfer from SEED specifications
   */
  @Test
  public void parseTransferType_D_defaultsToLaplacian(){
    assertEquals(TransferFunction.LAPLACIAN, InstrumentResponse.parseTransferType("D"));
  }

  @Test
  public void parseUnitType_velocity_checkCaseSensitivity() throws Exception{
    assertEquals(Unit.VELOCITY, InstrumentResponse.parseUnitType("M/s"));
  }

  @Test
  public void parseUnitType_velocity() throws Exception{
    assertEquals(Unit.VELOCITY, InstrumentResponse.parseUnitType("m/s"));
  }

  @Test
  public void parseUnitType_acceleration() throws Exception{
    assertEquals(Unit.ACCELERATION, InstrumentResponse.parseUnitType("m/s**2"));
  }

  @Test (expected = IOException.class)
  public void parseUnitType_defaultThrowsException() throws Exception{
    InstrumentResponse.parseUnitType("Amperage");
  }

  @Test
  public void skipToSelectedEpoch_methodReturnsEvenIfEpochMissing() throws Exception {
    BufferedReader reader = new BufferedReader(new InputStreamReader(InstrumentResponseTest.class.getResourceAsStream("/seismic-test-data/RESPs/RESP_with_empty_lines")));
    InstrumentResponse.skipToSelectedEpoch(reader, Instant.MAX);
    //Should be empty since it was fully read
    assertFalse(reader.ready());
  }

  @Test
  public void skipToSelectedEpoch_doesItStopOnCorrectEpochAfterSkippingOverEpoch_doesItSkipEmptyLines() throws Exception {
    LocalDate date = LocalDate.parse("2012-01-01");
    Instant epoch = date.atStartOfDay(ZoneOffset.UTC).toInstant();

    BufferedReader reader = new BufferedReader(new InputStreamReader(InstrumentResponseTest.class.getResourceAsStream("/seismic-test-data/RESPs/RESP_with_empty_lines")));
    InstrumentResponse.skipToSelectedEpoch(reader, epoch);
    assertTrue(reader.ready());

    //Verify that we are in the correct place
    assertEquals("B052F23     End date:    No Ending Time", reader.readLine());
  }

  @Test
  public void skipToSelectedEpoch_doesItStopOnFirstEpoch() throws Exception {
    LocalDate date = LocalDate.parse("2006-01-01");
    Instant epoch = date.atStartOfDay(ZoneOffset.UTC).toInstant();

    BufferedReader reader = new BufferedReader(new InputStreamReader(InstrumentResponseTest.class.getResourceAsStream("/seismic-test-data/RESPs/RESP_with_empty_lines")));
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
  public void parseInstrumentList() throws IOException{
    Set<String> respFiles = InstrumentResponse.parseInstrumentList();
    assertEquals(12, respFiles.size());
    for (String respFile: respFiles) {
      InstrumentResponse response = InstrumentResponse.loadEmbeddedResponse(respFile);
      assertNotNull(response);
    }
  }

  /**
   * This test tests that each field was parsed correctly and the correct epoch was parsed.
   */
  @Test
  public void parserDriver_notNullEpochParsesCorrectEpoch_firstEpoch() throws Exception {
    LocalDate date = LocalDate.parse("2006-01-01");
    Instant epochStart = date.atStartOfDay(ZoneOffset.UTC).toInstant();
    URL file = InstrumentResponseTest.class.getResource("/seismic-test-data/RESPs/RESP_with_empty_lines");

    InstrumentResponse response = new InstrumentResponse(Paths.get(file.toURI()).toString(), epochStart);


    assertEquals(epochStart, response.getEpochStart());

    LocalDate endDate = LocalDate.parse("2012-01-01");
    Instant epochEnd = endDate.atStartOfDay(ZoneOffset.UTC).toInstant();
    assertEquals(epochEnd, response.getEpochEnd());

    assertEquals(TransferFunction.LAPLACIAN, response.getTransferFunction());
    assertEquals(Unit.VELOCITY, response.getUnits());
    assertEquals(3.948580E3, response.getNormalization(), 1E-6);
    assertEquals(2E-2, response.getNormalizationFrequency(), 1E-6);
    assertEquals(3, response.getNumStages());

    double[] gain = {4.026530e+09, 2.400000e+03, 1.677721e+06, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    assertArrayEquals(gain, response.getGain(), 1E-6);

    List<Complex> zeros = new ArrayList<>();
    zeros.add(new Complex(0.000000e+00, 0));
    zeros.add(new Complex(0.000000e+00, 0));
    assertEquals(zeros, response.getZeros());

    List<Complex> poles = new ArrayList<>();
    poles.add(new Complex(-11.234000e-02, 15.234000e-02));
    poles.add(new Complex(-11.234000e-02, -12.234000e-02));
    poles.add(new Complex(-31.918000e+01, 5.912000e+01));
    poles.add(new Complex(-31.918000e+01, -5.912000e+01));
    assertEquals(poles, response.getPoles());

  }

  /**
   * This test tests that each field was parsed correctly and the correct epoch was parsed.
   */
  @Test
  public void parserDriver_notNullEpochParsesCorrectEpoch_lastEpoch() throws Exception{
    LocalDate date = LocalDate.parse("2012-01-01");
    Instant epochStart = date.atStartOfDay(ZoneOffset.UTC).toInstant();
    URL file = InstrumentResponseTest.class.getResource("/seismic-test-data/RESPs/RESP_with_empty_lines");

    InstrumentResponse response = new InstrumentResponse(Paths.get(file.toURI()).toString(), epochStart);


    assertEquals(epochStart, response.getEpochStart());

    assertNull(response.getEpochEnd());

    assertEquals(TransferFunction.LINEAR, response.getTransferFunction());
    assertEquals(Unit.ACCELERATION, response.getUnits());
    assertEquals(2.948580E3, response.getNormalization(), 1E-6);
    assertEquals(1E-2, response.getNormalizationFrequency(), 1E-6);
    assertEquals(6, response.getNumStages());

    double[] gain = {1.0, 1.400000e+02, 1.677721e+06, 1.0, 1.0, 1.026530e+09, 1.0, 1.0, 1.0, 1.0};
    assertArrayEquals(gain, response.getGain(), 1E-6);

    List<Complex> zeros = new ArrayList<>();
    zeros.add(new Complex(3.000000e+00, 0));
    assertEquals(zeros, response.getZeros());

    List<Complex> poles = new ArrayList<>();
    poles.add(new Complex(-2.234000e-02, 1.234000e-02));
    poles.add(new Complex(-2.234000e-02, -1.234000e-02));
    poles.add(new Complex(-1.918000e+01, 4.912000e+01));
    assertEquals(poles, response.getPoles());
  }

  /**
   * This test tests that each field was parsed correctly and the correct epoch was parsed.
   */
  @Test
  public void parserDriver_nullEpochParsesLastEpoch() throws Exception{
    LocalDate date = LocalDate.parse("2012-01-01");
    Instant epochStart = date.atStartOfDay(ZoneOffset.UTC).toInstant();
    URL file = InstrumentResponseTest.class.getResource("/seismic-test-data/RESPs/RESP_with_empty_lines");

    InstrumentResponse response = new InstrumentResponse(Paths.get(file.toURI()).toString());


    assertEquals(epochStart, response.getEpochStart());

    assertNull(response.getEpochEnd());

    assertEquals(TransferFunction.LINEAR, response.getTransferFunction());
    assertEquals(Unit.ACCELERATION, response.getUnits());
    assertEquals(2.948580E3, response.getNormalization(), 1E-6);
    assertEquals(1E-2, response.getNormalizationFrequency(), 1E-6);
    assertEquals(6, response.getNumStages());

    double[] gain = {1.0, 1.400000e+02, 1.677721e+06, 1.0, 1.0, 1.026530e+09, 1.0, 1.0, 1.0, 1.0};
    assertArrayEquals(gain, response.getGain(), 1E-6);

    List<Complex> zeros = new ArrayList<>();
    zeros.add(new Complex(3.000000e+00, 0));
    assertEquals(zeros, response.getZeros());

    List<Complex> poles = new ArrayList<>();
    poles.add(new Complex(-2.234000e-02, 1.234000e-02));
    poles.add(new Complex(-2.234000e-02, -1.234000e-02));
    poles.add(new Complex(-1.918000e+01, 4.912000e+01));
    assertEquals(poles, response.getPoles());
  }

}
