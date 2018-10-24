package asl.sensor.input;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.text.NumberFormat;
import java.time.Instant;
import java.time.LocalDate;
import java.time.LocalDateTime;
import java.time.ZoneOffset;
import java.time.format.DateTimeFormatter;
import java.time.format.DateTimeParseException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.complex.ComplexFormat;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;
import asl.sensor.gui.InputPanel;
import asl.sensor.utils.NumericUtils;

/**
 * This class is used to read in and store data from instrument response files
 * such as those found on the IRIS RESP database. This includes functions used
 * to parse in relevant data from those files to get values such as poles and
 * zeros and the gain stages of the detector/sensor setup.
 * See also: http://ds.iris.edu/ds/nodes/dmc/data/formats/resp/
 *
 * @author akearns
 */
public class InstrumentResponse {

  /**
   * Maximum proportion of nyquist rate of a signal response to fit (as in RandomizedExperiment)
   */
  public static final double PEAK_MULTIPLIER = 0.9;
  /**
   * Julian date parser for expected format of response files
   */
  public static final DateTimeFormatter RESP_DT_FORMAT =
      DateTimeFormatter.ofPattern("uuuu,DDD,HH:mm:ss").withZone(ZoneOffset.UTC);
  /**
   * Initial value for size of array of gain stages; actual number of gain stages can be more or
   * less than the value given here
   */
  private static final int INIT_GAIN_SIZE = 10;
  private TransferFunction transferType;
  // gain values, indexed by stage
  private double[] gain;
  private int numStages;
  // poles and zeros
  private List<Pair<Complex, Integer>> zeros;
  private List<Pair<Complex, Integer>> poles;
  private String name;
  private Unit unitType;
  private double normalization; // A0 normalization factor
  private double normalFreq; // cuz she's a normalFreq, normalFreq
  private Instant epochStart, epochEnd;
  /**
   * Reads in a response from an already-accessed bufferedreader handle
   * and assigns it to the name given (used with embedded response files)
   * Only the last epoch of a multi-epoch response file is used.
   *
   * @param br Handle to a buffered reader of a given RESP file
   * @param name Name of RESP file to be used internally
   */
  private InstrumentResponse(BufferedReader br, String name) throws IOException {

    this.name = name;

    parserDriver(br, null);
  }

  /**
   * Create a copy of an existing response object
   *
   * @param responseIn The response object to be copied
   */
  public InstrumentResponse(InstrumentResponse responseIn) {
    transferType = responseIn.getTransferFunction();

    gain = responseIn.getGain();
    numStages = responseIn.getNumStages();

    zeros = new ArrayList<>(responseIn.getZerosList());
    poles = new ArrayList<>(responseIn.getPolesList());

    unitType = responseIn.getUnits();

    normalization = responseIn.getNormalization();
    normalFreq = responseIn.getNormalizationFrequency();

    name = responseIn.getName();
    epochStart = responseIn.getEpochStart();
    epochEnd = responseIn.getEpochEnd();
  }
  /**
   * Reads in an instrument response from a RESP file with the specific epoch
   * If the RESP file has multiple epochs, only the last one is used.
   *
   * @param filename full path of the RESP file
   */
  public InstrumentResponse(String filename, Instant epoch) throws IOException {
    name = new File(filename).getName();
    parseResponseFile(filename, epoch);
  }

  /**
   * Reads in an instrument response from a RESP file and gets the last epoch
   * If the RESP file has multiple epochs, only the last one is used.
   *
   * @param filename full path of the RESP file
   */
  public InstrumentResponse(String filename) throws IOException {
    this(filename, null);

  }

  /**
   * Get one of the response files embedded in the program
   *
   * @return response file embedded into the program
   * @throws IOException If no file with the given name exists (this may happen
   * if a file listed in the responses.txt file does not exist in that location
   * which means it was likely improperly modified or a response file deleted)
   */
  public static InstrumentResponse loadEmbeddedResponse(String fname)
      throws IOException {

    ClassLoader cl = InputPanel.class.getClassLoader();
    InputStream is = cl.getResourceAsStream("resps/" + fname);
    BufferedReader fr = new BufferedReader(new InputStreamReader(is));
    return new InstrumentResponse(fr, fname);
  }

  /**
   * Get list of all responses embedded into the program by iterating over the list of files
   *
   * @return Set of strings representing response filenames
   */
  public static Set<String> parseInstrumentList() {
    Set<String> respFilenames = new HashSet<>();
    ClassLoader cl = InstrumentResponse.class.getClassLoader();

    InputStream respRead = cl.getResourceAsStream("responses.txt");
    BufferedReader respBuff =
        new BufferedReader( new InputStreamReader(respRead) );

    try {
      String name;
      name = respBuff.readLine();
      while (name != null) {
        respFilenames.add(name);
        name = respBuff.readLine();
      }
      respBuff.close();
    } catch (IOException e2) {
      e2.printStackTrace();
    }
    return respFilenames;
  }

  /**
   * Take in a string representing the path to a RESP file and outputs the epochs contained in it.
   * Used to select a specific epoch for a given file via the GUI.
   *
   * @param filename Location of a given RESP file
   * @return List of epochs (pair of start and end instances)
   * @throws IOException If there is a failure to read the resp file
   */
  public static List<Pair<Instant, Instant>> getRespFileEpochs(String filename)
      throws IOException {
    try (BufferedReader br = new BufferedReader(new FileReader(filename))) {
      return getRespFileEpochs(br);
    }
  }

  /**
   * Returns the RESP epoch start time under which data of interest is held. If data does not exist
   * in an epoch defined in this RESP file it will either return the first or last epoch depending
   * on if the time region ends before the first epoch or starts after the last.
   * Note that this function assumes the given RESP file does not have a gap over the time region
   * of interest, and also assumes that both start and end exist within the same RESP epoch.
   * @param filename Name of RESP file to get epoch data from
   * @param start Start time of data region under interest (epoch milliseconds)
   * @param end End time of data region under interest
   * @return Enclosing epoch's start time, or first or last
   * @throws IOException if the file is unable to be read
   */
  public static Instant getRespFileClosestEpoch(String filename, long start, long end)
      throws IOException {
    List<Pair<Instant, Instant>> epochBounds;
    try (BufferedReader br = new BufferedReader(new FileReader(filename))) {
      epochBounds = getRespFileEpochs(br);

      // get initial best-guess estimate (first point), check if earliest epoch post-dates range
      Instant closestStart = epochBounds.get(0).getFirst();
      if (epochBounds.size() == 1 || end < closestStart.toEpochMilli()) {
        return closestStart;
      }

      long diff = Long.MAX_VALUE; // used to determine closest epoch start to inputted time
      for (Pair<Instant, Instant> oneEpoch : epochBounds) {
        Instant epochStart = oneEpoch.getFirst();
        Instant epochEnd = oneEpoch.getSecond();
        // if range falls inside resp epoch (expected) then use that one
        if (epochStart.toEpochMilli() < start && end < epochEnd.toEpochMilli()) {
          return epochStart;
        } else if (Math.abs(epochStart.toEpochMilli() - start) < diff) {
          // this is only of concern if the specified range doesn't somehow fall in a resp epoch
          closestStart = epochStart;
          diff = Math.abs(epochStart.toEpochMilli() - start);
        }
      }
      return closestStart;
    }
  }

  /**
   * Take in a file object and output the epochs contained in it. Used to select a specific
   * epoch in a given file from the GUI.
   *
   * @param br BufferedReader wrapper for a resp file
   * @return List of epochs (pair of start and end instances)
   * @throws IOException If there is a failure to read the resp file
   */
  private static List<Pair<Instant, Instant>> getRespFileEpochs(BufferedReader br)
      throws IOException {
    List<Pair<Instant, Instant>> epochList = new ArrayList<>();
    String line = br.readLine();

    while (line != null) {
      if (line.length() == 0) {
        // empty line? need to skip it
        line = br.readLine();
        continue;
      }

      if (line.charAt(0) != '#') {
        // not a comment - don't skip
        // words array is the components of each line, assuming split by 2 or more spaces
        String[] words = line.split("\\s\\s+");
        String hexIdentifier = words[0];
        switch (hexIdentifier) {
          case "B052F22":
            Instant start = parseTermAsDate(line);
            // read line to get end date
            line = br.readLine();
            Instant end;
            if (line.split("\\s+")[0].equals("B052F23")) {
              end = parseTermAsDate(line);
            } else {
              end = null;
            }
            epochList.add(new Pair<>(start, end));
            break;
        }
      }
      line = br.readLine();

    }

    return epochList;
  }

  /**
   * Read the line defining a date (epoch start or end) as an instant
   *
   * @param line Line of a response file that should define an epoch
   * @return Instant parsed from the given line or null if not able to parse
   */
  static Instant parseTermAsDate(String line) {
    // reparse the line
    String[] words = line.split("\\s+");
    // index 0 is the identifier for the field types (used in switch-stmt)
    // index 1 is the word 'Start'
    // index 2 is the word 'End'
    // index 3 is comma-delimited date string (should be last thing in list)
    String time = words[words.length - 1];
    // response epoch definition MUST have year and julian day
    StringBuilder format = new StringBuilder("uuuu,DDD");
    // length of time string too short -- much not have year and julian day
    if (time.length() < format.length()) {
      return null;
    }

    // now check to make sure at least year and julian day are formatted correctly
    DateTimeFormatter respDTFormat =
        DateTimeFormatter.ofPattern(format.toString()).withZone(ZoneOffset.UTC);
    String timeTest = time.substring(0, format.length());
    Instant bestParsedDate;
    try {
      bestParsedDate = LocalDate.parse(timeTest, respDTFormat).atStartOfDay().toInstant(ZoneOffset.UTC);
    } catch (DateTimeParseException malformed) {
      // this is expected to happen if this epoch has an open end time
      return null;
    }

    // can parse past this level -- does epoch define hours, minutes, seconds?
    // do it in a loop because all of these components are optional but, say,
    // minutes can't exist if hours doesn't, etc.
    String[] parseLevels = new String[]{",HH",":mm",":ss"};
    for (String nextLevel : parseLevels) {
      format.append(nextLevel);
      if (time.length() < format.length()) {
        break; // return result from previous stage in loop
      }
      // keep trimming -- also handles case where time string has milliseconds defined
      timeTest = time.substring(0, format.length());
      // datetimeformatter is the same, but the string defining the format has changed
      respDTFormat =
          DateTimeFormatter.ofPattern(format.toString()).withZone(ZoneOffset.UTC);
      try {
        bestParsedDate = LocalDateTime.parse(timeTest, respDTFormat).toInstant(ZoneOffset.UTC);
      } catch (DateTimeParseException invalidSublevel) {
        // also return best result from previous loop if the new stage has an error
        break;
      }
    }

    return bestParsedDate;
  }

  /**
   * Extract the real and imaginary terms from a pole or zero in a RESP file
   *
   * @param line the line the zero or pole is found on in the file
   * @param array the array of zeros and poles the term will be added to
   */
  static void parseTermAsComplex(String line, Complex[] array) {
    // reparse the line. why are we doing this? well,
    // if a number is negative, only one space between it and prev. number
    // and the previous split operation assumed > 2 spaces between numbers
    String[] words = line.split("\\s+");

    // index 0 is the identifier for the field types (used in switch-stmt)
    // index 1 is where in the list this zero or pole is
    // index 2 is the real part, and index 3 the imaginary
    // indices 4 and 5 are error terms (ignored)
    int index = Integer.parseInt(words[1]);
    double realPart = Double.parseDouble(words[2]);
    double imagPart = Double.parseDouble(words[3]);
    array[index] = new Complex(realPart, imagPart);
  }

  private static List<Pair<Complex, Integer>> setComponentValues(Complex[] pzArr) {
    // other response code assumes all poles/zeros are listed in order of increasing magnitude
    // (i.e., lower frequency poles/zeros are listed first, and then increasingly higher)
    NumericUtils.complexMagnitudeSorter(pzArr);
    // first, sort matching P/Z values into bins, count up the number of each
    Map<Complex, Integer> values = new LinkedHashMap<>();
    for (Complex c : pzArr) {
      if (values.keySet().contains(c)) {
        int count = values.get(c) + 1;
        values.put(c, count);
      } else {
        values.put(c, 1);
      }
    }
    // initialize list here
    List<Pair<Complex, Integer>> valueList = new ArrayList<>();
    // since second value in pair keeps track of item count in array, use set to prevent duplicates
    // (note that we don't want to store in a set because order matters for this data)
    Set<Complex> inList = new HashSet<>();
    // iterate from the array; it should have the values in the correct order
    for (Complex c : pzArr) {
      if (!inList.contains(c)) {
        Pair<Complex, Integer> toAdd = new Pair<>(c, values.get(c));
        valueList.add(toAdd);
        inList.add(c);
      }
    }
    return valueList;
  }

  /**
   * Skip to the desired epoch start Instant.
   * If it does not exist, the parser will throw an exception later.
   *
   * @param reader the shared BufferedReader with RESP file
   * @param epoch the desired epoch's start Instant.
   * @throws IOException on file read errors.
   */
  static void skipToSelectedEpoch(BufferedReader reader, Instant epoch) throws IOException {
    String line = reader.readLine();
    while (line != null) {
      if (line.startsWith("B052F22")) {
        Instant start = parseTermAsDate(line);
        if (epoch.equals(start)) {
          break;
        }
      }
      line = reader.readLine();
    }
  }

  static TransferFunction parseTransferType(String word) {
    switch (word.charAt(0)) {
      case 'A':
        return TransferFunction.LAPLACIAN;
      case 'B':
        return TransferFunction.LINEAR;
      default:
        // defaulting to LAPLACIAN if type is different from a or b
        // which is likely to be more correct
        return TransferFunction.LAPLACIAN;
    }
  }

  static Unit parseUnitType(String unit) throws IOException {
    switch (unit.toLowerCase()) {
      case "m/s":
        return Unit.VELOCITY;
      case "m/s**2":
        return Unit.ACCELERATION;
      default:
        String e = "Unit type was given as " + unit + ".\n";
        e += "Nonstandard unit, or not a velocity or acceleration";
        throw new IOException(e);
    }
  }

  /**
   * Apply the values of this response object to a list of frequencies and
   * return the resulting (complex) frequencies
   * The response curve produced is in units of velocity. Some results
   * will need to have the produced response curve have acceleration units,
   * which can be done by multiplying by the integration factor defined here.
   *
   * @param frequencies inputted list of frequencies, such as FFT windows
   * @return application of the response to those frequencies
   */
  public Complex[] applyResponseToInput(double[] frequencies) {

    Complex[] resps = new Complex[frequencies.length];

    // precalculate gain for scaling the response
    double scale = 1.;
    // stage 0 is sensitivity (supposed to be product of all gains)
    // we will get scale by multiplying all gain stages except for it
    for (int i = 1; i < gain.length; ++i) {
      scale *= gain[i];
    }

    // how many times do we need to do differentiation?
    // outUnits (acceleration) - inUnits
    // i.e., if the units of this response are acceleration, we integrate once
    int diffs = Unit.VELOCITY.getDifferentiations(unitType);
    // unlike s (see below) this is always 2Pi
    double integConstant = NumericUtils.TAU;

    for (int i = 0; i < frequencies.length; ++i) {
      double deltaFrq = frequencies[i];

      // pole-zero expansion
      Complex s = new Complex(0, deltaFrq * transferType.getFunction());

      Complex numerator = Complex.ONE;
      Complex denominator = Complex.ONE;

      // generate the list structures and get the values with respect to transfer function
      // unfolding the paired list like this is as reasonable as most other options
      for (Complex zero : getZeros()) {
        numerator = numerator.multiply(s.subtract(zero));
      }

      for (Complex pole : getPoles()) {
        denominator = denominator.multiply(s.subtract(pole));
      }

      resps[i] = numerator.multiply(normalization).divide(denominator);

      if (diffs < 0) {
        // a negative number of differentiations
        // is a positive number of integrations
        // i*omega; integration is I(w) x (iw)^n
        Complex iw = new Complex(0.0, integConstant * deltaFrq);
        // start at 1 in these loops because we do mult. at least once
        for (int j = 1; j < Math.abs(diffs); j++) {
          iw = iw.multiply(iw);
        }
        resps[i] = resps[i].multiply(iw);
      } else if (diffs > 0) {
        // differentiation is I(w) / (-i/w)^n
        Complex iw = new Complex(0.0, -1.0 / (integConstant * deltaFrq));
        for (int j = 1; j < Math.abs(diffs); j++) {
          iw = iw.multiply(iw);
        }
        resps[i] = iw.multiply(resps[i]);
      }

      // lastly, scale by the scale we chose (gain0 or gain1*gain2)
      resps[i] = resps[i].multiply(scale);
    }

    return resps;
  }

  /**
   * Given a best-fit vector, build the poles and zeros to use the ones
   * defined by that vector. Imaginary values that are non-zero are constrained
   * to be (implicitly) defining a complex conjugate of the pole/zero they
   * are built into. Similarly, any poles or zeros in the initial list that
   * are duplicated are all set to have the same new value.
   *
   * @param params Array of real and imaginary component values of poles
   * and zeros
   * @param lowFreq True if the fit values are for low-frequency components
   * @param numZeros How much of the input parameter array is zero components
   * @return New InstrumentResponse with the fit values applied to it
   * @see #zerosToVector and #polesToVector
   */
  public InstrumentResponse
  buildResponseFromFitVector(double[] params, boolean lowFreq,
      int numZeros) {

    // get all distinct pole values as lists
    // List<Complex> zList = getDistinctZeros();
    // List<Complex> pList = getDistinctPoles();

    // first covert poles and zeros back to complex values to make this easier
    List<Complex> zerosAsComplex = new ArrayList<>();
    for (int i = 0; i < numZeros; i += 2) {
      Complex c = new Complex(params[i], params[i + 1]);
      zerosAsComplex.add(c);
    }

    List<Complex> polesAsComplex = new ArrayList<>();
    for (int i = numZeros; i < params.length; i += 2) {
      Complex c = new Complex(params[i], params[i + 1]);
      polesAsComplex.add(c);
    }

    // fit the zeros
    List<Pair<Complex, Integer>> builtZeros = new ArrayList<>();

    // first, add the literally zero values; these are never fit
    // (NOTE: we expect count to never be more than 2)
    int start;
    for (start = 0; start < zeros.size(); ++start) {
      Pair<Complex, Integer> valueAndCount = zeros.get(start);
      Complex c = valueAndCount.getFirst();
      // if it's not zero, it might need to get replaced
      if (c.abs() > 0) {
        break;
      }
      builtZeros.add(valueAndCount);
    }

    // add the low-frequency zeros from source if they're not being fit
    if (!lowFreq) {
      // add zeros until they reach the high-freq cutoff point
      // start from current index of data
      for (int i = start; i < zeros.size(); ++i) {
        Pair<Complex, Integer> valueAndCount = zeros.get(i);
        Complex zero = valueAndCount.getFirst();
        if (zero.abs() / NumericUtils.TAU > 1.) {
          // zeros after this point are high-frequency
          break;
        }
        builtZeros.add(valueAndCount);
      }
    }

    // now add the zeros under consideration for fit
    // these are the high-frequency zeros if we're doing high-frequency cal
    // or the low-frequency zeros otherwise
    int offset;
    offset = builtZeros.size();
    for (int i = 0; i < zerosAsComplex.size(); ++i) {
      // get the number of times the original value appeared
      int count = zeros.get(i + offset).getSecond();
      Complex zero = zerosAsComplex.get(i);
      builtZeros.add(new Pair<>(zero, count));

      // add conjugate if it has one
      if (zero.getImaginary() != 0.) {
        builtZeros.add(new Pair<>(zero.conjugate(), count));
        ++offset; // skipping over the original conjugate pair
      }
    }

    // now add in all remaining zeros
    for (int i = builtZeros.size(); i < zeros.size(); ++i) {
      builtZeros.add(zeros.get(i));
    }

    // now do the same thing as the zeros but for the poles
    List<Pair<Complex, Integer>> builtPoles = new ArrayList<>();

    // low frequency poles not being fit added first (keeps list sorted)
    if (!lowFreq) {
      // first add low-frequency poles not getting fit by high-freq cal
      for (Pair<Complex, Integer> valueAndCount : poles) {
        Complex pole = valueAndCount.getFirst();
        if (pole.abs() / NumericUtils.TAU > 1.) {
          break;
        }
        builtPoles.add(valueAndCount);
      }
    } else if (isKS54000()) {
      // used in the odd KS54000 case, we don't fit the low-freq damping pole
      // there should only be the one, because the KS54000 is a weird one
      builtPoles.add(poles.get(0));
    }

    offset = builtPoles.size();
    // now add the poles under consideration for fit as with zeros
    for (int i = 0; i < polesAsComplex.size(); ++i) {
      int count = poles.get(i + offset).getSecond();
      Complex pole = polesAsComplex.get(i);
      builtPoles.add(new Pair<>(pole, count));

      // add conjugate if it has one
      if (pole.getImaginary() != 0.) {
        builtPoles.add(new Pair<>(pole.conjugate(), count));
        ++offset;
      }
    }

    // now add the poles that remain
    for (int i = builtPoles.size(); i < poles.size(); ++i) {
      builtPoles.add(poles.get(i));
    }

    // create a copy of this instrument response and set the new values
    InstrumentResponse out = new InstrumentResponse(this);
    out.setZerosPairs(builtZeros);
    out.setPolesList(builtPoles);
    return out;

  }

  /**
   * Get the gain stages of the RESP file. Stage x is at index x. That is,
   * the sensitivity is at 0, the sensor gain is at 1, and the digitizer
   * gain is at 2.
   *
   * @return Array of all gain stages found in resp file, including stage 0
   */
  public double[] getGain() {
    return gain;
  }

  /**
   * Return the name of this response file (i.e., STS-1Q330 or similar)
   * Used primarily for identifying response curve on plots
   *
   * @return String containing ID of current response
   */
  public String getName() {
    return name;
  }

  /**
   * Set name of response file, used in some plot and report generation
   *
   * @param newName New name to give this response
   */
  public void setName(String newName) {
    name = newName;
  }

  /**
   * Get the normalization of the response
   *
   * @return normalization constant
   */
  public double getNormalization() {
    return normalization;
  }

  /**
   * Get the normalization frequency
   *
   * @return normalization frequency (Hz)
   */
  public double getNormalizationFrequency() {
    return normalFreq;
  }

  /**
   * Return the number of gain stages this response has
   *
   * @return number of gain stages
   */
  public int getNumStages() {
    return numStages;
  }

  /**
   * Return the list of poles in the RESP file, not including error terms
   *
   * @return List of complex numbers; index y is the yth pole in response list
   */
  public List<Complex> getPoles() {
    List<Complex> out = new ArrayList<>();
    for (Pair<Complex, Integer> valueAndCount : poles) {
      Complex value = valueAndCount.getFirst();
      Integer count = valueAndCount.getSecond();
      // count is number of times pole appears, add it that many times to the list
      for (int i = 0; i < count; ++i) {
        out.add(value);
      }
    }
    return out;
  }

  /**
   * Replace the current poles of this response with new ones from a list
   *
   * @param poleList List of poles to replace the current response poles with (repeated poles listed
   * the number of times they appear)
   */
  public void setPoles(List<Complex> poleList) {
    Complex[] poleArr = poleList.toArray(new Complex[]{});
    setPoles(poleArr);
  }

  private List<Pair<Complex, Integer>> getPolesList() {
    return poles;
  }

  /**
   * Apply new poles to this response object
   *
   * @param newPoles List of paired values; first entry is a unique pole in response, and
   * second is the number of times it appears
   */
  private void setPolesList(List<Pair<Complex, Integer>> newPoles) {
    poles = newPoles;
  }

  /**
   * Get the transfer function of this response file (laplacian, linear)
   *
   * @return transfer type as an enumeration (can get factor as numeric type
   * by calling getFunction() on the returned value)
   */
  public TransferFunction getTransferFunction() {
    return transferType;
  }

  /**
   * Gives the unit type of the RESP file (displacement, velocity, acceleration)
   *
   * @return Unit type as enumeration, distance measures of meters and time of
   * seconds (i.e., acceleration is m/s^2)
   */
  public Unit getUnits() {
    return unitType;
  }

  /**
   * Return the list of zeros in the RESP file, not including error terms
   *
   * @return List of complex numbers; index y is the yth zero in response list
   */
  public List<Complex> getZeros() {
    List<Complex> out = new ArrayList<>();
    for (Pair<Complex, Integer> valueAndCount : zeros) {
      Complex value = valueAndCount.getFirst();
      Integer count = valueAndCount.getSecond();
      // count is number of times zero appears, add it that many times to the list
      for (int i = 0; i < count; ++i) {
        out.add(value);
      }
    }
    return out;
  }

  private List<Pair<Complex, Integer>> getZerosList() {
    return zeros;
  }

  /**
   * Determines whether or not the first pole is too low to fit even with
   * a low-frequency random cal solver operation. Mainly an issue for
   * the KS54000 RESP, since the damping curve is unusual.
   *
   * @return True if the initial pole is too low-frequency to add to fit.
   */
  private boolean isKS54000() {
    final double CUTOFF = 1. / 1000.;
    List<Complex> pList = getPoles();
    return (pList.get(0).abs() / NumericUtils.TAU) < CUTOFF;
  }

  /**
   * Read in each line of a response and parse and store relevant lines
   * according to the hex value at the start of the line
   *
   * @param reader reader of a given file to be parse
   * @param epoch starting Instant of epoch that is needed
   * @throws IOException if the reader cannot read the given file
   */
  private void parserDriver(BufferedReader reader, Instant epoch) throws IOException {
    // read in lines until the epoch is found
    if (epoch != null) {
      skipToSelectedEpoch(reader, epoch);
      epochStart = epoch;
    }

    numStages = 0;
    double[] gains = new double[INIT_GAIN_SIZE];
    for (int i = 0; i < gains.length; ++i) {
      gains[i] = 1;
    }
    normalization = 0;
    normalFreq = 0;
    int gainStage = -1;
    Complex[] polesArr = null;
    Complex[] zerosArr = null;

    String line = reader.readLine();

    epochLoop:
    while (line != null) {
      if (line.length() != 0 && line.charAt(0) != '#') {
        // the components of each line, assuming split by 2 or more spaces
        String[] words = line.split("\\s\\s+");
        String hexIdentifier = words[0];

        switch (hexIdentifier) {
          case "B052F22":
            if (epoch != null) {
              // we already read to the desired epoch, so we can just return here
              // otherwise, go for the last epoch and reset out the data from prev. epochs
              break epochLoop;
            }
            epochStart = parseTermAsDate(line);
            numStages = 0;
            gains = new double[INIT_GAIN_SIZE];
            for (int i = 0; i < gains.length; ++i) {
              gains[i] = 1;
            }
            normalization = 0;
            normalFreq = 0;
            gainStage = -1;
            polesArr = null;
            zerosArr = null;
            break;
          case "B052F23":
            epochEnd = parseTermAsDate(line);
            break;
          case "B053F03":
            // transfer function type specified
            // first character of third component of words
            transferType = parseTransferType(words[2]);
            break;
          case "B053F05":
            // parse the units of the transfer function (usually velocity)
            // first *word* of the third component of words
            String[] unitString = words[2].split("\\s");
            unitType = parseUnitType(unitString[0]);
            break;
          case "B053F07":
            // this is the normalization factor A0
            // this is the entire third word of the line, as a double
            normalization = Double.parseDouble(words[2]);
            break;
          case "B053F08":
            // this is the normalization frequency
            // once again the entire third word of the line as double
            normalFreq = Double.parseDouble(words[2]);
            break;
          case "B053F09":
            // the number of zeros listed in reponse pole/zero lines
            // again, this is the entire third word, as an int
            int numZero = Integer.parseInt(words[2]);
            zerosArr = new Complex[numZero];
            break;
          case "B053F14":
            // same as above line but for the number of poles
            int numPole = Integer.parseInt(words[2]);
            polesArr = new Complex[numPole];
            break;
          case "B053F10-13":
            // these are the lists of response zeros, in order with index,
            // real (double), imaginary (double), & corresponding error terms
            parseTermAsComplex(line, zerosArr);
            break;
          case "B053F15-18":
            // as above but for poles
            parseTermAsComplex(line, polesArr);
            break;
          case "B058F03":
            // gain stage sequence number; again, full third word as int
            // this is used to map the gain value to an index
            gainStage = Integer.parseInt(words[2]);
            numStages = Math.max(numStages, gainStage);
            break;
          case "B058F04":
            // should come immediately and only after the gain sequence number
            // again, it's the third full word, this time as a double
            // map allows us to read in the stages in whatever order
            // in the event they're not sorted in the response file
            // and allows us to have basically arbitrarily many stages
            if (gainStage >= gains.length) {
              double[] temp = new double[gains.length * 2];
              for (int i = 0; i < temp.length; ++i) {
                if (i < gains.length) {
                  temp[i] = gains[i];
                } else {
                  temp[i] = 1;
                }
              }
              gains = temp;
            }
            gains[gainStage] = Double.parseDouble(words[2]);

            // reset the stage to prevent data being overwritten
            gainStage = -1;
            break;
        }
      }
      line = reader.readLine();
    }


    ++numStages; // offset by 1 to represent size of stored gain stages
    gain = Arrays.copyOfRange(gains, 0, numStages);

    // turn pole/zero arrays into maps from pole values to # times repeated
    setZerosFromComplex(zerosArr);
    setPoles(polesArr);
  }

  /**
   * Parses a response file of the sort found on the Iris Nominal Response
   * Library. These files can be found at http://ds.iris.edu/NRL/
   * This function currently does not parse a full response file, but instead
   * only examines fields relevant to self-noise calculations.
   *
   * @param filename Full path to the response file
   */
  private void parseResponseFile(String filename, Instant epoch) throws IOException {

    // response files have a very nice format that is not so nice as something
    // like JSON but still quite easy to parse
    // lines either begin with a hex value or a '#'
    // '#' marks comments while the hex value is used to interpret values
    // each line with a hex value appears to have the following format
    // [hex value] [whitespace] [name] [whitespace] [value]
    // where name is the human-readable explanation of what a value represents
    // in some cases, 'value' may not be just a raw value but also include
    // some information about the value, such as verbose unit specifications

    // there is one exception, the actual pole/zero fields, which have 5
    // components after the hex identifier

    try (BufferedReader br = new BufferedReader(new FileReader(filename))) {
      parserDriver(br, epoch);
    }
  }

  /**
   * Create a RealVector of the components in the list of poles for this
   * response, used for finding best-fit poles from a random calibration
   * The complex numbers each representing a pole are split into their
   * real and imaginary components (each a double), with the real components
   * being set on even indices and imaginary components on the odd indices;
   * each even-odd pair (i.e., 0-1, 2-3, 4-5) of values in the vector define a
   * single pole. Poles with non-zero imaginary components do not have their
   * conjugate included in this vector in order to maintain constraints.
   *
   * @param lowFreq True if the low-frequency poles are to be fit
   * @param peak Peak frequency of zero to add to fit set (upper bound on high-freq poles to fit)
   * @return RealVector with fittable pole values
   */
  public RealVector polesToVector(boolean lowFreq, double peak) {
    // create a list of doubles that are the non-conjugate elements from list
    // of poles, to convert to array and then vector format
    List<Double> componentList = new ArrayList<>();

    // starting index for poles, shift up one if lowest-freq pole is TOO low
    int start = 0;
    if (isKS54000()) {
      start = 1;
    }

    for (int i = start; i < poles.size(); ++i) {
      Pair<Complex, Integer> pairAndValue = poles.get(i);
      Complex p = pairAndValue.getFirst();

      double frq = p.abs() / NumericUtils.TAU;
      if (!lowFreq && (frq < 1.)) {
        // don't include poles below 1Hz in high-frequency calibration
        continue;
      }
      if (lowFreq && (frq > 1.)) {
        // only do low frequency calibrations on poles up to
        break;
      }
      if (!lowFreq && (frq >= peak)) {
        // don't fit poles above fraction of nyquist rate of sensor output
        break;
      }

      // a complex is just two doubles representing real and imaginary lengths
      double realPart = p.getReal();
      double imagPart = p.getImaginary();

      componentList.add(realPart);
      componentList.add(imagPart);

      if (imagPart != 0.) {
        ++i; // skip complex conjuagte
      }
    }

    // turn into array to be turned into vector
    // can't use toArray because List doesn't use primitive double objects
    double[] responseVariables = new double[componentList.size()];
    for (int i = 0; i < responseVariables.length; ++i) {
      responseVariables[i] = componentList.get(i);
    }

    return MatrixUtils.createRealVector(responseVariables);
  }

  /**
   * Replace the current poles of this response with new ones from an array
   *
   * @param poleList Array of poles to replace the current response poles with (repeated poles
   * listed by the number of times they appear)
   */
  private void setPoles(Complex[] poleList) {
    poles = setComponentValues(poleList);
  }

  /**
   * Set the list of zeros to a new array, such as after fitting from random cal
   *
   * @param zeroList Array of zeros to replace the current response poles with (repeated zeros
   * listed every time they appear)
   */
  private void setZerosFromComplex(Complex[] zeroList) {
    zeros = setComponentValues(zeroList);
  }

  /**
   * Apply new zeros to this response object
   *
   * @param newZeros List of paired values; first entry is a unique zero in response, and
   * second is the number of times it appears
   */
  private void setZerosPairs(List<Pair<Complex, Integer>> newZeros) {
    zeros = newZeros;
  }

  /**
   * Output text report of this response file. Not same format as IRIS RESP.
   */
  @Override
  public String toString() {
    StringBuilder sb = new StringBuilder();

    NumberFormat nf = NumberFormat.getInstance();
    nf.setMaximumFractionDigits(4);
    ComplexFormat cf = new ComplexFormat(nf);

    sb.append("Response name: ");
    sb.append(name);
    sb.append('\n');
    sb.append("Epoch start: ");
    sb.append(RESP_DT_FORMAT.format(epochStart));
    sb.append('\n');
    if (epochEnd != null) {
      sb.append("Epoch end: ");
      sb.append(RESP_DT_FORMAT.format(epochEnd));
      sb.append('\n');
    }

    sb.append("Gain stage values: ");
    sb.append('\n');

    for (int i = 0; i < numStages; ++i) {
      sb.append(i);
      sb.append(": ");
      sb.append(nf.format(gain[i]));
      sb.append("\n");
    }

    sb.append("Normalization: ");
    sb.append(normalization);
    sb.append('\n');
    sb.append("Normalization frequency (Hz): ");
    sb.append(normalFreq);
    sb.append('\n');

    sb.append("Transfer function ");
    if (transferType == TransferFunction.LAPLACIAN) {
      sb.append("is LAPLACIAN");
    } else {
      sb.append("is LINEAR");
    }
    sb.append('\n');

    sb.append("Response input units: ");
    switch (unitType) {
      case DISPLACEMENT:
        sb.append("displacement (m)");
        break;
      case VELOCITY:
        sb.append("velocity (m/s)");
        break;
      case ACCELERATION:
        sb.append("acceleration (m/s^2)");
        break;
    }
    sb.append('\n');

    sb.append("Response zeros: ");
    sb.append('\n');
    List<Complex> zList = getZeros();
    for (int i = 0; i < zList.size(); ++i) {
      sb.append(i);
      sb.append(": ");
      sb.append(cf.format(zList.get(i)));
      sb.append("\n");
    }

    sb.append("Response poles: ");
    sb.append('\n');
    List<Complex> pList = getPoles();
    for (int i = 0; i < pList.size(); ++i) {
      sb.append(i);
      sb.append(": ");
      sb.append(cf.format(pList.get(i)));
      sb.append("\n");
    }

    return sb.toString();
  }

  /**
   * Get the start time of the epoch of this response data
   *
   * @return Epoch expressed as an instant
   */
  public Instant getEpochStart() {
    return epochStart;
  }

  /**
   * Get the end time of the epoch of this response data
   *
   * @return Epoch expressed as an instant (can be null)
   */
  Instant getEpochEnd() {
    return epochEnd;
  }

  /**
   * Create a RealVector of the components in the list of zeros for this
   * response, used for finding best-fit zeros from a random calibration
   * The complex numbers each representing a zero are split into their
   * real and imaginary components (each a double), with the real components
   * being set on even indices and imaginary components on the odd indices;
   * each even-odd pair (i.e., 0-1, 2-3, 4-5) of values in the vector define a
   * single zero. Zeros with non-zero imaginary components do not have their
   * conjugate included in this vector in order to maintain constraints.
   *
   * @param lowFreq True if the low-frequency zeros are to be fit
   * @param peak Peak frequency of zero to add to fit set (upper bound on high-freq zeros to fit)
   * @return RealVector with fittable zero values
   */
  public RealVector zerosToVector(boolean lowFreq, double peak) {

    // create a list of doubles that are the non-conjugate elements from list
    // of poles, to convert to array and then vector format
    List<Double> componentList = new ArrayList<>();

    for (int i = 0; i < zeros.size(); ++i) {

      Pair<Complex, Integer> valueAndCount = zeros.get(i);
      Complex c = valueAndCount.getFirst();

      if (c.abs() == 0.) {
        // ignore zeros that are literally zero-valued
        continue;
      }

      double cutoffChecker = c.abs() / NumericUtils.TAU;

      if (lowFreq && (cutoffChecker > 1.)) {
        // only do low frequency calibrations on zeros up to 1Hz
        break;
      }
      if (!lowFreq && (cutoffChecker < 1.)) {
        // don't include zeros > 1Hz in high-frequency calibration
        continue;
      }
      if (!lowFreq && (cutoffChecker > peak)) {
        // don't fit zeros above 80% nyquist rate of sensor output
        break;
      }

      double realPart = c.getReal();
      double imagPart = c.getImaginary();
      componentList.add(realPart);
      componentList.add(imagPart);

      if (imagPart != 0.) {
        ++i;
      }

    }

    // turn into array to be turned into vector
    // can't use toArray because List doesn't use primitive double objects
    double[] responseVariables = new double[componentList.size()];
    for (int i = 0; i < responseVariables.length; ++i) {
      responseVariables[i] = componentList.get(i);
    }

    return MatrixUtils.createRealVector(responseVariables);
  }

}


