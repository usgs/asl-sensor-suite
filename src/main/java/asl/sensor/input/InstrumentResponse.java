package asl.sensor.input;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.complex.Complex;

/**
 * This class is used to read in and store data from instrument response files
 * such as those found on the IRIS RESP database. This includes functions used
 * to parse in relevant data from those files to get values such as poles and
 * zeros and the gain stages of the detector/sensor setup.
 * See also: http://ds.iris.edu/ds/nodes/dmc/data/formats/resp/
 * @author akearns
 *
 */
public class InstrumentResponse {

  
  /**
   * Extract the real and imaginary terms from a pole or zero in a RESP file
   * @param line the line the zero or pole is found on in the file
   * @param array the array of zeros and poles the term will be added to
   */
  private static void parseTermAsComplex(String line, Complex[] array) {
    // reparse the line. why are we doing this?
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
  
  private TransferFunction transferType;
  
  // gain values, indexed by stage
  private List<Double> gain;
  private List<Complex> zeros;
  
  private List<Complex> poles;
  
  private String name;
  
  private Unit unitType;
  private double normalization;
  
  private double normalFreq; // ♫ cuz she's a normalFreq, normalFreq 
  
  /**
   * Create a copy of an existing response object
   * @param responseIn The response object to be copied
   */
  public InstrumentResponse(InstrumentResponse responseIn) {
    transferType = responseIn.getTransferFunction();
    
    gain = new ArrayList<Double>( responseIn.getGain() );

    zeros = new ArrayList<Complex>( responseIn.getZeros() );
    poles = new ArrayList<Complex>( responseIn.getPoles() );
    
    unitType = responseIn.getUnits();
    
    normalization = Double.valueOf( responseIn.getNormalization() );
    normalFreq = Double.valueOf( responseIn.getNormalizationFrequency() );
    
    name = responseIn.getName();
  }
  
  /**
   * Reads in an instrument response from a RESP file
   * @param filename full path of the RESP file
   * @throws IOException
   */
  public InstrumentResponse(String filename) throws IOException {
    name = new File(filename).getName();
    parseResponseFile(filename);
  }
  
  /**
   * Reads in a response from an already-accessed bufferedreader handle
   * and assigns it to the name given
   * @param br Handle to a buffered reader of a given RESP file
   * @param name Name of RESP file to be used internally 
   * @throws IOException
   */
  public InstrumentResponse(BufferedReader br, String name) throws IOException {
    
    this.name = name;
    
    parserDriver(br);
  }
  
  /**
   * Return the name of this response file (i.e., STS-1Q330 or similar)
   * Used primarily for identifying response curve on plots
   * @return String containing ID of current response
   */
  public String getName() {
    return name;
  }
  
  /**
   * Apply the values of this response object to a list of frequencies and
   * return the resulting (complex) frequencies
   * @param frequencies inputted list of frequencies, such as FFT windows
   * @return application of the response to those frequencies
   */
  public Complex[] applyResponseToInput(double[] frequencies) {
   
    Complex[] resps = new Complex[frequencies.length];
    
    // precalculate gain for scaling the response
    double scale = 1.;
    // stage 0 is sensitivity (supposed to be product of all gains)
    // we will get scale by multiplying all gain stages except for it
    for (int i = 1; i < gain.size(); ++i) {
      scale *= gain.get(i);
    }
    
    // how many times do we need to do integration?
    // outUnits (acceleration) - inUnits
    // i.e., if the units of this response are velocity, we integrate once
    int integrations = Unit.ACCELERATION.getDifferentiations(unitType);
    // unlike s (see below) this is always 2Pi
    double integConstant = 2*Math.PI;
    
    for (int i = 0; i < frequencies.length; ++i) {
     
      double deltaFrq = frequencies[i];
      
      // pole-zero expansion
      Complex s = new Complex( 0, deltaFrq*transferType.getFunction() );
      
      Complex numerator = Complex.ONE;
      Complex denominator = Complex.ONE;
      
      for (Complex zero : zeros) {
        numerator = numerator.multiply( s.subtract(zero) );
      }
      
      for (Complex pole : poles) {
        denominator = denominator.multiply( s.subtract(pole) );
      }
      
      resps[i] = numerator.multiply(normalization).divide(denominator);
      
      if (integrations > 0) {
        // i*omega; integration is I(w) x (iw)^n
        Complex iw = new Complex(0.0, integConstant*deltaFrq);
        for (int j = 1; j < Math.abs(integrations); j++){
          iw = iw.multiply(iw);
        }
        resps[i] = resps[i].multiply(iw);
      } else if (integrations < 0) { 
        // a negative number of integrations 
        // is a positive number of differentiations
        // differentiation is I(w) / (-i/w)^n
        Complex iw = new Complex(0.0, -1.0 / (integConstant*deltaFrq) );
        for (int j = 1; j < Math.abs(integrations); j++){
          iw = iw.multiply(iw);
        }
        resps[i] = iw.multiply(resps[i]);
      }
      
      
      // lastly, scale by the scale we chose (gain0 or gain1*gain2)
      resps[i] = resps[i].multiply(scale);
      
      // unit conversion back out of acceleration
      /*
      double conversion =  frequencies[i]*2*Math.PI;
      resps[i] = resps[i].multiply( resps[i].conjugate() );
      resps[i] = resps[i].multiply( Math.pow(conversion, 2) );
      */
      
    }
    
    return resps;
  }
  
  /**
   * Get the gain stages of the RESP file. Stage x is at index x. That is,
   * the sensitivity is at 0, the sensor gain is at 1, and the digitizer
   * gain is at 2.
   * @return
   */
  public List<Double> getGain() {
    return gain;
  }
  
  /**
   * Get the normalization of the response
   * @return
   */
  public double getNormalization() {
    return normalization;
  }
  
  /**
   * Get the normalization frequency
   * @return
   */
  public double getNormalizationFrequency() {
    return normalFreq;
  }
  
  /**
   * Return the list of poles in the RESP file, not including error terms
   * @return List of complex numbers; index y is the yth pole in response list
   */
  public List<Complex> getPoles() {
    return poles;
  }
  
  /**
   * Get the transfer function of this response file (laplacian, linear)
   * @return
   */
  public TransferFunction getTransferFunction() {
    return transferType;
  }
  
  /**
   * Gives the unit type of the RESP file (displacement, velocity, acceleration)
   * @return
   */
  public Unit getUnits() {
    return unitType;
  }
  
  /**
   * Return the list of zeros in the RESP file, not including error terms
   * @return List of complex numbers; index y is the yth zero in response list
   */
  public List<Complex> getZeros() {
    return zeros;
  }
  
  /**
   * Parses a response file of the sort found on the Iris Nominal Response
   * Library. These files can be found at http://ds.iris.edu/NRL/
   * This function currently does not parse a full response file, but instead
   * only examines fields relevant to self-noise calculations.
   * @param filename Full path to the response file
   */
  private void parseResponseFile(String filename) throws IOException {
    
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
    
    BufferedReader br;
    try {
      br = new BufferedReader( new FileReader(filename) );
      parserDriver(br);
      br.close();
    } catch (FileNotFoundException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
    
  }
  
  private void parserDriver(BufferedReader br) throws IOException {
    
    String line = br.readLine();
    
    // <gain stage, gain value>
    Map<Integer, Double> gainMap = new HashMap<Integer, Double>();
    int gainStage = -1;
    Complex[] polesArr = null;
    Complex[] zerosArr = null;
    
    while (line != null) {
      
      if( line.length() == 0 ) {
        // empty line? need to skip it
        line = br.readLine();
        continue;
      }
      
      if (line.charAt(0) == '#') {
        // comment -- skip
        line = br.readLine();
        continue;
      } else {
        // the components of each line, assuming split by 2 or more spaces
        String[] words = line.split("\\s\\s+");
        String hexIdentifier = words[0];
        
        switch (hexIdentifier) {
        case "B053F03":
          // transfer function type specified
          // first character of third component of words
          switch ( words[2].charAt(0) ) {
          case 'A':
            transferType = TransferFunction.LAPLACIAN;
            break;
          case 'B':
            transferType = TransferFunction.LINEAR;
            break;
          default:
            // TODO: want to make sure that we're setting this right
            // because if a Fourier is specified, we should set to LAPLACIAN
            transferType = TransferFunction.LAPLACIAN;
          }
          break;
        case "B053F05":
          // parse the units of the transfer function (usually velocity)
          // first *word* of the third component of words
          // TODO support more unit types
          String[] unitString = words[2].split("\\s");
          String unit = unitString[0];
          switch (unit.toLowerCase()) {
          case "m/s":
            unitType = Unit.VELOCITY;
            break;
          case "m/s**2":
            unitType = Unit.ACCELERATION;
            break;
          default:
            String e = "Unit type was given as " + unit + ".\n";
            e += "Nonstandard unit, or not a velocity or acceleration";
            throw new IOException(e);
          }
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
          parseTermAsComplex(line,zerosArr);
          break;
        case "B053F15-18":
          // as above but for poles
          parseTermAsComplex(line, polesArr);
          break;
        case "B058F03":
          // gain stage sequence number; again, full third word as int
          // this is used to map the gain value to an index
          gainStage = Integer.parseInt(words[2]);
          break;
        case "B058F04":
          
          // should come immediately and only after the gain sequence number
          // again, it's the third full word, this time as a double
          // map allows us to read in the stages in whatever order
          // in the event they're not sorted in the response file
          // and allows us to have basically arbitrarily many stages
          gainMap.put( gainStage, Double.parseDouble(words[2]) );
          
          // reset the stage to prevent data being overwritten
          gainStage = -1;
          break;
        }
        
        line = br.readLine();
      } // else
      
    } // end of file-read loop (EOF reached, line is null)
    
    // turn map of gain stages into list
    List<Integer> stages = new ArrayList<Integer>( gainMap.keySet() );
    Collections.sort( stages );
    gain = new ArrayList<Double>();
    for (int stage : stages) {
      gain.add( gainMap.get(stage) );
    }
    // turn pole/zero arrays into lists
    zeros = Arrays.asList(zerosArr);
    poles = Arrays.asList(polesArr);
  }

  public void setPoles(List<Complex> poleList) {
    poles = poleList;
  }
  
}


