package asl.sensor.input;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import org.apache.log4j.Logger;

/**
 * Data representing a noise model for station data. These are expected to be a CSV with
 * each line having a period value, mean noise value, median noise value, 10th percentile value, and
 * 90th percentile value.
 */
public class NoiseModel {

  public double[] periods;
  public double[] meanLevels;
  public double[] medianLevels;
  public double[] percentile10;
  public double[] percentile90;
  public String name;

  Logger logger = Logger.getLogger(NoiseModel.class);

  /**
   * Construct a new noise model from a CSV file
   * @param filename CSV file to read in
   * @throws IOException If the file cannot be read
   * @throws NoiseModelFormatException If the file does not have the expected field number
   */
  public NoiseModel (String filename) throws IOException, NoiseModelFormatException {
    this(new File(filename));
  }

  /**
   * Construct a new noise model from a CSV file
   * @param file CSV file to read in
   * @throws IOException If the file cannot be read
   * @throws NoiseModelFormatException If the file does not have the expected field number
   */
  public NoiseModel (File file) throws IOException, NoiseModelFormatException {
    name = file.getName();
    List<Double> periods = new ArrayList<>();
    List<Double> means = new ArrayList<>();
    List<Double> medians = new ArrayList<>();
    List<Double> tenths = new ArrayList<>();
    List<Double> ninetieths = new ArrayList<>();
    int lineLength = 0;

    try (BufferedReader br = new BufferedReader(new FileReader(file))) {
      String line = br.readLine();
      String[] args = line.trim().split("\\s+");
      lineLength = args.length;


      // if this first line is a format description header skip it
      try {
        Double.valueOf(args[0].trim());
      } catch (NumberFormatException e) {
        // skip to the next line, this first one has no value
        line = br.readLine();
        // also get the number of parameters expected
        lineLength = line.trim().split("\\s+").length;
      }

      do {
        args = line.trim().split(",\\s+");
        // if a line doesn't match others in the file for length, or doesn't have at least
        // frequency and a mean value
        if (args.length != lineLength || args.length < 2) {
          throw new NoiseModelFormatException("Incorrect number of parameters per line in file "
              + file.getName());
        }
        // expect a common format for all data -- period, mean, median, 10th pct, 90th pct
        for (int i = 0; i < args.length; ++i) {
          try{
            double value = Double.parseDouble(args[i].trim());
            // they say the for-case is bad, but we're doing some very basic parsing so it's fine
            // after all what is parsing but an implicit state machine anyway
            switch(i) {
              case 0:
                periods.add(value);
                break;
              case 1:
                means.add(value);
                break;
              case 2:
                medians.add(value);
                break;
              case 3:
                tenths.add(value);
                break;
              case 4:
                ninetieths.add(value);
                break;
              default:
                break;
            }
          } catch (NumberFormatException e) {
            logger.error("Values in file " + file.getName() + " could not be parsed correctly,"
                + " specifically this value: " + args[i]);
          }
        }
      } while ((line = br.readLine()) != null);
    }
    int length = periods.size(); // this should be equal for all the arrays in question
    assert(length == means.size());
    this.periods = new double[length];
    this.meanLevels = new double[length];
    // these values are optional so we do no checks for them here
    this.medianLevels = new double[medians.size()];
    this.percentile10 = new double[tenths.size()];
    this.percentile90 = new double[ninetieths.size()];
    for (int i = 0; i < length; ++i) {
      this.periods[i] = periods.get(i);
      this.meanLevels[i] = means.get(i);
      if (i < medianLevels.length) {
        medianLevels[i] = medians.get(i);
      }
      if (i < percentile10.length) {
        percentile10[i] = tenths.get(i);
      }
      if (i < percentile90.length) {
        percentile90[i] = ninetieths.get(i);
      }
    }
  }

  /**
   * Return the filename (this is used to distinguish the noise model in plots)
   * @return filename of noise model source
   */
  public String getName() {
    return name;
  }

  /**
   * Return the period values and means as nested arrays. First dimension index 0 is the period
   * array and index 1 is the mean values. Second dimension k gives the kth value in the list, so
   * [0][k] and [1][k] give the kth period and mean value pair.
   * @return Nested array of period and mean level values.
   */
  public double[][] getPeriodAndMean() {
    return new double[][]{periods, meanLevels};
  }

  /**
   * Noise model exception thrown when the lines don't have the data expected.
   */
  public static class NoiseModelFormatException extends Exception {

    public NoiseModelFormatException(String s) {
      super(s);
    }
  }
}
