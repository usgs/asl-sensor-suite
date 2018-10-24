package asl.sensor.experiment;

import asl.sensor.input.DataBlock;
import asl.sensor.input.DataStore;
import java.util.ArrayList;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

/**
 * Performs the calculations of a 10 volt test
 */
public class VoltageExperiment extends Experiment {

  private static final int PANELS_MAX = 3; // allow up to 3 inputs of data

  //private double sensitivity, gain;
  private double[] sensitivity;
  private double[] gain;
  private String[] dataNames;
  private int loadedAmount;

  public VoltageExperiment() {
    sensitivity = new double[]{};
    gain = new double[]{};
    dataNames = new String[]{};
    loadedAmount = 0;
  }

  @Override
  public String[] getDataStrings() {
    StringBuilder returnValue = new StringBuilder();
    for (int i = 0; i < loadedAmount; ++i) {
      returnValue.append(dataNames[i]).append(":")
          .append("\nRESP. Stage 2 gain: ").append(DECIMAL_FORMAT.get().format(gain[i]))
          .append("\nEstimated sensitivity: ").append(DECIMAL_FORMAT.get().format(sensitivity[i]))
          .append("\nPercent difference: ")
          .append(DECIMAL_FORMAT.get().format(getPercentDifference(i))).append("\n\n");
    }
    // remove the last newline character
    returnValue.delete(returnValue.length() - 2, returnValue.length() - 1);
    return new String[]{returnValue.toString()};
  }

  @Override
  protected void backend(DataStore dataStore) {

    // array here used to list what data store indices have data
    // that is actually loaded -- i.e., [0], [0, 2], [1, 2], [0, 1, 2] etc.
    int[] loadedData = new int[PANELS_MAX];
    loadedAmount = 0;
    for (int i = 0; i < PANELS_MAX; ++i) {
      if (dataStore.bothComponentsSet(i)) {
        loadedData[loadedAmount] = i;
        ++loadedAmount;
      }
    }


    gain = new double[loadedAmount];
    sensitivity = new double[loadedAmount];
    dataNames = new String[loadedAmount];

    xySeriesData = new ArrayList<>();
    XYSeriesCollection xysc = new XYSeriesCollection();

    for (int i = 0; i < loadedAmount; ++i) {
      int indexUnderAnalysis = loadedData[i];
      fireStateChange("Calculating input " + (i+1)
          + " of " + (loadedAmount) + " (slot " + (indexUnderAnalysis+1) + ")...");

      DataBlock currentBlock = dataStore.getBlock(indexUnderAnalysis);
      dataNames[i] = currentBlock.getName();

      // stage 2 gain is the value from the digitizer, i.e., 2^24/40 or 2^26/40
      // depending on the digitizer's bit-depth
      gain[i] = dataStore.getResponse(indexUnderAnalysis).getGain()[2];

      double[] data = currentBlock.getData();
      double sampleRate = currentBlock.getSampleRate();
      long startTime = currentBlock.getStartTime();
      double min = data[0];
      double max = data[0];
      int minIndex = 0; // track value indices to get ~1s range on each side
      int maxIndex = 0;

      for (int j = 3; j < data.length; ++j) {
        // make sure the extremes are in a roughly flat part of the signal
        // i.e., both min and max values should be in the flat part of a pulse
        double diff = Math.abs(data[j-3] - data[j]);
        if (data[j] < min && diff < 1000) {
          min = data[j];
          minIndex = j;
        }
        if (data[j] > max && diff < 1000) {
          max = data[j];
          maxIndex = j;
        }
      }

      // to reduce sucsceptibility to noise, get a few samples over
      double avgMin = 0.;
      double avgMax = 0.;
      XYSeries xys = new XYSeries(dataNames[i]);
      // get the 5 points centered around the max/min value
      for (int j = -2; j < 3; ++j) {
        int currentMinLookup = minIndex + j;
        int currentMaxLookup = maxIndex + j;
        long minTime = startTime + (long) sampleRate * currentMinLookup;
        long maxTime = startTime + (long) sampleRate * currentMaxLookup;
        xys.add(minTime, data[currentMinLookup]);
        xys.add(maxTime, data[currentMaxLookup]);
        avgMin += data[currentMinLookup];
        avgMax += data[maxIndex + i];
      }

      xysc.addSeries(xys);

      avgMin /= 5.;
      avgMax /= 5.;

      // take mean, divide by 10 (volts) -- average is the sensitivity (counts/volt)
      sensitivity[i] = (Math.abs(avgMin) + Math.abs(avgMax)) / 20.;
    }

    xySeriesData.add(xysc);

  }

  /**
   * Provide the estimated sensitivity taken from averaging the max and min data values
   * @param i Index of data to be gotten
   * @return Estimated sensitivity, expected to be close to
   */
  public double getSensitivity(int i) {
    return sensitivity[i];
  }

  /**
   * Provide the estimated sensitivity taken from averaging the max and min data values
   * for all inputs.
   * @return Array of estimated sensitivities, expected to be close to digitizer gain values
   * from RESPs
   */
  public double[] getAllSensitivities() {
    return sensitivity;
  }

  /**
   * Nominal value of digitizer gain, which sensitivity calculation will be compared against.
   *
   * This is specifically the stage 2 gain, usually (2^24)/40 or (2^26)/40 depending on whether
   * the sensor resolution is 24-bit or 26-bit respectively.
   * @param i Index of data to be gotten
   * @return Stage-2 gain taken from RESP file
   */
  public double getDigitizerGain(int i) {
    return gain[i];
  }

  /**
   * Nominal value of digitizer gain for each input, which sensitivity calculation will be compared
   * against.
   * This is specifically the stage 2 gain, usually (2^24)/40 or (2^26)/40 depending on whether
   * the sensor resolution is 24-bit or 26-bit respectively.
   * @return Array of Stage-2 gains taken from each input's RESP file
   */
  public double[] getAllGainValues() {
    return gain;
  }

  /**
   * Get the percent difference between nominal and calculated sensitivity.
   *
   * These values are ideally within +/- 0.5% of each other.
   * @param i Index of data to be gotten
   * @return Percent difference of sensitivity values
   */
  public double getPercentDifference(int i) {
    return (1. - sensitivity[i]/gain[i]) * 100;
  }

  /**
   * Get the percent difference between nominal and calculated sensitivity for all data loaded in.
   *
   * The percent difference is ideally within a range of +/- 0.5%.
   * @return Percent difference of sensitivity values
   */
  public double[] getPercentDifferences() {
    double[] differences = new double[loadedAmount];
    for (int i = 0; i < loadedAmount; ++i) {
      differences[i] = getPercentDifference(i);
    }
    return differences;
  }

  @Override
  public int blocksNeeded() {
    return 1;
  }

  @Override
  public boolean hasEnoughData(DataStore dataStore) {
    boolean somethingIsLoaded = false;
    for (int i = 0; i < PANELS_MAX; ++i) {
      somethingIsLoaded |= dataStore.bothComponentsSet(i);
    }
    return somethingIsLoaded;
  }


}
