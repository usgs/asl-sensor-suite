package asl.sensor.experiment;

import asl.sensor.input.DataStore;
import asl.utils.input.DataBlock;
import java.util.ArrayList;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

/**
 * Performs the calculations of a 10 volt test to determine sensitivity values of a given channel.
 */
public class VoltageExperiment extends Experiment {

  private static final int PANELS_MAX = 3; // allow up to 3 inputs of data

  private double[] sensitivity;
  private double[] gain;
  private int[] gainStage;
  private int loadedAmount;

  public VoltageExperiment() {
    sensitivity = new double[]{};
    gain = new double[]{};
    // relevant gain stage is 3 in some newer installations, but default to 2
    gainStage = new int[]{2, 2, 2};
    loadedAmount = 0;
  }

  @Override
  public String[] getDataStrings() {
    String[] returnValue = new String[loadedAmount];
    for (int i = 0; i < loadedAmount; ++i) {

      returnValue[i] = dataNames.get(i)
          + ":\nRESP. Stage " + gainStage[i] + " gain: " + DECIMAL_FORMAT.get().format(gain[i])
          + "\nEstimated sensitivity: " + DECIMAL_FORMAT.get().format(sensitivity[i])
          + "\nPercent difference: "
          + DECIMAL_FORMAT.get().format(getPercentDifference(i));
    }
    return returnValue;
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
    // intialize as stage-2 gain
    gainStage = new int[]{2, 2, 2};

    xySeriesData = new ArrayList<>();
    XYSeriesCollection xysc = new XYSeriesCollection();

    for (int i = 0; i < loadedAmount; ++i) {
      int indexUnderAnalysis = loadedData[i];
      fireStateChange("Calculating input " + (i+1)
          + " of " + (loadedAmount) + " (slot " + (indexUnderAnalysis+1) + ")...");

      DataBlock currentBlock = dataStore.getBlock(indexUnderAnalysis);
      dataNames.add(currentBlock.getName());

      // stage 2 gain is the value from the digitizer, i.e., 2^24/40 or 2^26/40
      // depending on the digitizer's bit-depth
      gain[i] = dataStore.getResponse(indexUnderAnalysis).getGain()[2];
      if (gain[i] == 1) {
        // N4 sensors, possibly others may have a unit-value stage 2 gain
        // if that is the case when we need to look at the stage after that
        gain[i] = dataStore.getResponse(indexUnderAnalysis).getGain()[3];
        gainStage[i] = 3;
      }

      double[] data = currentBlock.getData();
      double min = data[0];
      double max = data[0];
      int minIndex = 0; // track value indices to get ~1s range on each side
      int maxIndex = 0;

      // offset by a second on either side to make sure that the data doesn't include ringing
      // artifacts -- i.e., data on either side of min/max should be flat relative to it
      int offset = (int) dataStore.getBlock(loadedData[i]).getSampleRate() + 1;

      outerLoop:
      for (int j = offset; j < data.length - offset; ++j) {

        // make sure the extremes are in a roughly flat part of the signal
        // i.e., all values to consider for analysis should have low variances
        for (int k = j - offset; k <= j + offset; ++k) {
          // represent error as magnitude difference between values
          double pctDiff = Math.abs(data[k] - data[j]) / Math.abs(data[j]) * 100;
          if (pctDiff > 5) {
            // if the percent error is too high, skip to next possible point
            continue outerLoop;
          }
        }

        if (data[j] <= min) {
          min = data[j];
          minIndex = j;
        }
        if (data[j] >= max) {
          max = data[j];
          maxIndex = j;
        }
      }

      // to reduce sucsceptibility to noise, get a few samples over
      double avgMin = 0.;
      double avgMax = 0.;
      int plotXPoint = 0;

      int startingPoint = -2; // start from 2 behind the minimum value of the data if possible
      while(minIndex + startingPoint < 0 || maxIndex + startingPoint < 0) {
        ++startingPoint;
      }

      XYSeries xys = new XYSeries(dataNames.get(i));
      // get the 5 points centered around the max/min value
      for (int j = startingPoint; j < (startingPoint + 5); ++j) {
        int currentMinLookup = minIndex + j;
        int currentMaxLookup = maxIndex + j;
        // x-value is just the given point in the set
        xys.add(++plotXPoint, Math.abs(data[currentMinLookup]));
        xys.add(plotXPoint + 5, Math.abs(data[currentMaxLookup]));
        avgMin += data[currentMinLookup];
        avgMax += data[maxIndex + i];
      }

      xysc.addSeries(xys);

      // take mean (div by 10 points), divide by 10 (volts) -- average is sensitivity (counts/volt)
      // we can divide by 2 to get mean because we took equal range of data from either side
      sensitivity[i] = (Math.abs(avgMin) + Math.abs(avgMax)) / (100.);
    }

    xySeriesData.add(xysc);

  }

  /**
   * Get an array representing the mean values of each trace's min and max values, for plotting for
   * each data loaded in.
   * @return array of values representing means of each set of inputs read in.
   */
  public double[] getMeanLines() {
    double[] meanLineValues = sensitivity.clone();
    for (int i = 0; i < meanLineValues.length; ++i) {
      meanLineValues[i] *= 10;
    }
    return meanLineValues;
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
