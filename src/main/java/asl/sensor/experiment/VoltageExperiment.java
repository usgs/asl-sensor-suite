package asl.sensor.experiment;

import asl.sensor.input.DataStore;
import java.util.ArrayList;
import java.util.List;
import org.apache.commons.math3.analysis.function.Exp;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

/**
 * Performs the calculations of a 10 volt test
 */
public class VoltageExperiment extends Experiment {

  private double sensitivity, gain;

  public VoltageExperiment() {
    sensitivity = 0.;
    gain = 0;
  }

  @Override
  public String[] getDataStrings() {
    String returnValue = "RESP. Stage 2 gain: " + DECIMAL_FORMAT.get().format(gain)
        + "\nEstimated sensitivity: " + DECIMAL_FORMAT.get().format(gain)
        + "\nPercent difference: " + DECIMAL_FORMAT.get().format(getPercentDifference());
    return new String[]{returnValue};
  }

  @Override
  protected void backend(DataStore dataStore) {

    List<XYSeriesCollection> xySeriesCollectionList = new ArrayList<>();


    fireStateChange("Calculating...");

    // stage 2 gain is the value from the digitizer, i.e., 2^24/40 or 2^26/40
    // depending on the digitizer's bit-depth
    gain = dataStore.getResponse(0).getGain()[2];

    double[] data = dataStore.getBlock(0).getData();
    double sampleRate = dataStore.getBlock(0).getSampleRate();
    long startTime = dataStore.getBlock(0).getStartTime();
    double min = data[0];
    double max = data[0];
    int minIndex = 0; // track value indices to get ~1s range on each side
    int maxIndex = 0;

    for (int i = 1; i < data.length; ++i) {
      if (data[i] < min) {
        min = data[i];
        minIndex = i;
      }
      if (data[i] > max) {
        max = data[i];
        maxIndex = i;
      }
    }

    // to reduce sucsceptibility to noise, get a few samples over
    double avgMin = 0.;
    double avgMax = 0.;
    XYSeries xys = new XYSeries(dataStore.getBlock(0).getName());
    // get the 5 points around the max/min value
    for (int i = -2; i < 3; ++i) {
      int currentMinLookup = minIndex + i;
      int currentMaxLookup = maxIndex + i;
      long minTime = startTime + (long) sampleRate * currentMinLookup;
      long maxTime = startTime + (long) sampleRate * currentMaxLookup;
      xys.add(minTime, data[currentMinLookup]);
      xys.add(maxTime, data[currentMaxLookup]);
      avgMin += data[currentMinLookup];
      avgMax += data[maxIndex + i];
    }

    XYSeriesCollection xysc = new XYSeriesCollection();
    xysc.addSeries(xys);

    xySeriesCollectionList.add(xysc);

    avgMin /= 5.;
    avgMax /= 5.;

    // take mean, divide by 10 (volts) -- average is the sensitivity (counts/volt)
    sensitivity = (Math.abs(avgMin) + Math.abs(avgMax)) / 20.;

    fireStateChange("Done!");
  }

  /**
   * Provide the estimated sensitivity taken from averaging the max and min data values
   * @return Estimated sensitivity, expected to be close to
   */
  public double getSensitivity() {
    return sensitivity;
  }

  /**
   * Nominal value of digitizer gain, which sensitivity calculation will be compared to.
   *
   * This is specifically the stage 2 gain, usually (2^24)/40 or (2^26)/40 depending on whether
   * the sensor resolution is 24-bit or 26-bit respectively.
   * @return Stage-2 gain taken from RESP file
   */
  public double getDigitizerGain() {
    return gain;
  }

  /**
   * Get the percent difference between nominal and calculated sensitivity.
   *
   * These values are ideally within +/- 0.5% of each other.
   * @return Percent difference of sensitivity values
   */
  public double getPercentDifference() {
    return (1. - sensitivity/gain) * 100;
  }

  @Override
  public int blocksNeeded() {
    return 1;
  }

  @Override
  public boolean hasEnoughData(DataStore dataStore) {
    return dataStore.bothComponentsSet(0);
  }


}
