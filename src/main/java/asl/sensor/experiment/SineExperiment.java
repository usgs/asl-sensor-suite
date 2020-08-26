package asl.sensor.experiment;

import static asl.utils.NumericUtils.demean;

import asl.sensor.input.DataStore;
import java.util.ArrayList;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

/**
 * This process compares the behavior for an input and output sine wave for a sensor. Unlike with
 * sine cals and step cals, a sine wave has a single (fixed) frequency, and thus no response is
 * used in this calculation. We take the RMS of the signals to get a ratio of their amplitudes,
 * and also calculate the distance between peaks of the sine wave to get an estimation of its
 * frequency. This estimation is likely to differ by a few Hz from the expected value in most cases,
 * and is mainly used as a quick way of checking that the calibration completed without significant
 * error making its results suspect (i.e., sine cals should have frequency in the sensor passband).
 * The main purpose of doing this cal is to get a midband sensitivity estimate for a seismometer,
 * under the assumption that the calibration components (i.e., coil abd drive) have fixed behavior.
 * For more details on the algorithm, see Ringler, Hutt, et al.,
 * "Obtaining Changes in Calibration-Coil to Seismometer Output Constants Using Sine Waves",
 * Bulletin of the Seismological Society of America, Vol 104 (Feb. 2014).
 *
 * @author akearns - KBRWyle
 */
public class SineExperiment extends Experiment {

  private double calSDev, outSDev, peakPeakFreq;

  public SineExperiment() {
    super();
    calSDev = 0.;
    outSDev = 0.;
    peakPeakFreq = 0.;
  }

  private String getResultData() {
    double calAmplitude = getCalAmplitude();
    double outAmplitude = getOutAmplitude();
    String calAmp = DECIMAL_FORMAT.get().format(calAmplitude);
    String outAmp = DECIMAL_FORMAT.get().format(outAmplitude);
    String ratio = DECIMAL_FORMAT.get().format(calAmplitude / outAmplitude);
    String estimatedFrequency = DECIMAL_FORMAT.get().format(getEstSineFreq());
    return "Calculated calibration amplitude: "
        + calAmp
        + "\nCalculated output amplitude: "
        + outAmp
        + "\nAmplitude ratio: "
        + ratio
        + "\nEstimated sine frequency: "
        + estimatedFrequency;
  }

  @Override
  public String[] getDataStrings() {
    return new String[]{getResultData()};
  }

  public double getCalAmplitude() {
    return calSDev;
  }

  public double getOutAmplitude() {
    return outSDev;
  }

  public double getEstSineFreq() {
    return peakPeakFreq;
  }

  @Override
  protected void backend(DataStore dataStore) {
    calSDev = 0.;
    outSDev = 0.;
    peakPeakFreq = 0.;

    double[] calTimeSeries = dataStore.getBlock(0).getData().clone();
    double[] outTimeSeries = dataStore.getBlock(1).getData().clone();

    dataNames.add(dataStore.getBlock(0).getName());
    dataNames.add(dataStore.getBlock(1).getName());

    // get the sine wave frequency by measuring wavelengths (get distance between peaks)
    int currentPeakDistance = 0;
    int totalPeakDistance = 0;
    int peakCount = -1; // start at negative 1 so first peak counts as 0
    for (int i = 1; i < outTimeSeries.length - 1; ++i) {
      if (outTimeSeries[i] > outTimeSeries[i - 1] && outTimeSeries[i] > outTimeSeries[i + 1]) {
        ++peakCount;
        totalPeakDistance += currentPeakDistance;
        currentPeakDistance = 0;
      } else {
        ++currentPeakDistance;
      }
    }
    if (peakCount > 0) {
      peakPeakFreq = totalPeakDistance / (double) peakCount;
    } else {
      peakPeakFreq = 0;
    }

    // we would like to have the signals centered on 0 for lining them up
    calTimeSeries = demean(calTimeSeries);
    outTimeSeries = demean(outTimeSeries);
    // standard deviation is a good
    calSDev = new DescriptiveStatistics(calTimeSeries).getStandardDeviation();
    outSDev = new DescriptiveStatistics(outTimeSeries).getStandardDeviation();

    // add plots sine waves
    XYSeriesCollection xysc = new XYSeriesCollection();
    XYSeries cal = new XYSeries(dataStore.getBlock(0).getName() + " [cal]");
    XYSeries out = new XYSeries(dataStore.getBlock(1).getName() + " [out, scaled]");
    double interval =
        dataStore.getBlock(0).getInterval();
    double start = getStart();
    for (int i = 0; i < calTimeSeries.length; ++i) {
      cal.add(start + interval * i, calTimeSeries[i]);
      out.add(start + interval * i, outTimeSeries[i] * calSDev / outSDev);
    }
    xysc.addSeries(cal);
    xysc.addSeries(out);
    xySeriesData = new ArrayList<>();
    xySeriesData.add(xysc);
    // produce linearity plots
    xysc = new XYSeriesCollection();
    // booleans: don't autosort, do allow duplicate values
    XYSeries lin = new XYSeries(dataStore.getBlock(1).getName() + " linearity", false, true);
    for (int i = 0; i < calTimeSeries.length; ++i) {
      lin.add(calTimeSeries[i], outTimeSeries[i]);
    }
    xysc.addSeries(lin);
    xySeriesData.add(xysc);
  }

  @Override
  public int blocksNeeded() {
    return 2;
  }

  @Override
  public boolean hasEnoughData(DataStore dataStore) {
    for (int i = 0; i < blocksNeeded(); ++i) {
      if (!dataStore.blockIsSet(i)) {
        return false;
      }
    }
    return true;
  }


}
