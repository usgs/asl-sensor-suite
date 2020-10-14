package asl.sensor.experiment;

import static asl.sensor.experiment.SineExperiment.*;
import static asl.utils.NumericUtils.demeanInPlace;
import static asl.utils.TimeSeriesUtils.ONE_HZ_INTERVAL;
import static java.lang.Math.sqrt;

import asl.sensor.input.DataStore;
import java.util.ArrayList;
import java.util.List;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

public class OrientedSineExperiment extends Experiment {

  /**
   * Cutoff value for amplitude between signals, which requires slighlty more
   */
  static final double AMP_CUTOFF = 5E-2;
  static final double FREQ_CUTOFF = 1E-2;
  static final int DATA_NEEDED = 3;
  private final double[] amplitude, frequency, phase;
  private double[] amplitudeError, frequencyError, phaseError;
  private double meanAmp, meanFreq, meanPhase;
  private boolean doRotation, isTrillium;

  private static final double[][] rotationMatrix = {
      {-sqrt(2. / 3.),        0.0, sqrt(1. / 3.)},
      { sqrt(1. / 6.),  sqrt(0.5), sqrt(1. / 3.)},
      { sqrt(1. / 6.), -sqrt(0.5), sqrt(1. / 3.)}
  };

  private static final double[][] rotationMatrixTrillium = {
      { sqrt(2. / 3.),        0.0, sqrt(1. / 3.)},
      {-sqrt(1. / 6.),  sqrt(0.5), sqrt(1. / 3.)},
      {-sqrt(1. / 6.), -sqrt(0.5), sqrt(1. / 3.)}
  };

  public OrientedSineExperiment() {
    super();
    amplitude = new double[DATA_NEEDED];
    frequency = new double[DATA_NEEDED];
    phase = new double[DATA_NEEDED];
    doRotation = false;
    isTrillium = false;
  }

  @Override
  protected void backend(DataStore dataStore) {
    // first, we need to get the components and (possibly) rotate them to UVW coordinates
    double[] north = dataStore.getBlock(0).getData();
    demeanInPlace(north);
    double[] east = dataStore.getBlock(1).getData();
    demeanInPlace(east);
    double[] vert = dataStore.getBlock(2).getData();
    demeanInPlace(vert);
    double testPoint = north[5];
    if (doRotation) {
      doRotation(north, east, vert); // rotation function modifies inputs in-place
    }
    dataNames.add(dataStore.getBlock(0).getName());
    dataNames.add(dataStore.getBlock(1).getName());
    dataNames.add(dataStore.getBlock(2).getName());
    double[][] traces = new double[][]{north, east, vert};
    xySeriesData = new ArrayList<>();
    // here we get the data into a plottable format
    {
      XYSeriesCollection xysc = new XYSeriesCollection();
      double interval = dataStore.getBlock(0).getInterval();
      double start = getStart();

      for (int j = 0; j < traces.length; ++j) {
        double[] trace = traces[j];
        XYSeries series = new XYSeries(dataNames.get(j));
        for (int i = 0; i < trace.length; ++i) {
          // TODO: add rotation steps here?
          series.add(start + interval * i, trace[i]);
        }
        xysc.addSeries(series);
      }
      xySeriesData.add(xysc);
    }
    // now that we have done that, we need to get the amplitude and frequency for each input
    Integer[][] allMaxPeakValues = new Integer[3][];
    for (int i = 0; i < DATA_NEEDED; ++i) {
      double[] trace = traces[i];
      Integer[][] values = getPeakLocations(trace);
      assert(values.length == 2);
      Integer[] maxes = values[0];
      Integer[] mins = values[1];
      allMaxPeakValues[i] = maxes;

      frequency[i] = calculateFrequency(trace);

      double amplitudes = 0;
      for (int j = 0; j < maxes.length; ++j) {
        if (j >= mins.length) break;
        // here we actually need the min and max peak values
        amplitudes += trace[maxes[j]] - trace[mins[j]];
      }
      amplitudes /= maxes.length;
      amplitude[i] = amplitudes;
    }

    // can't fit more than the fewest number of peaks
    int minLength =
        Math.min(Math.min(allMaxPeakValues[0].length, allMaxPeakValues[1].length),
                 allMaxPeakValues[2].length);
    double firstSecondDistance = 0;
    double firstThirdDistance = 0;
    double secondThirdDistance = 0;
    long interval = dataStore.getBlock(0).getInterval() / ONE_HZ_INTERVAL;
    for (int i = 0; i < minLength; ++i) {
      firstSecondDistance += Math.abs(allMaxPeakValues[0][i] - allMaxPeakValues[1][i]) * interval;
      firstThirdDistance += Math.abs(allMaxPeakValues[0][i] - allMaxPeakValues[2][i]) * interval;
      secondThirdDistance += Math.abs(allMaxPeakValues[1][i] - allMaxPeakValues[2][i]) * interval;
    }
    firstSecondDistance /= minLength;
    firstThirdDistance /= minLength;
    secondThirdDistance /= minLength;
    phase[0] = firstSecondDistance;
    phase[1] = firstThirdDistance;
    phase[2] = secondThirdDistance;

    // now that we have the values, do any of them have an error worse than our cutoff value?
    // (note: because CUTOFF is a decimal value and not percentage, don't need to scale by 100)
    // so we compare them to the average value for each and check the percent error
    meanAmp = 0;
    meanFreq = 0;
    meanPhase = 0;
    for (int i = 0; i < DATA_NEEDED; ++i) {
      meanAmp += (amplitude[i] / DATA_NEEDED);
      meanFreq += (frequency[i] / DATA_NEEDED);
      meanPhase += (phase[i] / DATA_NEEDED);
    }

    amplitudeError = checkError(amplitude, meanAmp);
    frequencyError = checkError(frequency, meanFreq);
    phaseError = checkError(phase, meanPhase);
    // values are all populated and the graphs were done before any of this processing
  }

  private void doRotation(double[] north, double[] east, double[] vert) {
    // performs rotation in-place
    double[][] rotation = isTrillium? rotationMatrixTrillium : rotationMatrix;
    for (int i = 0; i < north.length; ++i) {
      // since we do this in place, get the un-rotated XYZ values first
      double eastPoint = east[i];
      double northPoint = north[i];
      double vertPoint = vert[i];
      // this one is correct, so why don't the others work?
      double[] rotateOneAxis = rotation[0];
      east[i] = rotateOneAxis[0] * eastPoint +
          rotateOneAxis[1] * northPoint + rotateOneAxis[2] * vertPoint;

      rotateOneAxis = // new double[]{rotation[0][1], rotation[1][1], rotation[2][1]};
          rotation[1];
      north[i] = rotateOneAxis[0] * eastPoint +
          rotateOneAxis[1] * northPoint + rotateOneAxis[2] * vertPoint;

      rotateOneAxis = // new double[]{rotation[0][2], rotation[1][2], rotation[2][2]};
          rotation[2];
      vert[i] = rotateOneAxis[0] * eastPoint +
          rotateOneAxis[1] * northPoint + rotateOneAxis[2] * vertPoint;
    }
  }

  /*
  computes the percent error between the data and its mean value
   */
  private double[] checkError(double[] values, double mean) {
    double[] error = new double[values.length];
    for (int i = 0; i < error.length; ++i) {
      error[i] = Math.abs(values[i] - mean) / Math.abs(mean);
    }
    return error;
  }

  public boolean hasError() {
    for (int i = 0; i < amplitudeError.length; ++i) {
      if (amplitudeError[i] > AMP_CUTOFF || frequencyError[i] > FREQ_CUTOFF) {
        return true;
      }
    }
    return false;
  }

  public void setDoRotation(boolean doRotation) {
    this.doRotation = doRotation;
  }

  public void setRotationTrillium(boolean isTrillium) {
    this.isTrillium = isTrillium;
  }

  public double[] getAmplitudes() {
    return amplitude;
  }

  public double[] getFrequencies() {
    return frequency;
  }

  public double[] getPhaseDiscrepancies() {
    return phase;
  }

  public double[] getAmplitudeErrors() {
    return amplitudeError;
  }

  public double[] getFrequencyErrors() {
    return frequencyError;
  }

  public double[] getPhaseErrors() {
    return phaseError;
  }

  public double getMeanAmplitude() {
    return meanAmp;
  }

  public double getMeanFrequency() {
    return meanFreq;
  }

  public double getMeanPhase() {
    return meanPhase;
  }

  /*
   Produce the values of sine curve max values (frequency, phase, and amplitude calculations)
   */
  private Integer[][] getPeakLocations(double[] outTimeSeries) {
    List<Integer> maxPeakLocations = new ArrayList<>();
    List<Integer> minPeakLocations = new ArrayList<>();
    for (int i = 1; i < outTimeSeries.length - 1; ++i) {
      if (outTimeSeries[i] > outTimeSeries[i - 1] && outTimeSeries[i] > outTimeSeries[i + 1]) {
        maxPeakLocations.add(i);
      }
      if (outTimeSeries[i] < outTimeSeries[i - 1] && outTimeSeries[i] < outTimeSeries[i + 1]) {
        minPeakLocations.add(i);
      }
    }
    return new Integer[][]{maxPeakLocations.toArray(new Integer[]{}),
        minPeakLocations.toArray(new Integer[]{})};
  }

  @Override
  public int blocksNeeded() {
    return DATA_NEEDED;
  }

  @Override
  public boolean hasEnoughData(DataStore dataStore) {
    for (int i = 0; i < blocksNeeded(); ++i) {
      if (!dataStore.blockIsSet(i)) return false;
    }
    return true;
  }

  @Override
  public String[] getDataStrings() {
    // final String ERROR_STRING = "RESULTS OUTSIDE ERROR BOUND\n";
    final String ERROR_MARKER = " *** ";
    StringBuilder ampResults = new StringBuilder();
    ampResults.append("Mean amplitude: ").append(meanAmp).append('\n');
    StringBuilder freqResults = new StringBuilder();
    freqResults.append("Mean frequency: ").append(
        DECIMAL_FORMAT.get().format(meanFreq)).append('\n');
    StringBuilder phaseResults = new StringBuilder();
    phaseResults.append("Mean phase value ").append(
        DECIMAL_FORMAT.get().format(meanPhase)).append('\n');
    String[] ordinals = new String[]{"N", "E", "Z"};
    for (int i = 0; i < DATA_NEEDED; ++i) {

      ampResults.append(ordinals[i]).append(" input amp.: ")
          .append(DECIMAL_FORMAT.get().format(amplitude[i]));
      freqResults.append(ordinals[i]).append(" input freq.: ")
          .append(DECIMAL_FORMAT.get().format(frequency[i]));
      // and also include the error terms
      ampResults.append(" (err: ").append(
          DECIMAL_FORMAT.get().format(amplitudeError[i] * 100)).append(")\n");
      ampResults.append('\n');
      freqResults.append(" (err: ").append(
          DECIMAL_FORMAT.get().format(frequencyError[i] * 100)).append(")\n");
      freqResults.append('\n');
    }
    phaseResults.append("Phase diff. of first and second input: ")
        .append(DECIMAL_FORMAT.get().format(phase[0])).append('\n');
    phaseResults.append("Phase diff. of first and third input: ")
        .append(DECIMAL_FORMAT.get().format(phase[1])).append('\n');
    phaseResults.append("Phase diff. of second and third: ")
        .append(DECIMAL_FORMAT.get().format(phase[2])).append('\n');

    return new String[]{ampResults.toString(), freqResults.toString(), phaseResults.toString()};
  }
}
