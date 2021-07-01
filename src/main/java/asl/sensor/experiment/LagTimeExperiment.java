package asl.sensor.experiment;

import static asl.utils.NumericUtils.demeanInPlace;
import static asl.utils.timeseries.TimeSeriesUtils.ONE_HZ_INTERVAL;
import static java.lang.Double.isFinite;
import static java.lang.Double.isInfinite;
import static java.lang.Double.isNaN;
import static java.lang.Math.signum;

import asl.sensor.input.DataStore;
import asl.utils.FFTResult;
import asl.utils.response.ChannelMetadata;
import asl.utils.timeseries.TimeSeriesUtils;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import org.apache.commons.math3.complex.Complex;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

/**
 * Estimates lag time between two signals by performing the 1D correlation of signals from two
 * separate traces. The traces have been trimmed to within an upper bound of 5000 points and
 * had their responses removed, then upsampled to 1 sample per ms. Lag is determined by the
 * index of max correlation -- this represents the timing difference in ms between the traces.
 *
 * The algorithm as presented here is based on Ringler, A. T.,
 * R. E. Anthony, D. C. Wilson, D. Auerbach,
 * S. Bargabus, P. Davis, M. Gunnels,
 * K. Hafner, J. F. Holland, A. Kearns, et al.
 * (2021). A Review of Timing Accuracy
 * across the Global Seismographic Network,
 * Seismol. Res. Lett. XX, 1-12, doi: 10.1785/0220200394
 */
public class LagTimeExperiment extends Experiment {

  // one thousand samples per second = one sample per millisecond
  private static final long THOUSAND_SPS_INTERVAL = ONE_HZ_INTERVAL / 1000;
  private static final int MAX_DATA_LIMIT = 5000; // sample limit for data

  private int lagTime; // shift index based on max value from correlation check
  // which is lag time in ms because samples are upscaled to 1000 sps

  public LagTimeExperiment() {
    lagTime = 0;
  }

  @Override
  public String[] getInsetStrings() {
    return new String[]{"Estimated lag time: " + getLagTime() + "ms"};
  }

  @Override
  protected boolean doMatchIntervals() {
    // we DO NOT match intervals in this method because we will upscale later anyway
    return false;
  }

  @Override
  protected void backend(DataStore dataStore) {
    final int testIndex = 0;
    final int refIndex = 1;
    double[] testData = dataStore.getBlock(testIndex).getData();
    double[] refData = dataStore.getBlock(refIndex).getData();

    String testName = dataStore.getBlock(testIndex).getName();
    String refName = dataStore.getBlock(refIndex).getName();

    dataNames.add(testName);
    dataNames.add(dataStore.getResponse(testIndex).getName());
    dataNames.add(refName);
    dataNames.add(dataStore.getResponse(refIndex).getName());

    // note that data is of course always coerced to match sample rates going in
    // (if we need to change this, we need to override the preprocessing routine too)
    long testInterval = dataStore.getBlock(testIndex).getInterval();
    long refInterval = dataStore.getBlock(refIndex).getInterval();

    try (PrintWriter out = new PrintWriter(new FileWriter("lagtime-before-response.csv"))) {
      StringBuilder outString = new StringBuilder(Arrays.toString(testData));
      out.write(testName + "\n");
      out.write(outString.substring(1, outString.length()-1));
      outString = new StringBuilder(Arrays.toString(refData));
      out.write("\n" + refName + "\n");
      out.write(outString.substring(1, outString.length()-1));
    } catch (IOException e) {
      e.printStackTrace();
    }

    fireStateChange("Removing responses from upsampled traces...");
    // next we get FFTs to remove the responses from each and convert back to the time domain
    {
      testData = deconvolveResponse(testData, dataStore.getResponse(testIndex));
      refData = deconvolveResponse(refData, dataStore.getResponse(refIndex));
    }

    try (PrintWriter out = new PrintWriter(new FileWriter("lagtime-after-response.csv"))) {
      StringBuilder outString = new StringBuilder(Arrays.toString(testData));
      out.write(testName + "\n");
      out.write(outString.substring(1, outString.length()-1));
      outString = new StringBuilder(Arrays.toString(refData));
      out.write("\n" + refName + "\n");
      out.write(outString.substring(1, outString.length()-1));
    } catch (IOException e) {
      e.printStackTrace();
    }

    fireStateChange("Interpolating data to 1 sample per ms...");
    // now, let's convert the samples to 1000 samples per second
    testData = weightedAverageSlopesInterp(testData, testInterval, THOUSAND_SPS_INTERVAL);
    refData = weightedAverageSlopesInterp(refData, refInterval, THOUSAND_SPS_INTERVAL);

    // now that we have upsampled the data we should go ahead and assert that the data is
    // of a matching length just in case
    {
      int trimLength = Math.min(testData.length, refData.length);
      // trimLength = Math.min(trimLength, MAX_DATA_LIMIT);
      testData = testData.length == trimLength ?
          testData : Arrays.copyOfRange(testData, 0, trimLength);
      demeanInPlace(testData);
      refData = refData.length == trimLength ?
          refData : Arrays.copyOfRange(refData, 0, trimLength);
      demeanInPlace(refData);
    }

    // demean in place ok here because doing FFT-invFFT creates new array of data anyway
    demeanInPlace(testData);
    demeanInPlace(refData);

    fireStateChange("Performing correlation calculation...");
    // now it is time to actually calculate the correlation, but first we must pad the first trace
    assert(testData.length == refData.length);

    XYSeries correlationPlottable = new XYSeries("Correlation: " + testName + " & " + refName);
    double maxValue = Double.NEGATIVE_INFINITY;
    double minValue = Double.POSITIVE_INFINITY;
    // get the index of max value, representing lag time in ms relative to midpoint of data
    lagTime = 0;
    double[] correlations = getCorrelation(refData, testData);
    int centeringTerm = (correlations.length / 2); // index of midpoint of data, point of 0ms lag
    for (int k = 0; k < correlations.length; ++k) {
      minValue = Math.min(minValue, correlations[k]);
      if (correlations[k] > maxValue) {
        // center the index around the middle value (where peak should be if lag is 0)
        lagTime = k - centeringTerm;
        maxValue = correlations[k];
      }
    }

    for (int k = 0; k < correlations.length; ++k) {
      if (isFinite(correlations[k])) {
        double scaledValue = (correlations[k] - minValue) / (maxValue - minValue);
        correlationPlottable.add(k - centeringTerm, scaledValue);
      }
    }

    // now we have the lag time and the plottable correlation
    xySeriesData = new ArrayList<>();
    xySeriesData.add(new XYSeriesCollection(correlationPlottable));

    // TODO: make sure this data gets actually plotted in the panel in its own chart
    XYSeries referencePlot = new XYSeries(dataStore.getBlock(refIndex).getName());
    XYSeries shiftedTestPlot =
        new XYSeries(dataStore.getBlock(testIndex).getName() + "-shifted");
    for (int i = 0; i < refData.length; ++i) {
      referencePlot.add(getStart() + (i * refInterval), refData[i]);
    }
    // plot from the interpolated data -- now each index is 1ms, so shift i by interval
    for (int i = 0; i < testData.length; i += testInterval) {
      // i is already a multiple of the interval, so we don't do a scaling step here
      shiftedTestPlot.add(getStart() + i + lagTime, testData[i]);
    }
    {
      XYSeriesCollection shiftedPlots = new XYSeriesCollection();
      shiftedPlots.addSeries(referencePlot);
      shiftedPlots.addSeries(shiftedTestPlot);
      xySeriesData.add(shiftedPlots);
    }
  }

  public int getLagTime() {
    return lagTime;
  }

  static double[] getCorrelation(double[] refData, double[] testData) {
    int shift = (2 * testData.length - 1) / 2;
    int pad = 2 * shift;
    testData = TimeSeriesUtils.concatAll(new double[pad], testData, new double[pad]);
    int length = pad + 1;
    double[] correlations = new double[length];
    for (int k = 0; k < correlations.length; ++k) {
      for (int n = 0; n < refData.length; ++n) {
        assert (!(isNaN(refData[n])));
        correlations[k] += testData[n + shift + k] * refData[n];
      }
    }
    return correlations;
  }

  // TODO: this stuff will eventually be moved as part of java utils

  /**
   * Performs interpolation meant to replicate the the obspy method weighted_average_slopes,
   * albeit without trimming the data down to a new time.
   * @param oldArray Array of data to interpolate
   * @param interval Interval of data (in ms)
   * @param newInterval Interval of new data (in ms)
   * @return Interpolated data over the same range of time as the original data
   */
  public static double[] weightedAverageSlopesInterp(double[] oldArray, long interval,
      long newInterval) {
    int interpolatedLength = (int) ((oldArray.length - 1) * interval / newInterval) + 1;
    double[] interpolated = new double[interpolatedLength];
    double[] slope = new double[oldArray.length]; // weighted derivative values
    // calculate the weighted slope array
    {
      double[] m = differentiate(oldArray, interval);
      double[] w = new double[m.length];
      final double epsilon = Math.ulp((float) 0.);
      for (int i = 0; i < m.length; ++i) {
        w[i] = Math.abs(m[i]);
        w[i] = Math.max(w[i], epsilon);
        w[i] = 1. / w[i];
        assert(!isInfinite(w[i]));
      }

      slope[0] = m[0];
      for (int i = 0; i < m.length - 1; ++i) {
        double currentW = w[i];
        double nextW = w[i + 1];
        double currentM = m[i];
        double nextM = m[i + 1];
        slope[i + 1] = ((currentW * currentM) + (nextW * nextM)) / (currentW + nextW);
        assert(!isNaN(slope[i]));
      }
      // last value of slope is last value of m
      slope[slope.length - 1] = m[m.length - 1];
      // now ensure that any sign changes become 0
      for (int i = 0; i < m.length - 1; ++i) {
        if (signum(m[i]) - signum(m[i + 1]) != 0) {
          slope[i + 1] = 0;
        }
      }
    }

    // now we do hermite interpolation to match obspy -- polynomial piecewise considered too
    // memory intensive for use in obspy
    double oldInterval = (double) interval / ONE_HZ_INTERVAL;
    double newIntervalDouble = (double) newInterval / ONE_HZ_INTERVAL;
    // the obspy hermite interpolation method is stolen here because the apache commons one
    // doesn't work the way we expect it to (throws a lot of NaNs).
    hermiteInterpolation(oldArray, slope, interpolated, oldInterval, newIntervalDouble);
    return interpolated;
  }

  /**
   * Differentiate data. Intended to match the numpy.diff method, producing an array of length
   * n-1 where n is the length of the original array.
   * @param oldArray Original set of data, assumed spaced on an equal interval
   * @param interval Spacing between values (i.e., if y=f(x) this value represents dx)
   * @return Differentiated data.
   */
  public static double[] differentiate(double[] oldArray, long interval) {
    double intervalDouble = (double) interval / ONE_HZ_INTERVAL;
    double[] differentiated = new double[oldArray.length - 1];
    for (int i = 0; i < differentiated.length; ++i) {
      differentiated[i] = (oldArray[i + 1] - oldArray[i]) / intervalDouble;
    }
    return differentiated;
  }


  /**
   *
   * Perform Hermite interpolation on a set of data. Assume the same start and end times with
   * different sampling intervals between the two sets of data.
   *
   * Based on C code used in obspy. Original code is (C) by Lion Krischer, 2014.
   * That original code, as part of obspy, is licensed under LGPL 3.0.
   * Original source is available here:
   * https://github.com/obspy/obspy/blob/9dd8ccc493ba6e0d70d26966246e5d7641c5a6d0/obspy/signal/src/hermite_interpolation.c
   *
   * @param initial Initial data for interpolation over
   * @param slopes Calculated slopes of the
   * @param output Array to be populated with interpolated data
   * @param oldInterval Interval of original data
   * @param newInterval Interval of new data
   */
  private static void hermiteInterpolation(double[] initial, double[] slopes, double[] output,
      double oldInterval, double newInterval) {

    double a_0, a_1, b_minus_1, b_plus_1, b_0, c_0, c_1, d_0;
    // note that because the two traces are expected to have the same start time, we choose our
    // time values to implicitly start at 0 and change based on the sample intervals
    for (int idx = 0; idx < output.length; idx++) {
      // what is the current point in time given the index in the new data?
      double interpPointTime = idx * newInterval;
      // what is the closest index to this in the original data?
      int originalCurrentIndex = (int) (interpPointTime / oldInterval);
      int originalNextIndex = originalCurrentIndex + 1;
      // and what is the time at the expected index?
      double originalTime = originalCurrentIndex * oldInterval;

      // No need to interpolate if exactly at start of the interval.
      if (interpPointTime == originalTime)  {
        output[idx] = initial[originalCurrentIndex];
        continue;
      }

      double timeDifference = interpPointTime - originalTime;

      // matrix manipulation, in effect
      a_0 = initial[originalCurrentIndex];
      a_1 = initial[originalNextIndex];
      b_minus_1 = oldInterval * slopes[originalCurrentIndex];
      b_plus_1 = oldInterval * slopes[originalNextIndex];
      b_0 = a_1 - a_0;
      c_0 = b_0 - b_minus_1;
      c_1 = b_plus_1 - b_0;
      d_0 = c_1 - c_0;

      output[idx] = a_0 + (b_0 + (c_0 + d_0 * timeDifference) * (timeDifference - 1.0)) * timeDifference;
    }
  }

  private static double[] deconvolveResponse(double[] data, ChannelMetadata metadata) {
    double sampleRate = (double) THOUSAND_SPS_INTERVAL / ONE_HZ_INTERVAL;
    FFTResult result1 = FFTResult.singleSidedFFT(data, sampleRate, false);
    Complex[] fft =  result1.getFFT();
    double[] frequencies = result1.getFreqs();
    // now get the responses
    Complex[] resp = metadata.applyResponseToInput(frequencies);
    int startIndex = 0;
    if (frequencies[0] == 0.) {
      fft[0] = Complex.ZERO;
      startIndex = 1;
    }
    for (int i = startIndex; i < fft.length; ++i) {
      fft[i] = fft[i].divide(resp[i]);
    }
    data = FFTResult.singleSidedInverseFFT(fft, data.length);
    assert(!isNaN(data[0]));
    return data;
  }

  @Override
  public int blocksNeeded() {
    return 2;
  }

  @Override
  public boolean hasEnoughData(DataStore dataStore) {
    for (int i = 0; i < blocksNeeded(); ++i) {
      if (!dataStore.bothComponentsSet(i)) {
        return false;
      }
    }
    return true;
  }
}
