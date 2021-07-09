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
import asl.utils.response.ChannelMetadata.ResponseStageException;
import asl.utils.timeseries.TimeSeriesUtils;
import java.util.ArrayList;
import org.apache.commons.math3.analysis.interpolation.HermiteInterpolator;
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

    fireStateChange("Removing responses from upsampled traces...");
    // next we get FFTs to remove the responses from each and convert back to the time domain
    try{
      testData = deconvolveResponse(testData, dataStore.getResponse(testIndex),
          dataStore.getBlock(testIndex).getSampleRate());
      refData = deconvolveResponse(refData, dataStore.getResponse(refIndex),
          dataStore.getBlock(refIndex).getSampleRate());
    } catch(ResponseStageException e) {
      System.err.println("Response is missing the expected pole-zero stage -- proceeding anyway");
    }

    // these values are what we will plot later, so we're keeping track of them separately
    double[] deconvolvedTest = testData;
    double[] deconvolvedRef = refData;

    fireStateChange("Interpolating data to 1 sample per ms...");
    // now, let's convert the samples to 1000 samples per second
    testData = weightedAverageSlopesInterp(testData, testInterval, THOUSAND_SPS_INTERVAL);
    refData = weightedAverageSlopesInterp(refData, refInterval, THOUSAND_SPS_INTERVAL);

    fireStateChange("Performing correlation calculation...");
    // assert(testData.length == refData.length);

    XYSeries correlationPlottable = new XYSeries("Correlation: " + testName + " & " + refName);
    // terms used for normalization
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

    // plot correlation curve normalized between 0 and 1
    for (int k = 0; k < correlations.length; ++k) {
      if (isFinite(correlations[k])) {
        double scaledValue = (correlations[k] - minValue) / (maxValue - minValue);
        correlationPlottable.add(k - centeringTerm, scaledValue);
      }
    }

    // now we have the lag time and the plottable correlation
    xySeriesData = new ArrayList<>();
    xySeriesData.add(new XYSeriesCollection(correlationPlottable));

    XYSeries referencePlot = new XYSeries(dataStore.getBlock(refIndex).getName());
    XYSeries shiftedTestPlot =
        new XYSeries(dataStore.getBlock(testIndex).getName() + "-shifted");
    // plot from the interpolated data -- now each index is 1ms, so shift i by interval
    for (int i = 0; i < deconvolvedRef.length; ++i) {
      referencePlot.add(getStart() + (i * refInterval), deconvolvedRef[i]);
    }
    for (int i = 0; i < deconvolvedTest.length; ++i) {
      shiftedTestPlot.add(getStart() + (i * testInterval) - lagTime, deconvolvedTest[i]);
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
    int length = refData.length + testData.length - 1;
    int shift = (2 * refData.length - 1) / 2;
    int pad = 2 * shift;
    testData = TimeSeriesUtils.concatAll(new double[pad], testData, new double[pad]);
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
        w[i] = 1. / Math.max(w[i], epsilon);
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
    hermiteInterpolation(oldArray, slope, interpolated, interval, newInterval);
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
    // double intervalDouble = (double) interval; // makes sure the division is floating-point
    double[] differentiated = new double[oldArray.length - 1];
    for (int i = 0; i < differentiated.length; ++i) {
      differentiated[i] = (oldArray[i + 1] - oldArray[i]) / interval;
    }
    return differentiated;
  }


  /**
   *
   * Perform Hermite spline interpolation on a set of data. Assume the same start and end times with
   * different sampling intervals between the two sets of data.
   * This routine is meant to replicate the obspy hermite interpolation routine (more of a
   * first-order spline), though it does not generate equivalent results.
   *
   * @param initial Initial data for interpolation over
   * @param slopes Calculated slopes of the
   * @param output Array to be populated with interpolated data
   * @param oldInterval Interval of original data
   * @param newInterval Interval of new data
   */
  private static void hermiteInterpolation(double[] initial, double[] slopes, double[] output,
      long oldInterval, long newInterval) {

    // note that because the two traces are expected to have the same start time, we choose our
    // time values to implicitly start at 0 and change based on the sample intervals
    int originalCurrentIndex = 0; // closest index before current time where we have an evaluation
    HermiteInterpolator interpolator = new HermiteInterpolator();
    // we will treat each pair of sampled points as defining an interval 0-1 over which we will
    // do interpolation on the points in between them (fractions between 0 and 1).
    // each sample point is treated as an array in and of itself, so first array is
    // the sampled value, followed by the first derivative

    // TODO: this is probably not the best way to structure the iteration
    output[0] = initial[0];
    interpolator.addSamplePoint(0, new double[]{initial[originalCurrentIndex]},
        new double[]{slopes[originalCurrentIndex]});
    interpolator.addSamplePoint(1, new double[]{initial[originalCurrentIndex + 1]},
        new double[]{slopes[originalCurrentIndex + 1]});


    for (int i = 1; i < output.length; i++) {
      // what is the current point in time given the index in the new data?
      long interpPointTime = i * newInterval;
      // if the closest index is new, then

      // what is the closest index to this in the original data?
      // if we're looking at a different sample, it's time to regenerate the interpolator
      if (originalCurrentIndex != (int) (interpPointTime / oldInterval)) {
        // we're at a new 0 point, which doesn't get interpolated as we have its value already
        originalCurrentIndex = (int) (interpPointTime / oldInterval);
        output[i] = initial[originalCurrentIndex];
        if (originalCurrentIndex + 1 == initial.length) {
          // if we hit this, we should have reached the end of the array already;
          assert(i + 1 == output.length);
          return;
        }
        interpolator = new HermiteInterpolator();
        interpolator.addSamplePoint(0, new double[]{initial[originalCurrentIndex]},
            new double[]{slopes[originalCurrentIndex]});
        interpolator.addSamplePoint(1, new double[]{initial[originalCurrentIndex + 1]},
            new double[]{slopes[originalCurrentIndex + 1]});
        continue;
      }

      // and what is the time at the expected index?
      long originalTime = originalCurrentIndex * oldInterval;

      // timeDifference is fractional relationship between points and old interval
      // it should always be between 0 and 1 -- otherwise the previous conditional should trigger
      double timeDifference = (interpPointTime - originalTime) / (double) oldInterval;
      assert(timeDifference > 0. && timeDifference < 1.);

      double term = interpolator.value(timeDifference)[0];
      output[i] = term;
    }
  }

  static double[] deconvolveResponse(double[] data, ChannelMetadata metadata,
      double sampleRate) throws ResponseStageException {
    demeanInPlace(data);
    FFTResult result1 = FFTResult.simpleFFT(data, sampleRate);
    Complex[] fft =  result1.getFFT();
    double[] frequencies = result1.getFreqs();
    // now get the responses
    Complex[] resp = metadata.applyResponseToInputUnscaled(frequencies, 1, 2);
    double normalization = metadata.getPoleZeroStage().getNormalizationFactor();

    // water level calculation -- since unscaled resp curve doesn't use A0, factor it in here
    resp[0] = resp[0].multiply(normalization);
    double max = resp[0].abs();
    for (int i = 1; i < fft.length; ++i) {
      resp[i] = resp[i].multiply(normalization);
      max = Math.max(resp[i].abs(), max);
    }
    double minLevel = max * Math.pow(10, (-0.001 / 20.0));
    for (int i = 0; i < fft.length; ++i) {
      // calculate water level correction and inversion of response
      if (resp[i].abs() > 0. && resp[i].abs() < minLevel) {
        resp[i] = resp[i].multiply(new Complex(minLevel / resp[i].abs()));
      }
      if (resp[i].abs() > 0.) {
        resp[i] = Complex.ONE.divide(resp[i]);
      }
      else if (resp[i].abs() == 0.) {
        resp[i] = Complex.ZERO;
      }
      fft[i] = fft[i].multiply(resp[i]);
    }
    data = FFTResult.simpleInverseFFT(fft, data.length);
    for (int i = 0; i < data.length; ++i) {
      data[i] /= metadata.getResponse().getInstrumentSensitivity().getSensitivityValue();
    }
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
