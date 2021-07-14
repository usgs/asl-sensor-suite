package asl.sensor.experiment;

import static asl.utils.NumericUtils.demeanInPlace;
import static asl.utils.NumericUtils.weightedAverageSlopesInterp;
import static asl.utils.timeseries.TimeSeriesUtils.ONE_HZ_INTERVAL;
import static java.lang.Double.isFinite;
import static java.lang.Double.isNaN;

import asl.sensor.input.DataStore;
import asl.utils.FFTResult;
import asl.utils.response.ChannelMetadata;
import asl.utils.response.ChannelMetadata.ResponseStageException;
import asl.utils.timeseries.TimeSeriesUtils;
import java.util.ArrayList;
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
