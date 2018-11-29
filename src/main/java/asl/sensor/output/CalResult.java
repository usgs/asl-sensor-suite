package asl.sensor.output;

import java.util.HashMap;
import java.util.Map;

/**
 * Abstract class that forms an easy interface by which external programs can interact with
 * all the data provided from a calibration solver experiment. CalResult contains two maps, one of
 * which is a map from string descriptors to a set of binary objects representing images as PNGs
 * (these are stored as byte arrays to be more easily imported into, say, a Django database backend)
 * and the other of which is a map from string descriptors to the variables fit by the solver,
 * given as a list of doubles (which has more than one entry in the case of, say, poles and zeros
 * returned by a randomized cal experiment).
 * Implementing classes don't need to add additional functions but must populate the maps with
 * actual data to be returned, which varies in content depending on the type of calibration done.
 * This class is not useful for the GUI interface, as the results there are contained within the
 * panel and saved to PDF reports as desired. In practice this is (currently) only used by the
 * ASL calibration tracking database.
 */
public class CalResult {

  /**
   * Get data from a random cal result
   * @param fitPoles fit poles from random cal data
   * @param fitZeros fit zeros from random cal data
   * @param initialPoles initial poles from random cal data
   * @param initialZeros initial zeros from random cal data
   * @param images plots converted to png-format images as byte arrays
   * @return object holding these values in easily-accessed maps with variable descriptions
   */
  public static CalResult buildRandomCalData(double[] fitPoles, double[] fitZeros,
      double[] initialPoles, double[] initialZeros, byte[][] images) {
    CalResult out = new CalResult();
    out.numerMap.put("Best_fit_poles", fitPoles);
    out.numerMap.put("Best_fit_zeros", fitZeros);
    out.numerMap.put("Initial_poles", initialPoles);
    out.numerMap.put("Initial_zeros", initialZeros);
    out.imageMap.put("Response_amplitudes", images[0]);
    out.imageMap.put("Response_phases", images[1]);
    out.imageMap.put("Amplitude_error", images[2]);
    out.imageMap.put("Phase_error", images[3]);
    return out;
  }

  /**
   * Get data from a step cal result
   * @param images plots converted to png-format images as byte arrays
   * @param initParams all initial parameters (corner, damping, residual)
   * @param fitParams all fit parameters
   * @return object holding these values in easily-accessed maps with variable descriptions
   */
  public static CalResult buildStepCalData(byte[][] images, double[] initParams,
      double[] fitParams) {
    double fitCorner = fitParams[0];
    double fitDamping = fitParams[1];
    double fitResid = fitParams[2];
    double initCorner = initParams[0];
    double initDamping = initParams[1];
    double initResid = initParams[2];
    CalResult out = new CalResult();
    out.numerMap.put("Fit_corner", new double[]{fitCorner});
    out.numerMap.put("Fit_damping", new double[]{fitDamping});
    out.numerMap.put("Fit_residual", new double[]{fitResid});
    out.numerMap.put("Initial_corner", new double[]{initCorner});
    out.numerMap.put("Initial_damping", new double[]{initDamping});
    out.numerMap.put("Initial_residual", new double[]{initResid});
    out.imageMap.put("Step_plot", images[0]);
    out.imageMap.put("Response_amplitudes", images[1]);
    out.imageMap.put("Response_phases", images[2]);
    return out;
  }

  /**
   * Get data from a sine cal result
   * @param images plots converted to png-format images as byte arrays
   * @param calAmp amplitude estimation of calibration signal
   * @param outAmp amplitude estimation of output signal
   * @param freq estimated frequency of sine wave
   * @param ratio ratio of signal amplitudes
   * @return object holding these values in easily-accessed maps with variable descriptions
   */
  public static CalResult buildSineCalData(byte[][] images, double calAmp, double outAmp,
      double freq, double ratio) {
    CalResult out = new CalResult();
    out.numerMap.put("Calibration_amplitude", new double[]{calAmp});
    out.numerMap.put("Output_signal_amplitude", new double[]{outAmp});
    out.numerMap.put("Estimated_signal_frequency", new double[]{freq});
    out.numerMap.put("Calibration_to_output_ratio", new double[]{ratio});
    out.imageMap.put("Sine_curves_plot", images[0]);
    out.imageMap.put("Linearity", images[1]);
    return out;
  }

  /**
   * Get data from a six-input gain result
   * @param firstAngleRef boolean that is true if first series being plotted is also angle ref
   * @param northAzimuth Angle over which north-facing data was rotated
   * @param eastAzimuth Angle over which east-facing data was rotated
   * @param statistics Result taken from six-input gain's
   * {@link asl.sensor.experiment.GainSixExperiment#getStatistics() getStatistics} function
   * @param images Plots of PSDs in each orientation
   * @return object holding these values in easily-accessed maps with variable descriptions
   */
  public static CalResult buildSixGainData(boolean firstAngleRef, double northAzimuth,
      double eastAzimuth, double[][] statistics, byte[][] images) {
    CalResult out = new CalResult();
    String rotatedSeries = "one";
    if (firstAngleRef) {
      rotatedSeries = "two";
    }
    out.numerMap.put("North_azimuth_of_series_" + rotatedSeries, new double[]{northAzimuth});
    out.numerMap.put("East_azimuth_of_series_" + rotatedSeries, new double[]{eastAzimuth});
    String[] orientations = {"North", "East", "Vertical"};
    // quick refresher: here is the order of fields from getStats
    // format: {mean, standard deviation, ref. gain, calc. gain, ref. A0 freq., calc. A0 freq.}
    for (int i = 0; i < orientations.length; ++i) {
      double[] orientationStats = statistics[i];
      double mean = orientationStats[0];
      out.numerMap.put(orientations[i] + "_mean_value", new double[]{mean});
      double standardDeviation = orientationStats[1];
      out.numerMap.put(orientations[i] + "_standard_deviation", new double[]{standardDeviation});
      double refGain = orientationStats[2];
      out.numerMap.put(orientations[i] + "_reference_gain", new double[]{refGain});
      double calcGain = orientationStats[3];
      out.numerMap.put(orientations[i] + "_calculated_gain", new double[]{calcGain});
      double refA0 = orientationStats[4];
      out.numerMap.put(orientations[i] + "_reference_a0_frequency", new double[]{refA0});
      double calcA0 = orientationStats[5];
      out.numerMap.put(orientations[i] + "_calculated_a0_frequency", new double[]{calcA0});
    }
    out.imageMap.put("North_plot", images[0]);
    out.imageMap.put("East_plot", images[1]);
    out.imageMap.put("Vertical_plot", images[2]);
    return out;
  }

  Map<String, double[]> numerMap;
  Map<String, byte[]> imageMap;

  private CalResult() {
    numerMap = new HashMap<>();
    imageMap = new HashMap<>();
  }

  /**
   * Get 10-volt test data
   * @param chart Chart including plotted sample data used to get estimations
   * @param gains Gain values from inputted data responses
   * @param sensitivities Estimated sensitivities as given by plotted data
   * @param differences Percent error of each pair of inputted gain and estimated sensitivity
   * @return object holding these values in easily-accessed maps with variable descriptions
   */
  public static CalResult buildVoltageData(byte[] chart, double[] gains, double[] sensitivities, double[] differences) {
    CalResult out = new CalResult();
    out.imageMap.put("Plot", chart);
    out.numerMap.put("Gain_values", gains);
    out.numerMap.put("Sensitivity_values", sensitivities);
    out.numerMap.put("Percent_differences", differences);

    return out;
  }

  /**
   * Return the map of images
   * @return map of byte arrays representing images, keyed by strings with image descriptions
   */
  public Map<String, byte[]> getImageMap() {
    return imageMap;
  }

  /**
   * Return the map of numeric data
   * @return map of double arrays representing calculation parameters, keyed by strings with
   * descriptions of the given numbers (i.e., fit poles, estimated sine wave frequency)
   */
  public Map<String, double[]> getNumerMap() {
    return numerMap;
  }
}
