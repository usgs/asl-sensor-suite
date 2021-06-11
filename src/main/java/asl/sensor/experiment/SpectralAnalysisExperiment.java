package asl.sensor.experiment;

import asl.sensor.input.NoiseModel;
import asl.sensor.input.NoiseModel.NoiseModelFormatException;
import java.io.File;
import java.io.IOException;
import org.jfree.data.xy.XYSeries;

/**
 * Abstract class for behavior common to all experiments based on PSD calculation/analysis,
 * mainly for keeping track of station noise models
 */
public abstract class SpectralAnalysisExperiment extends Experiment {
  protected NoiseModel stationNoiseModel;

  /**
   * Insert noise model
   * @param filename Name of noise model file to load
   * @throws IOException if the file cannot be read or does not exist
   * @throws NoiseModelFormatException if the noise model does not fit the expected format
   */
  public void loadNoiseModel(String filename) throws IOException, NoiseModelFormatException {
    stationNoiseModel = new NoiseModel(filename);
  }

  public void loadNoiseModel(File file) throws IOException, NoiseModelFormatException {
    stationNoiseModel = new NoiseModel(file);
  }

  public void clearNoiseModel() {
    stationNoiseModel = null;
  }

  public NoiseModel getNoiseModel() {
    return stationNoiseModel;
  }

  public boolean noiseModelLoaded() {
    return stationNoiseModel != null;
  }

  /**
   * Return the noise model to compare with the results of this experiment, plotting
   * period (x-axis) vs. the mean value of the noise measured (y-axis).
   * This is not added to the primary plot data structure as part of backend calculations because
   * it may be kept null, and so that it can be added to the respected graph at any time
   * (i.e., not requiring a full replot of the experiment in order to perform it).
   * @param useHzUnits true if the data should be plotted with x-axis in Hz; if false, the data
   *                   will be plotted with x-axis representing sample interval in seconds
   * @return XYSeries (plottable) of the period and estimated noise level
   */
  public XYSeries getPlottableNoiseModelData(boolean useHzUnits) {
    if (stationNoiseModel == null) {
      return null;
    }

    XYSeries plot = new XYSeries("Model: " + stationNoiseModel.getName());
    double[][] periodAndMean = stationNoiseModel.getPeriodAndMean();
    int length = periodAndMean[0].length;
    assert(length == periodAndMean[1].length);
    for (int i = 0; i < length; ++i) {
      double x = useHzUnits ? 1. / periodAndMean[0][1] : periodAndMean[0][1];
      plot.add(x, periodAndMean[1][i]);
    }

    return plot;
  }
}
