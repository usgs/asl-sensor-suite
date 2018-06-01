package asl.sensor.output;

/**
 * Output format for sine calibration solver as used in ASL's automated test tracking database
 */
public class SineData extends CalResult {

  public SineData(byte[][] images, double calAmp, double outAmp, double freq, double ratio) {
    super();
    numerMap.put("Calibration_amplitude", new double[]{calAmp});
    numerMap.put("Output_signal_amplitude", new double[]{outAmp});
    numerMap.put("Estimated_signal_frequency", new double[]{freq});
    numerMap.put("Calibration_to_output_ratio", new double[]{ratio});
    imageMap.put("Sine_curves_plot", images[0]);
    imageMap.put("Linearity", images[1]);
  }
}