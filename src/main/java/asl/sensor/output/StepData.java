package asl.sensor.output;

/**
 * Output format for step calibration solver as used in ASL's automated test tracking database
 */
public class StepData extends CalResult {

  // constuctor to be used with step calibrations
  public StepData(byte[][] images, double[] initParams, double[] fitParams) {
    super();
    double fitCorner = fitParams[0];
    double fitDamping = fitParams[1];
    double fitResid = fitParams[2];
    double initCorner = initParams[0];
    double initDamping = initParams[1];
    double initResid = initParams[2];
    numerMap.put("Fit_corner", new double[]{fitCorner});
    numerMap.put("Fit_damping", new double[]{fitDamping});
    numerMap.put("Fit_residual", new double[]{fitResid});
    numerMap.put("Initial_corner", new double[]{initCorner});
    numerMap.put("Initial_damping", new double[]{initDamping});
    numerMap.put("Initial_residual", new double[]{initResid});
    imageMap.put("Step_plot", images[0]);
    imageMap.put("Response_amplitudes", images[1]);
    imageMap.put("Response_phases", images[2]);
  }
}
