package asl.sensor.utils;

import org.apache.commons.math3.complex.Complex;

/**
 * Class containing methods to serve as math functions, mainly for angle calcs.
 * @author akearns
 *
 */
public class NumericUtils {

  /**
   * Get the two-component arctan of a complex number. A simpler way of calling
   * arctan2 using the real and imaginary components of a complex number
   * (see Math.atan2 for more details). This is the angle of the complex number
   * along the positive real axis.
   * @param c Complex number to get atan value of
   * @return atan, between -pi and pi
   */
  public static double atanc(Complex c) {
    return Math.atan2( c.getReal(), c.getImaginary() );
  }
  
  /**
   * Given a plot of data within a 2 * Pi range, check that a point is
   * continuous with a previous value.
   * @param phi Angle to fit within range of previous value (radians)
   * @param prevPhi Angle to check discontinuity against (radians)
   * @return New angle, with distance < Pi rad from the previous value
   */
  public static double unwrap(double phi, double prevPhi) {
    // sets range to [0,TAU], this syntax used because Java mod is weird
    phi = ( (phi % TAU) + TAU) % TAU;
    
    while ( Math.abs(prevPhi - phi) > Math.PI ) {
      if (prevPhi < phi) {
        phi -= TAU;
      } else {
        phi += TAU;
      }
    }
    
    return phi;
  }
  
  /**
   * Given a list of doubles representing the curve of a function with output
   * in radians over a range of 2 * pi, create a new curve that removes any
   * discontinuities in the plot. The starting point will be set to be as close
   * to zero as possible.
   * @param angles Array of input angles to make continuous (radians)
   * @return New array of angles where each pair of continuous points is
   * within distance < Pi rad from the previous value
   */
  public static double[] unwrapList(double[] angles) {
    double[] out = new double[angles.length];
    double prevPhi = 0.;
    
    for (int i = 0; i < out.length; ++i) {
      out[i] = unwrap(angles[i], prevPhi);
      prevPhi = out[i];
    }
    
    return out;
  }

  /**
   * 2 * Pi, sometimes also referred to as Tau. 
   * The number of radians in a full circle.
   */
  public final static double TAU = Math.PI * 2; // radians in full circle
}
