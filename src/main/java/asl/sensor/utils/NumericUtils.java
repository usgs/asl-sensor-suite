package asl.sensor.utils;

import org.apache.commons.math3.complex.Complex;

public class NumericUtils {

  
  public static double atanc(Complex c) {
    return Math.atan2( c.getReal(), c.getImaginary() );
  }
  
  public static double unwrap(double phi, double prevPhi) {
    
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

  public final static double TAU = Math.PI * 2; // radians in full circle
}
