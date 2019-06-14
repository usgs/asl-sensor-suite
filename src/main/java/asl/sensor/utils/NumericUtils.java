package asl.sensor.utils;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

/**
 * Class containing methods to serve as math functions, mainly for angle calcs
 * and operations (stats) on complex numbers
 *
 * @author akearns
 */
public class NumericUtils {

  /**
   * 2 * Pi, sometimes also referred to as Tau.
   * The number of radians in a full circle.
   */
  public final static double TAU = Math.PI * 2; // radians in full circle

  public static final ThreadLocal<DecimalFormat> DECIMAL_FORMAT =
      ThreadLocal.withInitial(() -> {
        DecimalFormat format = new DecimalFormat("#.###");
        setInfinityPrintable(format);
        return format;
      });

  /**
   * Numerical Recipes cubic spline interpolation (spline.c and splint.c)
   * Expects arrays with +1 offset: x[1,...,n], etc. - we will pass it arrays
   * with 0 offset: x[0,1,...,n] and ignore the first points.
   *
   * @param x x
   * @param y y
   * @param n n
   * @param yp1 yp1
   * @param ypn ypn
   * @param y2 y2
   */
  private static void spline(double[] x, double[] y, int n, double yp1, double ypn, double[] y2) {

    double p, qn, sig, un;
    p = qn = sig = un = 0;

    // u=vector(1,n-1);
    double[] u = new double[n];

    if (yp1 > 0.99e30)
      y2[1] = u[1] = 0.0;
    else {
      y2[1] = -0.5;
      u[1] = (3.0 / (x[2] - x[1])) * ((y[2] - y[1]) / (x[2] - x[1]) - yp1);
    }
    for (int i = 2; i <= n - 1; i++) {
      sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]);
      p = sig * y2[i - 1] + 2.0;
      y2[i] = (sig - 1.0) / p;
      u[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]) - (y[i] - y[i - 1]) / (x[i] - x[i - 1]);
      u[i] = (6.0 * u[i] / (x[i + 1] - x[i - 1]) - sig * u[i - 1]) / p;
    }
    if (ypn > 0.99e30)
      qn = un = 0.0;
    else {
      qn = 0.5;
      un = (3.0 / (x[n] - x[n - 1])) * (ypn - (y[n] - y[n - 1]) / (x[n] - x[n - 1]));
    }
    y2[n] = (un - qn * u[n - 1]) / (qn * y2[n - 1] + 1.0);
    for (int k = n - 1; k >= 1; k--)
      y2[k] = y2[k] * y2[k + 1] + u[k];

  }

  /**
   * Same as above (+1 offset arrays) & y=double[1] is used to pass out the
   * interpolated value (y=f(x)).
   *
   * @param xa xa
   * @param ya ya
   * @param y2a y2a
   * @param n n
   * @param x x
   * @param y y
   */
  private static void splint(double[] xa, double[] ya, double[] y2a, int n, double x, double[] y) {

    int klo, khi, k;
    double h, b, a;

    klo = 1;
    khi = n;
    while (khi - klo > 1) {
      k = (khi + klo) >> 1;
      if (xa[k] > x)
        khi = k;
      else
        klo = k;
    }
    h = xa[khi] - xa[klo];
    // if (h == 0.0) nrerror("Bad XA input to routine SPLINT");
    if (h == 0.0)
      System.out.format("Bad XA input to routine SPLINT\n");

    a = (xa[khi] - x) / h;
    b = (x - xa[klo]) / h;
    // *y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
    y[0] = a * ya[klo] + b * ya[khi] + ((a * a * a - a) * y2a[klo] + (b * b * b - b) * y2a[khi]) * (h * h) / 6.0;
  }

  /**
   * Interpolate using a cubic spline.
   *
   * Interpolate measured Y[X] to the Y[Z]<br>
   * We know Y[X] = Y at values of X <br>
   * We want Y[Z] = Y interpolated to values of Z<br>
   *
   * @param X Measured X series
   * @param Y Measured Y series
   * @param Z Desired X coordinates to interpolate
   * @return interpolated array of Y values corresponding to input Z
   */
  public static double[] interpolate(double[] X, double[] Y, double[] Z) {

    double[] interpolatedValues = new double[Z.length];

    int n = X.length;

    double[] tmpY = new double[n + 1];
    double[] tmpX = new double[n + 1];

    // Create offset (+1) arrays to use with Num Recipes interpolation
    // (spline.c)
    for (int i = 0; i < n; i++) {
      tmpY[i + 1] = Y[i];
      tmpX[i + 1] = X[i];
    }
    double[] y2 = new double[n + 1];
    spline(tmpX, tmpY, n, 0., 0., y2);

    double[] y = new double[1];

    for (int i = 0; i < Z.length; i++) {
      splint(tmpX, tmpY, y2, n, Z[i], y);
      interpolatedValues[i] = y[0];
    }
    return interpolatedValues;
  }

  /**
   * Get the two-component arctan of a complex number. A simpler way of calling
   * arctan2 using the real and imaginary components of a complex number
   * (see Math.atan2 for more details). This is the angle of the complex number
   * along the positive real axis.
   *
   * @param c Complex number to get atan value of
   * @return atan, between -pi and pi
   */
  public static double atanc(Complex c) {
      
    // Setting this to a very small number to avoid jumps in phase response (e.g. RandomizedExperiment)
    // when the complex number is small. 
    // However, atan2 is discontinuous at 0, so we do not want that to be a possible input
  
    final double CUTOFF = 1e-16;  

    if (c.abs() < CUTOFF) {
      return 0.;
    }

    return Math.atan2(c.getImaginary(), c.getReal());
  }

  /**
   * Sort a list of complex values so that pure real numbers appear first.
   * The real and complex numbers are then each sorted according to their real
   * value, and then by their imaginary value.
   * @param complexes A list of complex numbers to sort
   */
  public static void complexRealsFirstSorter(List<Complex> complexes) {
    complexes.sort(CpxRealComparator.instance);
  }

  /**
   * Sort a list of complex values according to magnitude (i.e., period value).
   * Pole and zero values are iterated through twi
   */
  public static void complexMagnitudeSorter(Complex[] complexes) {
    Arrays.sort(complexes, CpxMagnitudeComparator.instance);
  }

  /**
   * Get the mean of the PSD calculation within the specified range
   *
   * @param psd PSD calculation (frequency space / FFT)
   * @param lower Starting index of window
   * @param higher Ending index of window
   * @return Mean of datapoints in the specific range
   */
  public static double getFFTMean(FFTResult psd, int lower, int higher) {
    double result = 0;
    int range = higher - lower;

    Complex[] data = psd.getFFT();

    for (int i = lower; i <= higher; ++i) {
      if (data[i].abs() >= Double.POSITIVE_INFINITY) {
        continue;
      }
      result += data[i].abs();
    }

    return result / range; // since result is a double, no cast needed?

  }

  /**
   * PSD Standard deviation calculation given mean ratio and specified range
   * (Get the result for mean calculation using the mean function before
   * calling this calculation; does not call the function directly)
   *
   * @param fft1 first PSD (numerator) of relative gain standard deviation
   * @param fft2 second PSD (denominator)
   * @param meanRatio mean(fft1)/mean(fft2)
   * @param lower Starting index of window (should be same used in mean calc)
   * @param higher Ending index of window (should be same used in mean calc)
   * @return double corresponding to standard deviation ratio over the window
   */
  public static double getFFTSDev(FFTResult fft1, FFTResult fft2, double meanRatio,
      int lower, int higher) {
    Complex[] density1 = fft1.getFFT();
    Complex[] density2 = fft2.getFFT();
    double sigma = 0.;
    for (int i = lower; i <= higher; ++i) {
      double value1 = density1[i].abs();
      double value2 = density2[i].abs();

      if (value1 >= Double.POSITIVE_INFINITY ||
          value2 >= Double.POSITIVE_INFINITY) {
        continue;
      }

      sigma += Math.pow((value1 / value2) - meanRatio, 2);
    }

    return Math.sqrt(sigma);
  }

  /**
   * Compute the standard deviation of a list of complex numbers. This will compute the deviation
   * on the real and imaginary parts separately, so that an error term can be generated for each.
   * @param cs List of complex numbers to do statistics over
   * @return Complex value representing standard deviation of real and imaginary part respectively
   */
  public static Complex getComplexSDev(Complex[] cs) {

    Complex summation = Complex.ZERO;
    for (Complex c : cs) {
      summation = summation.add(c);
    }
    summation = summation.divide(cs.length);

    double realSigma = 0.;
    double imaginarySigma = 0.;

    // now get the standard deviation of the real parts and the imaginary parts
    for (Complex c : cs) {
      double real = c.getReal();
      double imag = c.getImaginary();

      realSigma += Math.pow(real - summation.getReal(), 2);
      imaginarySigma += Math.pow(imag - summation.getImaginary(), 2);
    }

    return new Complex(Math.sqrt(realSigma / cs.length), Math.sqrt(imaginarySigma / cs.length));
  }

  /**
   * Perform a moving average on complex data, using the specified number of points to average at
   * each point in the input
   *
   * @param nums Complex data to be smoothed by use of moving average
   * @param points Number of points to include in moving average
   * @param forwardScan true if average should start from first point, false if starting from last
   * (this affects the shape of the curve slightly as the first few points are averaged with 0)
   * @return Smoothed data resulting from performing the moving average on input data.
   */
  public static Complex[] multipointMovingAverage(Complex[] nums, int points, boolean forwardScan) {
    if (points == 0) {
      return nums.clone();
    }
    DescriptiveStatistics realSide = new DescriptiveStatistics(points);
    DescriptiveStatistics imagSide = new DescriptiveStatistics(points);
    Complex[] out = new Complex[nums.length];
    for (int i = 0; i < nums.length; ++i) {
      int idx = i;
      if (!forwardScan) {
        idx = nums.length - (i + 1);
      }
      realSide.addValue(nums[idx].getReal());
      imagSide.addValue(nums[idx].getImaginary());
      Complex temp = new Complex(realSide.getMean(), imagSide.getMean());
      out[idx] = temp;
    }
    return out;
  }

  /**
   * Perform a moving average on real-val. data, using the specified number of points to average at
   * each point in the input
   *
   * @param nums Numeric data to be smoothed by use of moving average
   * @param points Number of points to include in moving average
   * @param forwardScan true if average should start from first point, false if starting from last
   * (this affects the shape of the curve slightly as the first few points are averaged with 0)
   * @return Smoothed data resulting from performing the moving average on input data.
   */
  public static double[] multipointMovingAverage(double[] nums, int points, boolean forwardScan) {
    DescriptiveStatistics windowStats = new DescriptiveStatistics(points);
    double[] out = new double[nums.length];

    for (int i = 0; i < nums.length; ++i) {
      int idx = i;
      if (!forwardScan) {
        idx = nums.length - (i + 1);
      }
      windowStats.addValue(nums[idx]);
      out[idx] = windowStats.getMean();
    }
    return out;
  }

  /**
   * Sets decimalformat object so that infinity can be printed in a PDF document
   *
   * @param df DecimalFormat object to change the infinity symbol value of
   */
  public static void setInfinityPrintable(DecimalFormat df) {
    DecimalFormatSymbols symbols = df.getDecimalFormatSymbols();
    symbols.setInfinity("Inf.");
    df.setDecimalFormatSymbols(symbols);
  }

  /**
   * Given a plot of data within a 2 * Pi range, check that a point is within
   * Pi of the previous value.
   *
   * @param phi Angle to fit within range of previous value (radians)
   * @param prevPhi Angle to check discontinuity against (radians)
   * @return New angle, with distance < Pi rad from the previous value
   */
  public static double unwrap(double phi, double prevPhi) {
    // sets range to [0,TAU]
    double newPhi = ((phi % TAU) + TAU) % TAU;

    while (Math.abs(prevPhi - newPhi) > Math.PI) {
      if (prevPhi < newPhi) {
        newPhi -= TAU;
      } else {
        newPhi += TAU;
      }
    }
    return newPhi;
  }

  /**
   * Given a list of doubles representing the curve of a function with output
   * in radians over a range of 2 * pi, create a new curve that removes any
   * discontinuities in the plot. The starting point will be set to be as close
   * to zero as possible.
   *
   * @param angles Array of input angles to make continuous (radians)
   * @return New array of angles where each pair of continuous points is
   * within distance < Pi rad from the previous value
   */
  public static double[] unwrapArray(double[] angles) {
    double[] out = new double[angles.length];
    double prevPhi = 0.;

    for (int i = 0; i < out.length; ++i) {
      out[i] = unwrap(angles[i], prevPhi);
      prevPhi = out[i];
    }

    return out;
  }

  /**
   * Wrap a degree to be between -180 and 180
   *
   * @param angle in degrees
   * @return same angle but between -180 and 180
   */
  public static double rewrapAngleDegrees(double angle) {
    while (angle < -180) {
      angle += 360;
    }
    while (angle > 180) {
      angle -= 360;
    }
    return angle;
  }

  /**
   * Complex comparator ordering by magnitude but listing pure-real values first
   * And then sorting complex numbers according to real value first, and then
   * by imaginary value if the real values match.
   *
   * @author akearns
   */
  static class CpxRealComparator implements Comparator<Complex> {

    static final CpxRealComparator instance = new CpxRealComparator();

    private CpxRealComparator() {

    }

    @Override
    public int compare(Complex c1, Complex c2) {

      if (null == c1 && null == c2) {
        return 0;
      } else if (null == c1) {
        return -1;
      } else if (null == c2) {
        return 1;
      }

      if (c1.getImaginary() == 0. && c2.getImaginary() == 0.) {
        return (int) Math.signum(c1.getReal() - c2.getReal());
      } else if (c1.getImaginary() == 0.) {
        return -1;
      } else if (c2.getImaginary() == 0.) {
        return 1;
      } else {
        if (c1.getReal() == c2.getReal()) {
          return (int) Math.signum(c1.getImaginary() - c2.getImaginary());
        } else {
          return (int) Math.signum(c1.getReal() - c2.getReal());
        }
      }
    }
  }

  /**
   * Complex comparator ordering by magnitude.
   *
   * @author akearns
   */
  static class CpxMagnitudeComparator implements Comparator<Complex> {

    static final CpxMagnitudeComparator instance = new CpxMagnitudeComparator();

    private CpxMagnitudeComparator() {

    }

    @Override
    public int compare(Complex c1, Complex c2) {

      if (null == c1 && null == c2) {
        return 0;
      } else if (null == c1) {
        return -1;
      } else if (null == c2) {
        return 1;
      }

      // assume
      if (c1.getReal() == c2.getReal()) {
        return (int) Math.signum(c1.getImaginary() - c2.getImaginary());
      }

      return (int) Math.signum(c1.abs() - c2.abs());
    }
  }
}


