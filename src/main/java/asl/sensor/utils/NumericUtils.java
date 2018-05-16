package asl.sensor.utils;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Collections;
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
   * Complex comparator ordering by magnitude but listing pure-real values first
   * And then sorting complex numbers according to real value first, and then
   * by imaginary value if the real values match.
   *
   * @author akearns
   */
  public static class CpxRealComparator implements Comparator<Complex> {

    public static final CpxRealComparator instance = new CpxRealComparator();

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
   * 2 * Pi, sometimes also referred to as Tau.
   * The number of radians in a full circle.
   */
  public final static double TAU = Math.PI * 2; // radians in full circle

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

    final double CUTOFF = 1. / 1000.;

    if (c.abs() < CUTOFF) {
      return 0.;
    }

    return Math.atan2(c.getImaginary(), c.getReal());
  }

  /**
   * Sort a list of complex values so that pure real numbers appear first.
   * The real and complex numbers are then each sorted according to their real
   * value, and then by their imaginary value.
   */
  public static void complexRealsFirstSorter(List<Complex> complexes) {
    complexes.sort(CpxRealComparator.instance);
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

}


