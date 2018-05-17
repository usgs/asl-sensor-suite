package asl.sensor.utils;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.List;
import org.apache.commons.math3.complex.Complex;
import org.junit.Test;

public class NumericUtilsTest {

  @Test
  public void movingAverageCorrectValues() {
    double[] averaged = new double[]{1., 1.5, 2.0, 3.0, 4.0};
    double[] init = new double[]{1., 2., 3., 4., 5.};
    double[] test = NumericUtils.multipointMovingAverage(init, 3, true);
    for (int i = 0; i < test.length; ++i) {
      assertEquals(averaged[i], test[i], 1E-25);
    }
  }

  @Test
  public void movingAverageCorrectValues_forwardScanFalse() {
    double[] averaged = new double[]{2.0, 3.0, 4.0, 4.5, 5.0};
    double[] init = new double[]{1., 2., 3., 4., 5.};
    double[] test = NumericUtils.multipointMovingAverage(init, 3, false);
    for (int i = 0; i < test.length; ++i) {
      assertEquals(averaged[i], test[i], 1E-25);
    }
  }

  @Test
  public void realSorterWorksRight() {
    List<Complex> cs = new ArrayList<>();
    cs.add(new Complex(1, 0));
    cs.add(new Complex(2, 0));
    cs.add(new Complex(-1, 0));
    cs.add(new Complex(-1, -1));
    cs.add(new Complex(1, -1));
    cs.add(new Complex(-3, -3));
    cs.add(new Complex(-3, 3));
    NumericUtils.complexRealsFirstSorter(cs);
    Complex[] csl = new Complex[cs.size()];
    csl[0] = new Complex(-1, 0);
    csl[1] = new Complex(1, 0);
    csl[2] = new Complex(2, 0);
    csl[3] = new Complex(-3, -3);
    csl[4] = new Complex(-3, 3);
    csl[5] = new Complex(-1, -1);
    csl[6] = new Complex(1, -1);
    for (int i = 0; i < csl.length; ++i) {
      assertTrue(cs.get(i).equals(csl[i]));
    }
  }

  @Test
  public void unwrap_testRangingTo0andTau_ChangesValues() {
    assertEquals(0, NumericUtils.unwrap(0 - NumericUtils.TAU, 0), 1E-12);
    assertEquals(0, NumericUtils.unwrap(NumericUtils.TAU + NumericUtils.TAU, 0), 1E-12);
    assertEquals(1, NumericUtils.unwrap(1 + NumericUtils.TAU, 0), 1E-12);
  }

  @Test
  public void unwrap_testRangingTo0andTau_UnchangedValues() {
    assertEquals(0, NumericUtils.unwrap(0, 0), 1E-12);
    assertEquals(0, NumericUtils.unwrap(NumericUtils.TAU, 0), 1E-12);
    assertEquals(1, NumericUtils.unwrap(1, 0), 1E-12);
  }

  @Test
  public void unwrap_prev0_nextCloseToTau() {
    assertEquals(NumericUtils.TAU * -0.01, NumericUtils.unwrap(NumericUtils.TAU * 0.99, 0), 1E-12);
    assertEquals(NumericUtils.TAU * 0.01, NumericUtils.unwrap(NumericUtils.TAU * 1.01, 0), 1E-12);
  }

  @Test
  public void unwrap_prevTau_nextCloseToTau() {
    assertEquals(NumericUtils.TAU * 0.99,
        NumericUtils.unwrap(NumericUtils.TAU * 0.99, NumericUtils.TAU), 1E-12);
    assertEquals(NumericUtils.TAU * 1.01,
        NumericUtils.unwrap(NumericUtils.TAU * 1.01, NumericUtils.TAU), 1E-12);
  }

  @Test
  public void unwrap_prevPi_nextCloseToTau() {
    assertEquals(NumericUtils.TAU * 0.99,
        NumericUtils.unwrap(NumericUtils.TAU * 0.99, NumericUtils.TAU / 2), 1E-12);
    assertEquals(NumericUtils.TAU * 0.01,
        NumericUtils.unwrap(NumericUtils.TAU * 1.01, NumericUtils.TAU / 2), 1E-12);
  }

  @Test
  public void unwrapList_StartsCloseToZero() {
    double[] input = {0 - NumericUtils.TAU, NumericUtils.TAU + NumericUtils.TAU,
        1 + NumericUtils.TAU, NumericUtils.TAU * 0.99};
    double[] expected = {0, 0, 1, NumericUtils.TAU * -0.01};
    assertArrayEquals(expected, NumericUtils.unwrapArray(input), 1E-12);
  }

  @Test
  public void rewrapAngleDegrees_unchanged() {
    assertEquals(179, NumericUtils.rewrapAngleDegrees(179), 1E-7);
    assertEquals(-179, NumericUtils.rewrapAngleDegrees(-179), 1E-7);
    assertEquals(0, NumericUtils.rewrapAngleDegrees(0), 1E-7);
    assertEquals(45, NumericUtils.rewrapAngleDegrees(45), 1E-7);
    assertEquals(180, NumericUtils.rewrapAngleDegrees(180), 1E-7);
    assertEquals(-180, NumericUtils.rewrapAngleDegrees(-180), 1E-7);
  }

  @Test
  public void rewrapAngleDegrees_higherthan180() {
    assertEquals(167, NumericUtils.rewrapAngleDegrees(527), 1E-7);
    assertEquals(-43, NumericUtils.rewrapAngleDegrees(677), 1E-7);
    assertEquals(0, NumericUtils.rewrapAngleDegrees(360), 1E-7);
    assertEquals(45, NumericUtils.rewrapAngleDegrees(405), 1E-7);
  }

  @Test
  public void rewrapAngleDegrees_lowerthanMinus180() {
    assertEquals(167, NumericUtils.rewrapAngleDegrees(-193), 1E-7);
    assertEquals(-43, NumericUtils.rewrapAngleDegrees(-763), 1E-7);
    assertEquals(0, NumericUtils.rewrapAngleDegrees(-360), 1E-7);
    assertEquals(45, NumericUtils.rewrapAngleDegrees(-1035), 1E-7);
  }
}

