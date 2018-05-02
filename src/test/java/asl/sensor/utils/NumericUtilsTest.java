package asl.sensor.utils;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import org.apache.commons.math3.complex.Complex;
import org.junit.Test;

public class NumericUtilsTest {

  @Test
  public void movingAverageCorrectValues() {
    double[] averaged = new double[]{1., 1.5, 2.0, 3.0, 4.0};
    double[] init = new double[]{1., 2., 3., 4., 5.};
    double[] test = NumericUtils.multipointMovingAverage(init, 3, true);
    System.out.println( Arrays.toString(test) );
    for (int i = 0; i < test.length; ++i) {
      assertEquals(averaged[i], test[i], 1E-25);
    }
  }

  @Test
  public void realSorterWorksRight() {
    List<Complex> cs = new ArrayList<>();
    cs.add( new Complex(1, 0) );
    cs.add( new Complex(2, 0) );
    cs.add( new Complex(-1, 0) );
    cs.add( new Complex(-1, -1) );
    cs.add( new Complex(1, -1) );
    cs.add( new Complex(-3, -3) );
    cs.add( new Complex(-3, 3) );
    NumericUtils.complexRealsFirstSorter(cs);
    Complex[] csl = new Complex[ cs.size() ];
    csl[0] = new Complex(-1, 0);
    csl[1] = new Complex(1, 0);
    csl[2] = new Complex(2, 0);
    csl[3] = new Complex(-3, -3);
    csl[4] = new Complex(-3, 3);
    csl[5] = new Complex(-1, -1);
    csl[6] = new Complex(1, -1);
    for (int i = 0; i < csl.length; ++i) {
      assertTrue( cs.get(i).equals(csl[i]) );
    }
  }

  @Test
  public void unwrapAngleDegrees_unchanged() {
    assertEquals(179, NumericUtils.rewrapAngleDegrees(179), 1E-7);
    assertEquals(-179, NumericUtils.rewrapAngleDegrees(-179), 1E-7);
    assertEquals(0, NumericUtils.rewrapAngleDegrees(0), 1E-7);
    assertEquals(45, NumericUtils.rewrapAngleDegrees(45), 1E-7);
    assertEquals(180, NumericUtils.rewrapAngleDegrees(180), 1E-7);
    assertEquals(-180, NumericUtils.rewrapAngleDegrees(-180), 1E-7);
  }

  @Test
  public void unwrapAngleDegrees_higherthan180() {
    assertEquals(167, NumericUtils.rewrapAngleDegrees(527), 1E-7);
    assertEquals(-43, NumericUtils.rewrapAngleDegrees(677), 1E-7);
    assertEquals(0, NumericUtils.rewrapAngleDegrees(360), 1E-7);
    assertEquals(45, NumericUtils.rewrapAngleDegrees(405), 1E-7);
  }

  @Test
  public void unwrapAngleDegrees_lowerthanMinus180() {
    assertEquals(167, NumericUtils.rewrapAngleDegrees(-193), 1E-7);
    assertEquals(-43, NumericUtils.rewrapAngleDegrees(-763), 1E-7);
    assertEquals(0, NumericUtils.rewrapAngleDegrees(-360), 1E-7);
    assertEquals(45, NumericUtils.rewrapAngleDegrees(-1035), 1E-7);
  }
}

