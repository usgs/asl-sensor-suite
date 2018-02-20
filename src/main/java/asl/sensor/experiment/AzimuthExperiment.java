package asl.sensor.experiment;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresBuilder;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresOptimizer;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresProblem;
import org.apache.commons.math3.fitting.leastsquares.LevenbergMarquardtOptimizer;
import
org.apache.commons.math3.fitting.leastsquares.MultivariateJacobianFunction;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math3.util.Pair;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import asl.sensor.input.DataBlock;
import asl.sensor.input.DataStore;
import asl.sensor.utils.FFTResult;
import asl.sensor.utils.NumericUtils;
import asl.sensor.utils.TimeSeriesUtils;

/**
 * More specific javadoc will be incoming, but for now a brief explanation
 * of the angle conventions used
 * The program attempts to fit known-orthogonal sensors of unknown azimuth to a
 * reference sensor assumed to be north. The rotation angle between the
 * reference sensor and the unknown components is solved for via least-squares
 * using the correlation calculation of the rotated and reference signal.
 * The resulting angle, then, is the clockwise rotation from the reference.
 * If the angle of the reference is zero (i.e., pointing directly north),
 * the result of this calculation SHOULD be the value of the azimuth, using
 * a clockwise rotation convention.
 * If the reference sensor is itself offset X degrees clockwise from
 * north, the azimuth is the sum of the estimated angle difference between
 * the sensors plus the offset from north.
 * This calculation is mostly based on Ringler, Edwards, et al.,
 * 'Relative azimuth inversion by way of damped maximum correlation estimates',
 * Elsevier Computers and Geosciences 43 (2012)
 * @author akearns
 *
 */
public class AzimuthExperiment extends Experiment {

  private double offset = 0.;

  private double latestCorrelation = 0.; // cache correlation estimates during windowing

  private double angle, uncert;
  // private double[] freqs;

  // private double[] coherence;
  private double[] correlations; // best-fit correlations used to find windows w/ good estimates
  private double minCorr; // lowest correlation value still used in angle estimation
  private boolean simpleCalc; // used for nine-noise calculation
  private boolean enoughPts; // enough points in range for estimation?

  public AzimuthExperiment() {
    super();
    simpleCalc = false;
  }

  /**
   * Entry point for this experiment to guarantee the use of the simple solver
   * and require less overhead, callable from another experiment
   * @param testNorth timeseries data from presumed north-facing test sensor
   * @param testEast timeseries data from presumed east-facing test sensor
   * @param refNorth timeseries data from known north-facing sensor
   * @param interval sampling interval of data
   * @param start start time of data
   * @param end end time of data
   */
  protected void alternateEntryPoint(
      double[] testNorth, double[] testEast,
      double[] refNorth, long interval, long start, long end) {
    dataNames = new ArrayList<String>();
    dataNames.add("N");
    dataNames.add("E");
    dataNames.add("R");
    simpleCalc = true;

    backend(testNorth.clone(), testEast.clone(), refNorth.clone(),
        interval, start, end);
  }

  @Override
  protected void backend(final DataStore ds) {
    // assume the first two are the reference and the second two are the test
    // we just need the timeseries, don't actually care about response

    DataBlock testNorthBlock = ds.getXthLoadedBlock(1);
    DataBlock testEastBlock = ds.getXthLoadedBlock(2);
    DataBlock refNorthBlock = ds.getXthLoadedBlock(3);

    dataNames = new ArrayList<String>();
    dataNames.add( testNorthBlock.getName() );
    dataNames.add( testEastBlock.getName() );
    dataNames.add( refNorthBlock.getName() );

    // resampling should already have been done when loading in data
    long interval = testNorthBlock.getInterval();
    long startTime = testNorthBlock.getStartTime();
    long endTime = testNorthBlock.getEndTime();

    double[] testNorth = testNorthBlock.getData();
    double[] testEast = testEastBlock.getData();
    double[] refNorth = refNorthBlock.getData();

    backend(testNorth, testEast, refNorth, interval, startTime, endTime);

  }

  /**
   * Backend library call for both datasets
   * @param testNorth North-facing data to find azimuth of
   * @param testEast East-facing data to find azimuth of
   * @param refNorth North-facing data to use as reference
   * @param interval Time in nanoseconds between data samples
   * @param startTime Start time of data
   * @param endTime End time of data
   */
  private void backend(
      double[] testNorth, double[] testEast,
      double[] refNorth, long interval, long startTime, long endTime) {

    minCorr = 0; // make sure to initialize
    double tau = NumericUtils.TAU; // 2 pi

    enoughPts = false;

    // does nothing if the data is already 1Hz sample rate
    testNorth = TimeSeriesUtils.decimate(testNorth, interval, TimeSeriesUtils.ONE_HZ_INTERVAL);
    testEast = TimeSeriesUtils.decimate(testEast, interval, TimeSeriesUtils.ONE_HZ_INTERVAL);
    refNorth = TimeSeriesUtils.decimate(refNorth, interval, TimeSeriesUtils.ONE_HZ_INTERVAL);
    // update the actual sample rate if data was above 1Hz sampling
    interval = Math.max(interval, TimeSeriesUtils.ONE_HZ_INTERVAL);

    double[] initTestNorth = TimeSeriesUtils.demean(testNorth);
    double[] initTestEast = TimeSeriesUtils.demean(testEast);
    double[] initRefNorth = TimeSeriesUtils.demean(refNorth);

    initTestNorth = TimeSeriesUtils.detrend(initTestNorth);
    initTestEast = TimeSeriesUtils.detrend(initTestEast);
    initRefNorth = TimeSeriesUtils.detrend(initRefNorth);

    // should there be a normalization step here?

    // data will be downsampled to 1 if > 1Hz rate, else will keep sample rate from input
    double sps = Math.min(1., TimeSeriesUtils.ONE_HZ_INTERVAL / interval);
    double low = 1./8; // filter from 8 seconds interval
    double high = 1./3; // up to 3 seconds interval

    initTestNorth = FFTResult.bandFilter(initTestNorth, sps, low, high);
    initTestEast = FFTResult.bandFilter(initTestEast, sps, low, high);
    initRefNorth = FFTResult.bandFilter(initRefNorth, sps, low, high);


    MultivariateJacobianFunction jacobian =
        getJacobianFunction(initTestNorth, initTestEast, initRefNorth);

    double initAngle = -Math.toRadians(offset);

    LeastSquaresProblem findAngleY = new LeastSquaresBuilder().
        start(new double[] {initAngle}).
        model(jacobian).
        target(new double[]{1}).
        maxEvaluations(Integer.MAX_VALUE).
        maxIterations(Integer.MAX_VALUE).
        lazyEvaluation(false).
        build();

    LeastSquaresOptimizer optimizer = new LevenbergMarquardtOptimizer().
        withCostRelativeTolerance(1E-5).
        withParameterRelativeTolerance(1E-5);

    LeastSquaresOptimizer.Optimum optimumY = optimizer.optimize(findAngleY);
    RealVector angleVector = optimumY.getPoint();
    double tempAngle = angleVector.getEntry(0);
    tempAngle = ( (tempAngle % NumericUtils.TAU) + NumericUtils.TAU)
        % NumericUtils.TAU;

    String newStatus = "Found initial guess for angle: " + tempAngle;
    fireStateChange(newStatus);

    // how much data we need (i.e., iteration length) to check 10 seconds
    // used when checking if alignment is off by 180 degrees
    // int tenSecondsLength = (int)  ( sps * 10 ) + 1;
    // int antipolarTrimLen = tenSecondsLength * 100; // thousand secs?

    if (simpleCalc) {

      // used for orthogonality & multi-component self-noise and gain
      // where a 'pretty good' estimate of the angle is all we need
      // just stop here, don't do windowing
      angle = tempAngle;

      return;
    }

    // angleVector is our new best guess for the azimuth
    // now let's cut the data into 1000-sec windows with 500-sec overlap
    // store the angle and resulting correlation of each window
    // and then take the best-correlation angles and average them
    long timeRange = endTime - startTime;

    // first double -- angle estimate over window
    // second double -- correlation from that estimate over the window
    Map<Long, Pair<Double,Double>> angleCorrelationMap =
        new HashMap<Long, Pair<Double, Double>> ();
    List<Double> sortedCorrelation = new ArrayList<Double>();

    // want (correlation-1+damping) to be as close to 0 as possible
    RealVector target = MatrixUtils.createRealVector(new double[]{0});

    // the best correlation and azimuth angle producing that correlation
    // for the purpose of providing damped estimates
    // (improves susceptibility to noise)
    double bestCorr = jacobian.value(angleVector).getFirst().getEntry(0);
    double bestTheta = tempAngle;
    final long twoThouSecs = 2000L * TimeSeriesUtils.ONE_HZ_INTERVAL;
    // 1000 ms per second, range length
    final long fiveHundSecs = twoThouSecs / 4L; // distance between windows
    int numWindows = (int) ( (timeRange - twoThouSecs) / fiveHundSecs);
    // look at 2000s windows, sliding over 500s of data at a time
    for (int i = 0; i < numWindows; ++i) {
      StringBuilder sb = new StringBuilder();
      sb.append("Fitting angle over data in window ");
      sb.append(i + 1);
      sb.append(" of ");
      sb.append(numWindows);
      newStatus = sb.toString();

      fireStateChange(newStatus);

      // get start and end indices from given times
      long wdStart = fiveHundSecs * i; // start of 500s-sliding window
      long wdEnd = wdStart + twoThouSecs; // end of window (2000s long)

      int startIdx = (int) (wdStart / interval);
      int endIdx = (int) (wdEnd / interval);

      double[] testNorthWin = Arrays.copyOfRange(initTestNorth, startIdx, endIdx);
      double[] testEastWin = Arrays.copyOfRange(initTestEast, startIdx, endIdx);
      double[] refNorthWin = Arrays.copyOfRange(initRefNorth, startIdx, endIdx);

      testNorthWin = TimeSeriesUtils.detrend(testNorthWin);
      testEastWin = TimeSeriesUtils.detrend(testEastWin);
      refNorthWin = TimeSeriesUtils.detrend(refNorthWin);

      testNorthWin = FFTResult.bandFilter(testNorthWin, sps, low, high);
      testEastWin = FFTResult.bandFilter(testEastWin, sps, low, high);
      refNorthWin = FFTResult.bandFilter(refNorthWin, sps, low, high);

      jacobian =
          getDampedJacobianFunction(testNorthWin, testEastWin, refNorthWin, bestCorr, bestTheta);

      LeastSquaresProblem findAngleWindow = new LeastSquaresBuilder().
          start(new double[]{bestTheta}).
          model(jacobian).
          target(target).
          maxEvaluations(Integer.MAX_VALUE).
          maxIterations(Integer.MAX_VALUE).
          lazyEvaluation(false).
          build();

      optimumY = optimizer.optimize(findAngleWindow);

      RealVector angleVectorWindow = optimumY.getPoint();
      findAngleWindow.evaluate(angleVectorWindow);
      // call to evaluate at best-fit point gives corresponding latestCorrelation as side effect
      double angleTemp = angleVectorWindow.getEntry(0);

      angleTemp = ( (angleTemp % tau) + tau ) % tau;

      double correlation = latestCorrelation;

      if (correlation > bestCorr) {
        bestCorr = correlation;
        bestTheta = angleTemp;
      }

      angleCorrelationMap.put(
          wdStart, new Pair<Double, Double>(angleTemp, correlation) );
      sortedCorrelation.add(correlation);
    }

    int minCorrelations = 5;
    double[] angles = new double[]{};
    correlations = new double[]{};
    // TODO: can refactor this to break out if numWindows is < 5
    if (angleCorrelationMap.size() < minCorrelations) {
      fireStateChange("Window size too small for good angle estimation...");
      angle = tempAngle % tau; // tempAngle is the initial estimate from before windowing occurs
    } else {
      // get the best-correlation estimations of angle and average them
      enoughPts = true;
      Collections.sort(sortedCorrelation);
      Collections.reverse(sortedCorrelation); // sort from best to worst
      // get either exactly top 5 or (if larger), top 15% of estimates
      int maxBoundary = Math.max(minCorrelations, sortedCorrelation.size() * 3 / 20);
      // start from 0 because sort is descending order
      sortedCorrelation = sortedCorrelation.subList(0, maxBoundary);
      minCorr = sortedCorrelation.get( sortedCorrelation.size() - 1);

      // store good values for use in std dev calculation
      List<Double> acceptedAngles = new ArrayList<Double>();

      // deal with wraparound issue
      correlations = new double[angleCorrelationMap.size()];
      angles = new double[angleCorrelationMap.size()];


      List<Long> times = new ArrayList<Long>( angleCorrelationMap.keySet() );
      Collections.sort(times);

      for (int i = 0; i < times.size(); ++i) {
        long time = times.get(i);
        Pair<Double, Double> angleAndCorrelation = angleCorrelationMap.get(time);
        angles[i] = angleAndCorrelation.getFirst();
        correlations[i] = angleAndCorrelation.getSecond();
      }

      // shift term here used to deal with potential discontinuities in the mean of the data
      double shift = angles[0] + Math.PI/4; // 45 degrees offset
      for (int i = 0; i < angles.length; ++i) {
        angles[i] = ( ( (angles[i] + shift) % tau ) + tau ) % tau;
      }

      // now get the average
      double averageAngle = 0.;
      for (int i = 0; i < angles.length; ++i) {
        angles[i] -= shift; // subtract shift term back out
        double correlation = correlations[i];
        if ( correlation >= minCorr && acceptedAngles.size() < maxBoundary ) {
          // don't keep adding angles once we've hit the max size
          // don't break out of the loop, though, so we can remove the shift from everything
          acceptedAngles.add(angles[i]);
          averageAngle += angles[i];
        }
      }
      averageAngle /= acceptedAngles.size();

      uncert = 0.;

      // now get standard deviation
      for (double angle : acceptedAngles) {
        uncert += Math.pow(angle - averageAngle, 2);
      }

      uncert = Math.sqrt( uncert / acceptedAngles.size() );
      uncert *= 2; // two-sigma gets us 95% confidence interval

      // do this calculation to get plot of freq/correlation, a side effect
      // of running evaluation at the given point; this will be plotted
      RealVector angleVec =
          MatrixUtils.createRealVector(new double[]{averageAngle});
      findAngleY.evaluate(angleVec);

      angle = ( (averageAngle % tau) + tau ) % tau;

    }

    fireStateChange("Solver completed! Producing plots...");

    double angleDeg = Math.toDegrees(angle);

    String northName = dataNames.get(0);
    String eastName = dataNames.get(1);
    String refName = dataNames.get(2);

    XYSeries ref = new XYSeries(northName + " rel. to reference");
    ref.add(offset + angleDeg, 0);
    ref.add(offset + angleDeg, 1);
    XYSeries set = new XYSeries(eastName + " rel. to reference");
    set.add(offset + angleDeg + 90, 1);
    set.add(offset + angleDeg + 90, 0);
    XYSeries fromNorth = new XYSeries (refName + " location");
    fromNorth.add(offset, 1);
    fromNorth.add(offset, 0);

    // xySeriesData = new XYSeriesCollection();
    XYSeriesCollection xysc = new XYSeriesCollection();
    xysc.addSeries(ref);
    xysc.addSeries(set);
    xysc.addSeries(fromNorth);
    xySeriesData.add(xysc);

    xysc = new XYSeriesCollection();
    XYSeries timeMapAngle = new XYSeries("Best-fit angle per window (not including ref. shift)");
    XYSeries timeMapCorrelation = new XYSeries("Correlation estimate per window");
    xysc.addSeries(timeMapAngle);
    xysc.addSeries(timeMapCorrelation);

    for (int i = 0; i < angles.length; ++i) {
        long xVal = i * 500;
        double angle = angles[i];
        // angle += Math.toRadians(offset);
        angle = (angle % tau);
        double correlation = correlations[i];
        timeMapCorrelation.add(xVal, correlation);
        timeMapAngle.add( xVal, Math.toDegrees(angle) );
    }


    xySeriesData.add( new XYSeriesCollection(timeMapAngle) );
    xySeriesData.add( new XYSeriesCollection(timeMapCorrelation) );
  }

  @Override
  public int blocksNeeded() {
    return 3;
  }

  /**
   * Return the fit angle calculated by the backend in degrees
   * @return angle result in degrees
   */
  public double getFitAngle() {
    return Math.toDegrees(angle);
  }

  /**
   * Return the fit angle calculated by the backend in radians
   * @return angle result in radians
   */
  public double getFitAngleRad() {
    return angle;
  }

  /**
   * Returns the minimum acceptable correlation used in angle estimates
   * @return correlation cut-off point
   */
  public double getMinCorr() {
    return minCorr;
  }

  /**
   * Returns the jacobian function for initial estimate given input timeseries data.
   * The timeseries are used as input to the rotation function.
   * We take the inputs as fixed and rotate copies of the data to find the
   * Jacobian of the data.
   * @param l1 Data from the test sensor's north-facing component
   * @param l2 Data from the test sensor's east-facing component
   * @param l3 Data from the known north-facing sensor
   * @return jacobian function to fit an angle of max correlation of this data
   */
  private MultivariateJacobianFunction
  getJacobianFunction(double[] l1, double[] l2, double[] l3) {

    // make my func the j-func, I want that func-y stuff
    // make my func the j-func, I want that func-y stuff
    MultivariateJacobianFunction jFunc = new MultivariateJacobianFunction() {

      final double[] finalTestNorth = l1;
      final double[] finalTestEast = l2;
      final double[] finalRefNorth = l3;

      @Override
      public Pair<RealVector, RealMatrix> value(final RealVector point) {
        return jacobian(point,
            finalRefNorth,
            finalTestNorth,
            finalTestEast);
      }
    };
    return jFunc;
  }

  // This is the damped jacobian function for windowed estimates
  // we use a different cost function for initial estimate since using the
  // squared correlation would make x, 180+x produce the same values
  private MultivariateJacobianFunction
  getDampedJacobianFunction(double[] l1, double[] l2, double[] l3, double cr, double th) {

    // make my func the j-func, I want that func-y stuff
    MultivariateJacobianFunction jFunc = new MultivariateJacobianFunction() {

      final double[] finalTestNorth = l1;
      final double[] finalTestEast = l2;
      final double[] finalRefNorth = l3;
      final double bestCorr = cr;
      final double bestTheta = th;

      @Override
      public Pair<RealVector, RealMatrix> value(final RealVector point) {
        return jacobian(point,
            finalRefNorth,
            finalTestNorth,
            finalTestEast,
            bestCorr,
            bestTheta);
      }

    };

    return jFunc;
  }

  public double getOffset() {
    return ( (offset % 360) + 360 ) % 360;
  }

  /**
   * Get the uncertainty of the angle
   * @return Uncertainty estimation of the current angle (from variance)
   */
  public double getUncertainty() {
    return Math.toDegrees(uncert);
  }

  @Override
  public boolean hasEnoughData(DataStore ds) {
    for (int i = 0; i < blocksNeeded(); ++i) {
      if ( !ds.blockIsSet(i) ) {
        return false;
      }
    }
    return true;
  }

  public double[] getCorrelations() {
    return correlations;
  }

  /**
   * Jacobian function for the azimuth solver. Takes in the directional
   * signal components (DataBlocks) and the angle to evaluate at and produces
   * a damped cost function based on best correlation (and the angle producing it) from previous
   * data windows, each one fixed for a given window (not changing on recursive calls for the same
   * set of input data).
   * @param point Current angle
   * @param refNorth Reference sensor, facing north
   * @param testNorth Test sensor, facing approximately north
   * @param testEast Test sensor, facing approximately east and orthogonal to
   * testNorth
   * @param bestCorr Most recent best-result value for the correlation (from previous windows)
   * @param bestTheta Most recent best-result value for the correlation (from previous windows)
   * @return Correlation (RealVector) and forward difference
   * approximation of the Jacobian (RealMatrix) at the current angle
   */
  private Pair<RealVector, RealMatrix> jacobian(
      final RealVector point,
      final double[] refNorth,
      final double[] testNorth,
      final double[] testEast,
      final double bestCorr,
      final double bestTheta) {

    double diff = 1E-12;

    double theta = ( point.getEntry(0) );
    double thetaDelta = theta + diff;

    // was the frequency range under examination (in Hz) when doing coherence
    // double lowFreq = 1./18.;
    // double highFreq = 1./3.;

    // angles of rotation are x, x+dx respectively
    double[] testRotated =
        TimeSeriesUtils.rotate(testNorth, testEast, theta);
    double[] rotatedDiff =
        TimeSeriesUtils.rotate(testNorth, testEast, thetaDelta);

    PearsonsCorrelation pc = new PearsonsCorrelation();
    double value = pc.correlation(refNorth, testRotated);
    latestCorrelation = value;
    double damping = (bestCorr - 1) * (theta - bestTheta);
    value = Math.pow(value - 1 + damping, 2);
    RealVector valueVec = MatrixUtils.createRealVector(new double[]{value});
    double deltaY = pc.correlation(refNorth, rotatedDiff);
    damping = (bestCorr - 1) * (thetaDelta - bestTheta);
    deltaY = Math.pow(deltaY - 1 + damping, 2);
    double change = (deltaY - value) / diff;
    double[][] jacobianArray = new double[][]{{change}};
    RealMatrix jacobian = MatrixUtils.createRealMatrix(jacobianArray);
    return new Pair<RealVector, RealMatrix>(valueVec, jacobian);

  }

  private Pair<RealVector, RealMatrix> jacobian(
      final RealVector point,
      final double[] refNorth,
      final double[] testNorth,
      final double[] testEast) {

    double diff = 1E-12;

    double theta = ( point.getEntry(0) );
    double thetaDelta = theta + diff;

    // angles of rotation are x, x+dx respectively
    double[] testRotated =
        TimeSeriesUtils.rotate(testNorth, testEast, theta);
    double[] rotatedDiff =
        TimeSeriesUtils.rotate(testNorth, testEast, thetaDelta);

    PearsonsCorrelation pc = new PearsonsCorrelation();
    double value = pc.correlation(refNorth, testRotated);
    RealVector valueVec = MatrixUtils.createRealVector(new double[]{value});
    double deltaY = pc.correlation(refNorth, rotatedDiff);
    double change = (deltaY - value) / diff;
    double[][] jacobianArray = new double[][]{{change}};
    RealMatrix jacobian = MatrixUtils.createRealMatrix(jacobianArray);
    return new Pair<RealVector, RealMatrix>(valueVec, jacobian);

  }

  /**
   * Set the angle offset for the reference sensor (degrees from north)
   * @param newOffset Degrees from north that the reference sensor points
   */
  public void setOffset(double newOffset) {
    offset = newOffset;
  }

  /**
   * Used to set a simple calculation of rotation angle, such as for
   * nine-input self-noise. This is the case where the additional windowing
   * is NOT done, and the initial least-squares guess gives us an answer.
   * When creating an instance of this object, this is set to false and only
   * needs to be explicitly set when a simple calculation is desired.
   * This is used primarily when aligning data is done as part of a multi-input experiment
   * (i.e., 9-input self-noise, 6-input relative gain)
   * @param isSimple True if a simple calculation should be done
   */
  public void setSimple(boolean isSimple) {
    simpleCalc = isSimple;
  }

  /**
   * Returns true if there were enough points to do correlation windowing step
   * @return Boolean that is true if correlation windows were taken
   */
  public boolean hadEnoughPoints() {
    return enoughPts;
  }
}