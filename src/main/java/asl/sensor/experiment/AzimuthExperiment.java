package asl.sensor.experiment;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresBuilder;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresOptimizer;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresProblem;
import org.apache.commons.math3.fitting.leastsquares.LevenbergMarquardtOptimizer;
import 
org.apache.commons.math3.fitting.leastsquares.MultivariateJacobianFunction;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import asl.sensor.input.DataBlock;
import asl.sensor.input.DataStore;
import asl.sensor.utils.FFTResult;
import asl.sensor.utils.TimeSeriesUtils;

/**
 * More specific javadoc will be incoming, but for now a brief explanation
 * of the angle conventions used
 * The program attempts to fit known-orthogonal sensors of unknown azimuth to a
 * reference sensor assumed to be north. The rotation angle between the
 * reference sensor and the unknown components is solved for via least-squares
 * using the coherence calculation of the rotated and reference signal.
 * The resulting angle, then, is the clockwise rotation from the reference.
 * If the angle of the reference is zero (i.e., pointing directly north),
 * the result of this calculation SHOULD be the value of the azimuth, using
 * a clockwise rotation convention.
 * If the reference sensor is itself offset X degrees clockwise from
 * north, the azimuth is the sum of the estimated angle difference between
 * the sensors plus the offset from north.
 * This calculation is mostly based on Ringler, Edwards, et al. (2012)
 * 'Relative azimuth inversion by way of damped maximum correlation estimates'
 * but using coherence maximization rather than correlation to find optimized
 * angles.
 * @author akearns
 *
 */
public class AzimuthExperiment extends Experiment {

  final static double TAU = Math.PI * 2;
  
  private double offset = 0.;
  private double angle;
  
  private double[] freqs;
  private double[] coherence;
  
  public AzimuthExperiment() {
    super();

  }
  
  @Override
  protected void backend(final DataStore ds) {
    
    // assume the first two are the reference and the second two are the test?
    
    // we just need four timeseries, don't actually care about response
    DataBlock testNorthBlock = new DataBlock( ds.getXthLoadedBlock(1) );
    DataBlock testEastBlock = new DataBlock( ds.getXthLoadedBlock(2) );
    DataBlock refNorthBlock = new DataBlock( ds.getXthLoadedBlock(3) );
    
    List<Number> testNorth = new ArrayList<Number>( testNorthBlock.getData() );
    String northName = testNorthBlock.getName();
    List<Number> testEast = new ArrayList<Number>( testEastBlock.getData() );
    String eastName = testEastBlock.getName();
    List<Number> refNorth = new ArrayList<Number>( refNorthBlock.getData() );
    String refName = refNorthBlock.getName();
    
    FFTResult.detrend(testNorth);
    FFTResult.detrend(testEast);
    FFTResult.detrend(refNorth);
    
    // originally had normalization step here, but that harmed the estimates
    
    double sps = TimeSeriesUtils.ONE_HZ_INTERVAL / testNorthBlock.getInterval();
    double low = 1./8;
    double high = 1./4;
    
    testNorth = FFTResult.bandFilter(testNorth, sps, low, high);
    testEast = FFTResult.bandFilter(testEast, sps, low, high);
    refNorth = FFTResult.bandFilter(refNorth, sps, low, high);
    
    testNorthBlock.setData(testNorth);
    testEastBlock.setData(testEast);
    refNorthBlock.setData(refNorth);
    
    MultivariateJacobianFunction jacobian = new MultivariateJacobianFunction() {
      
      final DataBlock finalTestNorthBlock = new DataBlock(testNorthBlock);
      final DataBlock finalTestEastBlock = new DataBlock(testEastBlock);
      final DataBlock finalRefNorthBlock = new DataBlock(refNorthBlock);
      
      public Pair<RealVector, RealMatrix> value(final RealVector point) {
        
        return jacobian(point, 
            finalRefNorthBlock, 
            finalTestNorthBlock, 
            finalTestEastBlock);
        
      }
    };
    
    // want mean coherence to be as close to 1 as possible
    RealVector target = MatrixUtils.createRealVector(new double[]{1.});
    
    
    /*
    // first is rel. tolerance, second is abs. tolerance
    ConvergenceChecker<LeastSquaresProblem.Evaluation> cv = 
        new EvaluationRmsChecker(1E-3, 1E-3);
    */
    
    LeastSquaresProblem findAngleY = new LeastSquaresBuilder().
        start(new double[] {0}).
        model(jacobian).
        target(target).
        maxEvaluations(Integer.MAX_VALUE).
        maxIterations(Integer.MAX_VALUE).
        lazyEvaluation(false).
        //checker(cv).
        build();
    
    LeastSquaresOptimizer optimizer = new LevenbergMarquardtOptimizer().
        withCostRelativeTolerance(1E-15).
        withParameterRelativeTolerance(1E-15);
    
    LeastSquaresOptimizer.Optimum optimumY = optimizer.optimize(findAngleY);
    RealVector angleVector = optimumY.getPoint();
    double tempAngle = angleVector.getEntry(0);
    
    System.out.println("Found initial guess for angle");
    
    // angleVector is our new best guess for the azimuth
    // now let's cut the data into 2000-sec windows with 500-sec overlap
    // store the angle and resulting correlation of each window
    // and then take the best-correlation angles and average them
    long start = testNorthBlock.getStartTime();
    long end = testNorthBlock.getEndTime();
    long timeRange = end - start;
    
    // first double -- angle estimate over window
    // second double -- coherence from that estimate over the window
    List<Pair<Double,Double>> angleCoherenceList = 
        new ArrayList<Pair<Double, Double>> ();
    List<Double> sortedCoherence = new ArrayList<Double>();
    
    final long twoThouSecs = 2000L * TimeSeriesUtils.ONE_HZ_INTERVAL; 
    // 1000 ms per second, range length
    final long fiveHundSecs = twoThouSecs / 4L; // distance between windows
    int numWindows = (int) ( (timeRange - twoThouSecs) / fiveHundSecs);
    
    System.out.println("Num. windows for better fit: " + numWindows);
    
    for (int i = 0; i < numWindows; ++i) {
      
      System.out.println("Fitting angle over data in window " + i);
      
      if (timeRange < 2 * twoThouSecs) {
        break;
      }
      
      long wdStart = fiveHundSecs * i + start; // start of 500s-sliding window
      long wdEnd = wdStart + twoThouSecs; // end of window (2000s long)
      
      DataBlock testNorthWindow = new DataBlock(testNorthBlock, wdStart, wdEnd);
      DataBlock testEastWindow = new DataBlock(testEastBlock, wdStart, wdEnd);
      DataBlock refNorthWindow = new DataBlock(refNorthBlock, wdStart, wdEnd);
      
      MultivariateJacobianFunction jbn2 = new MultivariateJacobianFunction() {

        final DataBlock finalTestNorthBlock = testNorthWindow;
        final DataBlock finalTestEastBlock = testEastWindow;
        final DataBlock finalRefNorthBlock = refNorthWindow;

        public Pair<RealVector, RealMatrix> value(final RealVector point) {

          return jacobian(point, 
              finalRefNorthBlock, 
              finalTestNorthBlock, 
              finalTestEastBlock);

        }
      };
      
      LeastSquaresProblem findAngleWindow = new LeastSquaresBuilder().
          start(new double[]{tempAngle}).
          model(jbn2).
          target(target).
          maxEvaluations(Integer.MAX_VALUE).
          maxIterations(Integer.MAX_VALUE).
          lazyEvaluation(false).
          // checker(cv).
          build();
            
      optimumY = optimizer.optimize(findAngleWindow);
      
      RealVector angleVectorWindow = optimumY.getPoint();
      double angleTemp = angleVectorWindow.getEntry(0);
      RealVector resi = optimumY.getResiduals();
      double errorEst = resi.getEntry(0);
      angleCoherenceList.add( new Pair<Double, Double>(angleTemp, errorEst) );
      sortedCoherence.add(errorEst);
    }
    
    if (angleCoherenceList.size() < 1) {
      System.out.println("Too little data for good angle estimation...");
      angle = Math.toDegrees( angleVector.getEntry(0) );
    } else {
      // get the best-coherence estimations of angle and average them
      Collections.sort(sortedCoherence);
      int maxBoundary = Math.max(5, sortedCoherence.size() * 3 / 20);
      sortedCoherence = sortedCoherence.subList(0, maxBoundary);
      Set<Double> acceptableCoherences = new HashSet<Double>(sortedCoherence);
      
      double averageAngle = 0.;
      int coherenceCount = 0;
      
      for (Pair<Double, Double> angCoherePair : angleCoherenceList) {
        double angleTemp = angCoherePair.getFirst();
        double coherence = angCoherePair.getSecond();
        if ( acceptableCoherences.contains(coherence) ) {
          averageAngle += angleTemp;
          ++coherenceCount;
        }
      }
      
      averageAngle /= coherenceCount;
      RealVector angleVec = 
          MatrixUtils.createRealVector(new double[]{averageAngle});
      findAngleY.evaluate(angleVec);
      
      angle = averageAngle;
      
    }

    System.out.println("Found angle");
    
    double angleDeg = Math.toDegrees(angle);
    angleDeg = ( (angleDeg % 360) + 360 ) % 360;
    
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
    
    XYSeries coherenceSeries = new XYSeries("COHERENCE");
    for (int i = 0; i < freqs.length; ++i) {
      coherenceSeries.add(freqs[i], coherence[i]);
    }
    
    xySeriesData.add( new XYSeriesCollection(coherenceSeries) );
    
    
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
  
  public double getOffset() {
    // TODO Auto-generated method stub
    return ( (offset % 360) + 360 ) % 360;
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

  /**
   * Jacobian function for the azimuth solver. Takes in the directional
   * signal components (DataBlocks) and the angle to evaluate at and produces
   * the coherence at that point and the forward difference
   * @param point Current angle
   * @param refNorth Reference sensor, facing north
   * @param testNorth Test sensor, facing approximately north
   * @param testEast Test sensor, facing approximately east and orthogonal to
   * testNorth
   * @return Coherence (RealVector) and forward difference 
   * approximation of the Jacobian (RealMatrix) at the current angle
   */
  public Pair<RealVector, RealMatrix> jacobian(
      final RealVector point, 
      final DataBlock refNorth,
      final DataBlock testNorth, 
      final DataBlock testEast) {
    
    double theta = ( point.getEntry(0) );
    
    double diff = 1E-2;
    
    double lowFreq = 1./18.;
    double highFreq = 1./3.;
    
    DataBlock testRotated = TimeSeriesUtils.rotate(testNorth, testEast, theta);
    
    FFTResult crossPower = FFTResult.spectralCalc(refNorth, testRotated);
    FFTResult rotatedPower = FFTResult.spectralCalc(testRotated, testRotated);
    FFTResult refPower = FFTResult.spectralCalc(refNorth, refNorth);
    
    freqs = crossPower.getFreqs();
    
    Complex[] crossPowerSeries = crossPower.getFFT();
    Complex[] rotatedSeries = rotatedPower.getFFT();
    Complex[] refSeries = refPower.getFFT();
    
    coherence = new double[crossPowerSeries.length];
    
    for (int i = 0; i < crossPowerSeries.length; ++i) {
      Complex conj = crossPowerSeries[i].conjugate();
      Complex numerator = crossPowerSeries[i].multiply(conj);
      Complex denom = rotatedSeries[i].multiply(refSeries[i]);
      coherence[i] = numerator.divide(denom).getReal();
    }
    
    double peakVal = Double.NEGATIVE_INFINITY;
    double peakFreq = 0;
    
    for (int i = 0; i < freqs.length; ++i) {
      if (freqs[i] < lowFreq) {
        continue;
      } else if (freqs[i] > highFreq) {
        break;
      }
      if (peakVal < coherence[i]) {
        peakVal = coherence[i];
        peakFreq = freqs[i];
      }
    }
    
    if (peakFreq / 2 > lowFreq) {
      lowFreq = peakFreq / 2.;
    }
    
    if (peakFreq * 2 < highFreq) {
      highFreq = peakFreq * 2.;
    }
    
    double meanCoherence = 0.;
    int samples = 0;
    
    for (int i = 0; i < freqs.length; ++i) {
      if (freqs[i] < highFreq && freqs[i] > lowFreq) {
        meanCoherence += coherence[i];
        ++samples;
      }
    }
    
    meanCoherence /= samples;
    
    RealVector curValue = 
        MatrixUtils.createRealVector(new double[]{meanCoherence});
    
    double thetaDelta = theta + diff;
    DataBlock rotateDelta = 
        TimeSeriesUtils.rotate(testNorth, testEast, thetaDelta);
    
    crossPower = FFTResult.spectralCalc(refNorth, rotateDelta);
    rotatedPower = FFTResult.spectralCalc(rotateDelta, rotateDelta);
    crossPowerSeries = crossPower.getFFT();
    rotatedSeries = rotatedPower.getFFT();
    
    double fwdMeanCoherence = 0.;
    samples = 0;
    double[] fwdCoherence = new double[crossPowerSeries.length];
    
    for (int i = 0; i < crossPowerSeries.length; ++i) {
      Complex conj = crossPowerSeries[i].conjugate();
      Complex numerator = crossPowerSeries[i].multiply(conj);
      Complex denom = rotatedSeries[i].multiply(refSeries[i]);
      fwdCoherence[i] = numerator.divide(denom).getReal();
      
      if (freqs[i] < highFreq && freqs[i] > lowFreq) {
        fwdMeanCoherence += fwdCoherence[i];
        ++ samples;
      }
      
    }
    
    fwdMeanCoherence /= (double) samples;
    double deltaMean = 0.;
    deltaMean = (fwdMeanCoherence - meanCoherence) / diff;
    
    // System.out.println(deltaMean);
    
    double[][] jacobianArray = new double[1][1];
    jacobianArray[0][0] = deltaMean;
    
    // we have only 1 variable, so jacobian is a matrix w/ single column
    RealMatrix jbn = MatrixUtils.createRealMatrix(jacobianArray);
    
    return new Pair<RealVector, RealMatrix>(curValue, jbn);
  }

  /**
   * Set the angle offset for the reference sensor (degrees from north)
   * @param newOffset Degrees from north that the reference sensor points
   */
  public void setOffset(double newOffset) {
    offset = newOffset;
  }
  

}