package asl.sensor.experiment;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresBuilder;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresOptimizer;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresProblem;
import org.apache.commons.math3.fitting.leastsquares.LevenbergMarquardtOptimizer;
import org.apache.commons.math3.fitting.leastsquares.MultivariateJacobianFunction;
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

public class AzimuthExperiment extends Experiment {

  final static double TAU = Math.PI * 2;
  
  double offset = 0.;
  double angle;
  
  public AzimuthExperiment() {
    super();

  }
  
  public void setOffset(double newOffset) {
    offset = newOffset;
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

    
    testNorth = TimeSeriesUtils.normalize(testNorth);
    testEast = TimeSeriesUtils.normalize(testEast);
    refNorth = TimeSeriesUtils.normalize(refNorth);
    
    double sps = TimeSeriesUtils.ONE_HZ_INTERVAL / testNorthBlock.getInterval();
    double low = 1./8;
    double high = 1./4;
    
    testNorth = FFTResult.bandFilter(testNorth, sps, low, high);
    testEast = FFTResult.bandFilter(testEast, sps, low, high);
    refNorth = FFTResult.bandFilter(refNorth, sps, low, high);
    
    int len = testNorth.size();
    
    double[] testYArr = new double[len];
    double[] testXArr = new double[len];
    double[] refYArr = new double[len];
    
    for (int i = 0; i < len; ++i) {
      testYArr[i] = testNorth.get(i).doubleValue();
      testXArr[i] = testEast.get(i).doubleValue();
      refYArr[i] = refNorth.get(i).doubleValue();
    }
    
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
            finalTestEastBlock, 
            len);
        
      }
    };
    
    // want mean coherence to be as close to 1 as possible
    RealVector target = MatrixUtils.createRealVector(new double[]{1.});
    
    LeastSquaresProblem findAngleY = new LeastSquaresBuilder().
        start(new double[] {0}).
        model(jacobian).
        target(target).
        maxEvaluations(Integer.MAX_VALUE).
        maxIterations(Integer.MAX_VALUE).
        lazyEvaluation(false).
        build();
        
    LeastSquaresOptimizer optimizer = new LevenbergMarquardtOptimizer().
        withCostRelativeTolerance(1.0E-10).
        withParameterRelativeTolerance(1.0E-10);
    
    LeastSquaresOptimizer.Optimum optimumY = optimizer.optimize(findAngleY);
    RealVector angleVector = optimumY.getPoint();
    
    angle = Math.toDegrees( angleVector.getEntry(0) );
    angle = ( (angle % 360) + 360 ) % 360;
    // allows us to have 
    if (angle > 180) {
      angle = (360-angle) % 360;
    }
    
    XYSeries ref = new XYSeries(northName + " rel. to reference");
    ref.add(offset - angle, 0);
    ref.add(offset - angle, 1);
    XYSeries set = new XYSeries(eastName + " rel. to reference");
    set.add(offset - angle + 90, 1);
    set.add(offset - angle + 90, 0);
    XYSeries fromNorth = new XYSeries (refName + " location");
    fromNorth.add(offset, 1);
    fromNorth.add(offset, 0);

    xySeriesData = new XYSeriesCollection();
    xySeriesData.addSeries(ref);
    xySeriesData.addSeries(set);
    xySeriesData.addSeries(fromNorth);
    
    /*
    XYSeries diffSrs = new XYSeries("Diff(" + testName + ", " + refName + ")");
    XYSeries diffRotSrs = new XYSeries("Diff(" + testName + ", Rotated Ref.)");
    
    
    RealVector diffLH1 = testY.subtract(refY);
    RealVector diffComponents = jacobian.value(angleVector).getFirst();
    
    for (int i = 0; i < len; ++i) {
      diffSrs.add (timeAtPoint, diffLH1.getEntry(i) );
      diffRotSrs.add( timeAtPoint, diffComponents.getEntry(i) );
      
      timeAtPoint += tick;
    }
    
    XYSeriesCollection xysc = (XYSeriesCollection) xySeriesData;
    
    xysc.addSeries(diffSrs);
    xysc.addSeries(diffRotSrs);
    */
  }
  
  public Pair<RealVector, RealMatrix> jacobian(
      final RealVector point, 
      final DataBlock refNorth,
      final DataBlock testNorth, 
      final DataBlock testEast, 
      int len) {
    
    
    double theta = ( point.getEntry(0) );
    
    double diff = 1E-12;
    
    double lowFreq = 1./18.;
    double highFreq = 1./3.;
    
    DataBlock testRotated = TimeSeriesUtils.rotate(testNorth, testEast, theta);
        
    FFTResult crossPower = FFTResult.spectralCalc(refNorth, testRotated);
    FFTResult rotatedPower = FFTResult.spectralCalc(testRotated, testRotated);
    FFTResult refPower = FFTResult.spectralCalc(refNorth, refNorth);
    
    double[] freqs = crossPower.getFreqs();
    
    Complex[] crossPowerSeries = crossPower.getFFT();
    Complex[] rotatedSeries = rotatedPower.getFFT();
    Complex[] refSeries = refPower.getFFT();
    
    double[] coherence = new double[crossPowerSeries.length];
    
    for (int i = 0; i < crossPowerSeries.length; ++i) {
      // 
      Complex conj = crossPowerSeries[i].conjugate();
      Complex numerator = crossPowerSeries[i].multiply(conj);
      Complex denom = rotatedSeries[i].multiply(refSeries[i]);
      coherence[i] = numerator.divide(denom).getReal();
    }
    
    double peakVal = Double.NEGATIVE_INFINITY;
    double peakFreq = 0;
    
    for (int i = 0; i < freqs.length; ++i) {
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
    double[] fwdCoherence = new double[crossPowerSeries.length];
    
    for (int i = 0; i < crossPowerSeries.length; ++i) {
      Complex conj = crossPowerSeries[i].conjugate();
      Complex numerator = crossPowerSeries[i].multiply(conj);
      Complex denom = rotatedSeries[i].multiply(refSeries[i]);
      fwdCoherence[i] = numerator.divide(denom).getReal();
      
      if (freqs[i] < highFreq && freqs[i] > lowFreq) {
        fwdMeanCoherence += fwdCoherence[i];
      }
      
    }
    
    fwdMeanCoherence /= (double) samples;
    double deltaMean = 0.;
    deltaMean = (fwdMeanCoherence - meanCoherence) / (thetaDelta - theta);
    
    double[][] jacobianArray = new double[1][1];
    jacobianArray[0][0] = deltaMean;
    
    // we have only 1 variable, so jacobian is a matrix w/ single column
    RealMatrix jbn = MatrixUtils.createRealMatrix(jacobianArray);
    
    return new Pair<RealVector, RealMatrix>(curValue, jbn);
  }
  
  public double getFitAngle() {
    return angle;
  }
  
  @Override
  public boolean hasEnoughData(DataStore ds) {
    for (int i = 0; i < 3; ++i) {
      if ( !ds.blockIsSet(i) ) {
        return false;
      }
    }
    return true;
  }
  

}