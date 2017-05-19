package asl.sensor.experiment;

import java.util.ArrayList;
import java.util.Arrays;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.fitting.leastsquares.EvaluationRmsChecker;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresBuilder;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresOptimizer;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresProblem;
import 
org.apache.commons.math3.fitting.leastsquares.LevenbergMarquardtOptimizer;
import 
org.apache.commons.math3.fitting.leastsquares.MultivariateJacobianFunction;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.optim.ConvergenceChecker;
import org.apache.commons.math3.util.Pair;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import asl.sensor.input.DataBlock;
import asl.sensor.input.DataStore;
import asl.sensor.input.InstrumentResponse;
import asl.sensor.utils.FFTResult;
import asl.sensor.utils.NumericUtils;
import asl.sensor.utils.TimeSeriesUtils;

/**
 * Basic outline of what this does
 * Given a seed file timeseries representing the calibration function passed
 * into a sensor (placed in first timeseries input in the inputpanel / 
 * datastore) and one of the output of that function by the sensor (in the
 * following timeseries input) with a corresponding response, this function
 * calculates corner frequency and damping parameters from the response,
 * deconvolves that from the sensor's output (using frequency space, thus
 * doable by simple division of the FFT of the sensor data divided by the
 * response FFT), and then produces a normalized function that should match
 * the step function.
 * In addition to finding the output produced by the response file, this
 * experiment also uses the Apache Commons least-squares solver with a
 * Levenberg-Marquardt optimizer in order to find the best-fit parameters for
 * the corner and damping. In order to do this, the function defines a further
 * backend function that calculates the forward-approximation of the derivative
 * of the function given a sample point.
 * In addition to the deconvolved and best-fit timeseries, this class is
 * also able to return the values for the corner frequency and damping that 
 * produce them, as well as their corresponding residual values.
 * @author akearns
 *
 */
public class StepExperiment extends Experiment{

  private double f, h; //corner and damping of output (uncorrected)
  private double fCorr, hCorr; // fit parameters to turn output into cal input
  private double initResid, fitResid; // residual values
  private double sps; // samples per second
  
  private int trimmedLength, cutAmount;
  private double[] freqs;
  private Complex[] sensorFFTSeries; // FFT of step cal from sensor
  private double[] stepCalSeries; // time series of raw step cal function
  
  private int sensorOutIdx;
  
  final double STEP_FACTOR = 1E-16;
  
  public StepExperiment() {
    super();
  }
  
  @Override
  protected void backend(final DataStore ds) {
    
    dataNames = new ArrayList<String>();
    
    // assume that the first block is the raw step calibration
    // the raw calibration is defined as not having an associated response
    DataBlock stepCalRaw = ds.getXthLoadedBlock(1);
    dataNames.add( stepCalRaw.getName() );
    
    stepCalSeries = new double[stepCalRaw.size()];
    
    // resulting inverse fft should be this length
    trimmedLength = stepCalRaw.size();
    long interval = stepCalRaw.getInterval();
    
    double[] stepCalUnfiltered = new double[stepCalRaw.size()];
    for (int i = 0; i < trimmedLength; ++i) {
      stepCalUnfiltered[i] = stepCalRaw.getData().get(i).doubleValue();
    }
    
    sps = TimeSeriesUtils.ONE_HZ_INTERVAL / interval;
    
    double[] stepCalFiltered = 
        FFTResult.bandFilter(stepCalUnfiltered, sps, 0., 0.1);
    
    // delete the last 10 seconds from the end (has filtering artifacts)
    
    cutAmount = (int) sps * 10;
    int highBound = stepCalFiltered.length - cutAmount;
    trimmedLength = highBound - cutAmount; // length - 2 * cutAmount
    
    double[] stepCalTrimmed = 
        Arrays.copyOfRange(stepCalFiltered, cutAmount, highBound);
    
    stepCalSeries = FFTResult.demean(stepCalTrimmed);
    
    stepCalSeries = TimeSeriesUtils.normalize(stepCalSeries);
    
    // FFTResult.detrend(toDetrend);
    
    // stepCalSeries = TimeSeriesUtils.normalize(filteredStepCal);
    // stepCalSeries = filteredStepCal;
    
    // but we want the response and the data of the cal result
    sensorOutIdx = ds.getXthFullyLoadedIndex(1);
    
    // if someone did load a raw cal with the response, then we wouldn't
    // get a different block with the second call above, so we get the 
    // next loaded block/response pair
    if ( ds.getBlock(sensorOutIdx).getName().equals( stepCalRaw.getName() ) ) {
      sensorOutIdx = ds.getXthFullyLoadedIndex(2);
    }

    // get data of the result of the step calibration
    DataBlock sensorOutput = ds.getBlock(sensorOutIdx);
    dataNames.add( sensorOutput.getName() );
    
    // long interval = sensorOutput.getInterval();
    InstrumentResponse ir = ds.getResponse(sensorOutIdx);
    dataNames.add( ir.getName() );
    Complex pole = ir.getPoles().get(0);
    
    f = 1. / (NumericUtils.TAU / pole.abs() ); // corner frequency
    h = Math.abs( pole.getReal() / pole.abs() ); // damping
    
    // these manually-set parameters were used in testing convergence
    // f = 0.002777;
    // h = 0.707107;
    
    double[] params = new double[]{f, h};
    
    // get FFT of datablock timeseries, deconvolve with response
    // single sided FFT includes lowpass filter, demean, and normalization
    FFTResult sensorsFFT = FFTResult.singleSidedFFT(sensorOutput);
    // these values used in calculating the response deconvolution
    sensorFFTSeries = sensorsFFT.getFFT();
    freqs = sensorsFFT.getFreqs();
    
    double[] toPlot = calculate(params);
    
    long start = 0L; // was db.startTime();
    long now = start;
    XYSeries xys = new XYSeries("STEP *^(-1) RESP");
    XYSeries scs = new XYSeries( stepCalRaw.getName() ); 
    for (double point : toPlot) {
      double seconds = (double) now / TimeSeriesUtils.ONE_HZ_INTERVAL;
      xys.add(seconds, point);
      now += interval;
    }
    now = start;
    for (Number point : stepCalSeries) {
      double seconds = (double) now / TimeSeriesUtils.ONE_HZ_INTERVAL;
      scs.add(seconds, point);
      now += interval;
    }
    
    // next we'll want to find the parameters to fit the plots
    // to the inputted data
    XYSeriesCollection xysc = new XYSeriesCollection();
    xysc.addSeries(scs);
    xysc.addSeries(xys);
    
    // next step: curve fitting
    RealVector startVector = MatrixUtils.createRealVector(params);
    RealVector observedComponents = MatrixUtils.createRealVector(stepCalSeries);
    
    ConvergenceChecker<LeastSquaresProblem.Evaluation> svc = 
        new EvaluationRmsChecker(1E-50, 1E-50);
    
    // used to fit parameters
    MultivariateJacobianFunction jbn = new MultivariateJacobianFunction() {
      
      public Pair<RealVector, RealMatrix> value(RealVector point) {
          return jacobian(point);
      }
      
    };
    
    LeastSquaresProblem lsp = new LeastSquaresBuilder().
        start(startVector).
        target(observedComponents).
        model(jbn).
        lazyEvaluation(false).
        maxEvaluations(Integer.MAX_VALUE).
        maxIterations(Integer.MAX_VALUE).
        checker(svc).
        build();
    
    LeastSquaresProblem.Evaluation initEval = lsp.evaluate(startVector);
    // System.out.println("INITIAL GUESS RESIDUAL: " +  initEval.getRMS() );
    
    initResid = initEval.getRMS() * 100;
    
    LeastSquaresOptimizer optimizer = new LevenbergMarquardtOptimizer().
        withCostRelativeTolerance(1.0E-15).
        withParameterRelativeTolerance(1.0E-15);
    
    LeastSquaresOptimizer.Optimum optimum = optimizer.optimize(lsp);
    
    // System.out.println("FIT PARAMS RESIDUAL: " +  optimum.getRMS() );
    
    fitResid = optimum.getRMS() * 100;
    
    double[] newParams = optimum.getPoint().toArray();
    
    fCorr = newParams[0];
    hCorr = newParams[1];
    
    double[] fitPlot = calculate(newParams);
    // fitPlot = TimeSeriesUtils.normalize(fitPlot);
    start = 0L; // was db.startTime();
    now = start;
    XYSeries bfs = new XYSeries("BEST FIT PLOT");
    for (double point : fitPlot) {
      double seconds = (double) now / TimeSeriesUtils.ONE_HZ_INTERVAL;
      bfs.add(seconds, point);
      now += interval;
    }
    
    xysc.addSeries(bfs);

    // add plot of step stuff
    xySeriesData.add(xysc);

    fireStateChange("Fit gotten. Getting Bode plots...");
    
    // now add the plots of response curve and magnitude from init & fit value
    
    // p1 = -(h+i*sqrt(1-h^2))*2*pi*f
    Complex p1 = new Complex( h, Math.sqrt( 1 - Math.pow(h, 2) ) );
    p1 = p1.multiply( -1 * NumericUtils.TAU * f );
    Complex p2 = new Complex( h, -1 * Math.sqrt( 1 - Math.pow(h, 2) ) );
    p2 = p2.multiply( -1 * NumericUtils.TAU * f);
    
    InstrumentResponse fitResp = new InstrumentResponse(ir);
    fitResp.setName( fitResp.getName() + " [FIT]" );
    
    Complex[] inputCurve = ir.applyResponseToInput(freqs);
    Complex[] fitCurve = fitResp.applyResponseToInput(freqs);
    
    XYSeries inMag = new XYSeries( ir.getName() + " " + " magnitude" );
    XYSeries inPhase = new XYSeries ( ir.getName() + " " + " phase" );
    XYSeries fitMag = new XYSeries( fitResp.getName() + " " + " magnitude" );
    XYSeries fitPhase = new XYSeries ( fitResp.getName() + " " + " phase" );
    
    double phiPrevIn = .0;
    double phiPrevFit = .0;
    for (int i = 0; i < freqs.length; ++i) {
      
      if (freqs[i] == 0.) {
        continue;
      }
      
      Complex tmpIn = inputCurve[i].divide(NumericUtils.TAU * freqs[i]);
      Complex tmpFit = fitCurve[i].divide(NumericUtils.TAU * freqs[i]);
      
      double phiIn = NumericUtils.atanc(tmpIn);
      phiIn = NumericUtils.unwrap(phiIn, phiPrevIn);
      phiPrevIn = phiIn;
      phiIn = Math.toDegrees(phiIn);
      
      double phiFit = NumericUtils.atanc(tmpIn);
      phiFit = NumericUtils.unwrap(phiFit, phiPrevFit);
      phiPrevFit = phiFit;
      phiFit = Math.toDegrees(phiFit);
      
      double magAccelIn = tmpIn.abs();
      inMag.add( freqs[i], 10 * Math.log10(magAccelIn) );
      inPhase.add(freqs[i], phiIn);
      
      double magAccelFit = tmpFit.abs();
      fitMag.add( freqs[i], 10 * Math.log10(magAccelFit) );
      fitPhase.add(freqs[i], phiFit);
    }
    
    xysc = new XYSeriesCollection();
    xysc.addSeries(inMag);
    xysc.addSeries(fitMag);
    xySeriesData.add(xysc);
    
    XYSeriesCollection phaseCollection = new XYSeriesCollection();
    phaseCollection.addSeries(inPhase);
    phaseCollection.addSeries(fitPhase);
    xySeriesData.add(phaseCollection);
    
  }
  
  @Override
  public int blocksNeeded() {
    return 2;
  }
  
  /**
   * Does the deconvolution of the response calculated from the corner freq. (f)
   * and damping (h) parameters passed in
   * @param beta Double array of form {f,h}
   * @return The timeseries resulting from deconvolution of the calculated
   * response from the sensor-input timeseries (done in frequency space)
   */
  public double[] calculate(double[] beta) {
    
    // the original length of the timeseries data we've gotten the FFT of
    int inverseTrim = trimmedLength + 2 * cutAmount;
    // the upper bound on the returned data
    int upperBound = inverseTrim - cutAmount;
    
    double f = beta[0];
    double h = beta[1];
    
    // use these terms to calculate the poles
    
    // term inside the square root in the calculations of p1, p2
    // (h^2-1)
    Complex tempResult = new Complex( Math.pow(h, 2) - 1 );
    
    double omega = 2 * Math.PI * f; // omega_0
    
    // casting h as a complex number to enable these calculations
    Complex hCast = new Complex(h);
    
    // - (h + sqrt(h^2-1))
    Complex pole1 = hCast.add( tempResult.sqrt() ).multiply(-1);
    pole1 = pole1.multiply(omega); 
      // cannot modify complex numbers "in-place"; they MUST be (re)assigned
    
    // - (h - sqrt(h^2-1))
    Complex pole2 = hCast.subtract( tempResult.sqrt() ).multiply(-1);
    pole2 = pole2.multiply(omega);
    
    // calculate the FFT of the response
    Complex[] respFFT = new Complex[sensorFFTSeries.length]; // array of resps
    double max = 0.0;
    // don't let denominator be zero
    respFFT[0] = Complex.ONE;
    for (int i = 1; i < respFFT.length; ++i) {

      // replaced freqs[i] with 
      // 2*pi*i*f
      Complex factor = new Complex(0, 2*Math.PI*freqs[i]); 
      
      // (2*pi*i*f - p1) * (2*pi*f*i - p2)
      Complex denom = factor.subtract(pole1).multiply( factor.subtract(pole2) );
      
      respFFT[i] = factor.divide(denom);
      
      
      if (respFFT[i].abs() > max) {
        max = respFFT[i].abs();
      }
      
    }
    
    Complex[] toDeconvolve = new Complex[sensorFFTSeries.length];
    
    // deconvolving response is dividing fft(signal) by fft(response)
    
    for (int i = 0; i < sensorFFTSeries.length; ++i) {
      // the conjugate of the response, used twice in deconvolution
      Complex conjResp = respFFT[i].conjugate();
      double aboveZero = 0.008*max;  // term to keep the denominator above 0
      
      // resp * conj(resp) + 0.008(max(|resp|))
      Complex denom = respFFT[i].multiply(conjResp).add(aboveZero);
      
      // deconvolving the response from output
      // fft * conj(resp) / (resp * conjResp)+0.008(max(|resp|))
      toDeconvolve[i] = sensorFFTSeries[i].multiply(conjResp).divide(denom);
    }
    
    double[] returnValue =  
        FFTResult.singleSidedInverseFFT(toDeconvolve, inverseTrim);
    
    returnValue = FFTResult.demean(returnValue);
    
    // attempt to filter out additional noise
    returnValue = 
        FFTResult.bandFilterWithCuts(returnValue, sps, 0.0, 0.1, 0.0, sps);
    
    // trim out the ends
    returnValue = Arrays.copyOfRange(returnValue, cutAmount, upperBound);
    
    // normalize in order to fit the plot    
    return TimeSeriesUtils.normalize(returnValue);
    
    
  }
  
  /**
   * Returns the corner, damping, and residual values from fit in that order
   * @return Double array with fit values of f, h, and residual
   */
  public double[] getFitParams() {
    return new double[]{fCorr, hCorr, fitResid};
  }
  
  /**
   * Returns the corner, damping, and residual from the initial parameters
   * calculated from the resp input as a double array in that order
   * @return Double array with the specified f, h, r values
   */
  public double[] getInitParams() {
    return new double[]{f, h, initResid};
  }
  
  /**
   * Return the residual values of initial and fit parameters
   * @return Array whose entries are, respectively, the initial and 
   * fit residuals
   */
  public double[] getResiduals() {
    return new double[]{initResid, fitResid};
  }

  @Override
  public boolean hasEnoughData(DataStore ds) {
    return ( ds.blockIsSet(0) && ds.bothComponentsSet(1) );
  }
  
  /**
   * Computes the forward change in value of the calculations for response
   * formed from a given corner and damping value
   * @param variables Vector with the corner and damping values from which the
   * derivatives are calculated
   * @return The result at the passed-in point plus the approximate derivative
   * of these points, as a vector and matrix respectively
   */
  private Pair<RealVector, RealMatrix> jacobian(RealVector variables) {
    
    // approximate through forward differences
    double[][] jacobian = new double[trimmedLength][2];
    
    double f1 = variables.getEntry(0);
    double h1 = variables.getEntry(1);
    double f2 = f1 + STEP_FACTOR;
    double h2 = h1 + STEP_FACTOR;
    
    double[] fInit = calculate(new double[]{f1, h1});
    double[] diffOnF = calculate(new double[]{f2, h1});
    double[] diffOnH = calculate(new double[]{f1, h2});

    for (int i = 0; i < trimmedLength; ++i) {
      jacobian[i][0] = (diffOnF[i] - fInit[i]) / (f2 - f1);
      jacobian[i][1] = (diffOnH[i] - fInit[i]) / (h2 - h1);
    }
    
    RealMatrix jMat = MatrixUtils.createRealMatrix(jacobian);
    RealVector fnc = MatrixUtils.createRealVector(fInit);
    
    return new Pair<RealVector, RealMatrix>(fnc, jMat);
  }
  
  @Override
  public int[] listActiveResponseIndices() {
    // NOTE: not used by corresponding panel, overrides with active indices
    // of components in the combo-box
    return new int[]{sensorOutIdx};
  }
  
}
