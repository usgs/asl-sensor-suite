package asl.sensor.experiment;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;

import org.apache.commons.math3.complex.Complex;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import asl.sensor.input.DataStore;
import asl.sensor.input.InstrumentResponse;
import asl.sensor.utils.NumericUtils;

/**
 * Produces plots of response curves' magnitudes (Bode plot) and angle of
 * rotation in complex space. 1-3 responses can be plotted at a time.
 * No timeseries data (that is, miniSEED) is used in this calculation.
 * Response curves can be plotted in either frequency or interval space
 * (units of Hz or seconds respectively).
 * @author akearns
 *
 */
public class ResponseExperiment extends Experiment {

  public static final String MAGNITUDE = "Response amplitude";
  public static final String ARGUMENT = "Response phase";
  
  private boolean freqSpace;
  
  private Set<InstrumentResponse> responses;
  
  public ResponseExperiment() {
    super();
    freqSpace = false;
  }
  
  @Override
  protected void backend(DataStore ds) {
    
    responses = new HashSet<InstrumentResponse>();
    
    double lowFreq = .0001;
    double highFreq = 10000;
    
    int pointCount = 100000;
    double linearChange = (highFreq - lowFreq) / pointCount;
    // find logarithmic parameters for linear components
    double b = Math.log10(lowFreq / highFreq) / (lowFreq - highFreq);
    double a = lowFreq / Math.pow(10, b * lowFreq);
    
    double[] freqArray = new double[pointCount];
    
    double currentFreq = lowFreq;
    for (int i = 0; i < freqArray.length; ++i) {
      freqArray[i] = currentFreq;
      // System.out.println(currentFreq);
      currentFreq = a * Math.pow(10, b * (i * linearChange) );
    }
    
    // used to prevent issues with duplicate response plotting / XYSeries names
    Set<String> respNames = new HashSet<String>();
    
    XYSeriesCollection args = new XYSeriesCollection();
    XYSeriesCollection mags = new XYSeriesCollection();
    
    for (int r = 0; r < 3; ++r) {
      if ( !ds.responseIsSet(r) ) {
        continue;
      }
      
      InstrumentResponse ir = ds.getResponse(r);
      
      if ( respNames.contains( ir.getName() ) ) {
        continue;
      } else {
        respNames.add( ir.getName() );
        responses.add(ir);
      }
       
      Complex[] result = ir.applyResponseToInput(freqArray);
      
      String name = ir.getName();
      
      double phiPrev = 0; // use with unwrapping
      XYSeries magnitude = new XYSeries(name + " " + MAGNITUDE);
      XYSeries argument = new XYSeries (name + " " + ARGUMENT);
      for (int i = 0; i < freqArray.length; ++i) {
        Complex scaleFactor = new Complex(0., NumericUtils.TAU * freqArray[i]);
        Complex tmp = result[i].divide(scaleFactor);
        // Complex tmp = result[i].divide(NumericUtils.TAU * freqArray[i]);
        
        double phi = NumericUtils.atanc(tmp);
        phi = NumericUtils.unwrap(phi, phiPrev);
        phiPrev = phi;
        phi = Math.toDegrees(phi);
        double magAccel = tmp.abs();
        if (freqSpace) {
          magnitude.add( freqArray[i], 10 * Math.log10(magAccel) );
          argument.add( freqArray[i], phi );
        } else {
          magnitude.add( 1/freqArray[i], 10 * Math.log10(magAccel) );
          argument.add( 1/freqArray[i], phi );
        }
      }
      
      mags.addSeries(magnitude);
      args.addSeries(argument);
      
    }
    
    xySeriesData.add(mags);
    xySeriesData.add(args);

    dataNames = new ArrayList<String>(respNames); 
    
  }
  
  @Override
  public int blocksNeeded() {
    return 0;
  }

  @Override
  public long getEnd() {
    return 0L;
  }

  @Override
  public long getStart() {
    // there is no actual timeseries data used in this experiment
    // only response data, so start and end times are not defined
    return 0L;
  }
  
  @Override
  public boolean hasEnoughData(DataStore ds) {
    for (int i = 0; i < 3; ++i) {
      if ( ds.responseIsSet(i) ) {
        return true;
      }
    }
    
    return false;
  }
  
  public InstrumentResponse[] getResponses() {
    return responses.toArray( new InstrumentResponse[responses.size()] );
  }
  
  /**
   * Used to set the x-axis over which the response curve is plotted,
   * either frequency (Hz) units or sample-interval (s) units
   * @param freqSpace True if the plot should use units of Hz
   */
  public void setFreqSpace(boolean freqSpace) {
    this.freqSpace = freqSpace;
  }

}
