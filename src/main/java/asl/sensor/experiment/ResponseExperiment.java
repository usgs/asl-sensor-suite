package asl.sensor.experiment;

import java.util.HashSet;
import java.util.Set;

import org.apache.commons.math3.complex.Complex;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import asl.sensor.input.DataStore;
import asl.sensor.input.InstrumentResponse;

public class ResponseExperiment extends Experiment {

  public static final String MAGNITUDE = "Response magnitude";
  public static final String ARGUMENT = "Response argument (phi)";
  
  private boolean freqSpace;
  
  public ResponseExperiment() {
    super();
    freqSpace = false;
  }
  
  public void setFreqSpace(boolean freqSpace) {
    this.freqSpace = freqSpace;
  }
  
  @Override
  protected void backend(DataStore ds) {
    
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
    
    Set<String> respNames = new HashSet<String>();
    xySeriesData = new XYSeriesCollection();
    
    for (int r = 0; r < 3; ++r) {
      if ( !ds.responseIsSet(r) ) {
        continue;
      }
      
      InstrumentResponse ir = ds.getResponse(r);
      
      if ( respNames.contains( ir.getName() ) ) {
        continue;
      } else {
        respNames.add( ir.getName() );
      }
       
      Complex[] result = ir.applyResponseToInput(freqArray);
      
      String name = ir.getName();
      
      XYSeries magnitude = new XYSeries(name + " " + MAGNITUDE);
      XYSeries argument = new XYSeries (name + " " + ARGUMENT);
      for (int i = 0; i < freqArray.length; ++i) {
        double phi = ( Math.toDegrees( result[i].getArgument() ) + 360 ) % 360;
        if (freqSpace) {
          double magAccel = result[i].abs() / (2*Math.PI*freqArray[i]);
          magnitude.add( freqArray[i], 10 * Math.log10(magAccel) );
          argument.add( freqArray[i], phi );
        } else {
          magnitude.add( 1/freqArray[i], 10 * Math.log10( result[i].abs() ) );
          argument.add( 1/freqArray[i], phi );
        }
      }
      
      xySeriesData.addSeries(magnitude);
      xySeriesData.addSeries(argument);
      
    }

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

}
