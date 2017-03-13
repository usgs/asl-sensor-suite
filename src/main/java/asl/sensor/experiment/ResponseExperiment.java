package asl.sensor.experiment;

import org.apache.commons.math3.complex.Complex;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import asl.sensor.input.DataStore;
import asl.sensor.input.InstrumentResponse;

public class ResponseExperiment extends Experiment {

  @Override
  protected void backend(DataStore ds, final boolean freqSpace) {
    
    xySeriesData = new XYSeriesCollection();
    
    InstrumentResponse ir = ds.getResponse(0);
    double lowFreq = .0001;
    double highFreq = 10000;
    
    int pointCount = 100000;
    double freqDelta = (highFreq - lowFreq) / pointCount;
    
    double[] freqArray = new double[pointCount];
    
    double currentFreq = lowFreq;
    for (int i = 0; i < freqArray.length; ++i) {
      freqArray[i] = currentFreq;
      currentFreq += freqDelta;
    }
    
    Complex[] result = ir.applyResponseToInput(freqArray);
    XYSeries magnitude = new XYSeries("Response magnitude");
    for (int i = 0; i < freqArray.length; ++i) {
      if (freqSpace) {
        magnitude.add( freqArray[i], 10 * Math.log10( result[i].abs() ) );
      } else {
        magnitude.add( 1/freqArray[i], 10 * Math.log10( result[i].abs() ) );
      }
    }
    
    
    ((XYSeriesCollection) xySeriesData).addSeries(magnitude);
  }

  @Override
  public boolean hasEnoughData(DataStore ds) {
    return ds.responseIsSet(0);
  }

}
