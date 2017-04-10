package asl.sensor.test;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.complex.Complex;
import org.junit.Test;

import asl.sensor.experiment.RandomizedExperiment;
import asl.sensor.input.InstrumentResponse;

public class RandomizedExperimentTest {

  @Test
  public void ResponseSetCorrectlyLowFreq() {
    
    String fname = "responses/TST5_response.txt";
    
    try {
      InstrumentResponse ir = new InstrumentResponse(fname);
      boolean lowFreq = true;
      List<Complex> poles = new ArrayList<Complex>( ir.getPoles() );
      double[] newPoles = new double[4];
      newPoles[0] = 0.;
      newPoles[1] = 1.;
      newPoles[2] = 0.;
      newPoles[3] = 1.;
      InstrumentResponse ir2 = 
          RandomizedExperiment.polesToResp(newPoles, ir, lowFreq);
      for (int i = 0; i < poles.size(); ++i) {
        Complex cmp, cmp2;
        if (i < 2) {
          cmp = ir2.getPoles().get(i);
          cmp2 = new Complex(0., 1.);
          
          assertFalse( cmp.equals( poles.get(i) ) );
          
        } else {
          cmp = ir2.getPoles().get(i);
          cmp2 = poles.get(i);
        }
        assertTrue( cmp.equals(cmp2) );
      }
    } catch (IOException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
      fail();
    }
    
  }
  
  @Test
  public void ResponseSetCorrectlyHighFreq() {
    String fname = "responses/TST5_response.txt";
    
    try {
      InstrumentResponse ir = new InstrumentResponse(fname);
      boolean lowFreq = false;
      List<Complex> poles = new ArrayList<Complex>( ir.getPoles() );
      double[] newPoles = new double[( poles.size() * 2 ) - 4];
      for (int i = 0; i < newPoles.length; ++i) {
        newPoles[i] = i;
      }
      InstrumentResponse ir2 = 
          RandomizedExperiment.polesToResp(newPoles, ir, lowFreq);
      for (int i = 0; i < poles.size(); ++i) {
        Complex cmp, cmp2;
        if (i < 2) {
          cmp = ir2.getPoles().get(i);
          cmp2 = poles.get(i);
        } else {
          int idx = (i - 2) * 2;
          cmp = ir2.getPoles().get(i);
          cmp2 = new Complex(idx, idx + 1);
          
          assertFalse( cmp.equals( poles.get(i) ) );
          
        }
        assertTrue( cmp.equals( cmp2 ) );
      }
    } catch (IOException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
      fail();
    }
  }
  
}
