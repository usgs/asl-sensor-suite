package asl.sensor.test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.math3.complex.Complex;
import org.apache.pdfbox.pdmodel.PDDocument;
import org.junit.Test;

import asl.sensor.input.InstrumentResponse;
import asl.sensor.input.TransferFunction;
import asl.sensor.input.Unit;
import asl.sensor.utils.ReportingUtils;

public class InstrumentResponseTest {

  @Test
  public void testFileParse() {
    String filename = "./responses/RESP.XX.NS087..BHZ.STS1.20.2400";
    
    try {
      InstrumentResponse ir = new InstrumentResponse(filename);
      
      assertEquals(TransferFunction.LAPLACIAN, ir.getTransferFunction() );
      
      double nml = Double.parseDouble("3.948580E+03");
      assertEquals( nml, ir.getNormalization(), 0.0001 );
      
      double nmf = Double.parseDouble("3.000000E-01");
      assertEquals( nmf, ir.getNormalizationFrequency(), 0.0001 );
      
      Double[] gn = {2.400000e+03, 2.400000e+03, 1.000000e+00};
      List<Double> gnL = Arrays.asList(gn);
      assertTrue( gnL.equals(ir.getGain() ) );
      
      assertEquals( Unit.VELOCITY, ir.getUnits() );
      
      List<Complex> zrs = new ArrayList<Complex>();
      zrs.add( new Complex(0.000000e+00, 0.000000e+00) );
      zrs.add( new Complex(0.000000e+00, 0.000000e+00) );
      assertEquals( zrs, ir.getZeros() );
      
      List<Complex> pls = new ArrayList<Complex>();
      pls.add( new Complex(-2.221000e-01,  2.221000e-01) );
      pls.add( new Complex(-2.221000e-01, -2.221000e-01) );
      pls.add( new Complex(-3.918000e+01,  4.912000e+01) );
      pls.add( new Complex(-3.918000e+01, -4.912000e+01) );
      assertEquals( pls, ir.getPoles() );
      
    } catch (IOException e) {
      // TODO Auto-generated catch block
      fail("Unexpected error trying to read response file");
    }
    
  }
  
  @Test
  public void testStringOutput() {
    
    String currentDir = System.getProperty("user.dir");
    String filename = currentDir + "/responses/RESP.XX.NS087..BHZ.STS1.20.2400";

    try {
      InstrumentResponse ir = new InstrumentResponse(filename);
      
      System.out.println(ir);
      
      PDDocument pdf = new PDDocument();
      
      ReportingUtils.textToPDFPage( ir.toString(), pdf );
      
      String testResultFolder = currentDir + "/testResultImages/";
      File dir = new File(testResultFolder);
      if ( !dir.exists() ) {
        dir.mkdir();
      }
      
      String testResult = testResultFolder + "response-report.pdf";
      pdf.save( new File(testResult) );
      pdf.close();

    } catch (IOException e) {
      // TODO Auto-generated catch block
      fail("Unexpected error trying to read response file");
    }
  }

}
