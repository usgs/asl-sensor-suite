package asl.sensor.test;

import static org.junit.Assert.*;

import java.awt.Font;
import java.io.File;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.List;
import java.util.TimeZone;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.linear.RealVector;
import org.apache.pdfbox.pdmodel.PDDocument;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.LogarithmicAxis;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.junit.Test;

import asl.sensor.experiment.ExperimentEnum;
import asl.sensor.experiment.ExperimentFactory;
import asl.sensor.experiment.RandomizedExperiment;
import asl.sensor.gui.InputPanel;
import asl.sensor.gui.RandomizedPanel;
import asl.sensor.input.DataBlock;
import asl.sensor.input.DataStore;
import asl.sensor.input.InstrumentResponse;
import asl.sensor.utils.ReportingUtils;
import asl.sensor.utils.TimeSeriesUtils;

public class RandomizedExperimentTest {

  
  @Test
  public void AutomateRandomCalTesting() {
    
    String currentDir = System.getProperty("user.dir");

    
    String respName = "responses/RESP.XX.NS088..BHZ.STS1.360.2400";
    String dataFolderName = "data/random_cal/"; 
    String calName =  dataFolderName + "_EC0.512.seed";
    String sensOutName = dataFolderName + "00_EHZ.512.seed";
    String metaName;
    
    boolean lowFreq = false;
    
    try {
      metaName = TimeSeriesUtils.getMplexNameList(calName).get(0);
      DataBlock calib = TimeSeriesUtils.getTimeSeries(calName, metaName);
      
      InstrumentResponse ir = new InstrumentResponse(respName);
      
      metaName = TimeSeriesUtils.getMplexNameList(sensOutName).get(0);
      DataBlock sensor = TimeSeriesUtils.getTimeSeries(sensOutName, metaName);

      DataStore ds = new DataStore();
      ds.setData(0, calib);
      ds.setData(1, sensor);
      ds.setResponse(1, ir);
      
      SimpleDateFormat sdf = InputPanel.SDF;
      sdf.setTimeZone( TimeZone.getTimeZone("UTC") );
      // sdf.setLenient(false);
      
      Calendar cCal = Calendar.getInstance( sdf.getTimeZone() );
      cCal.setTimeInMillis( ds.getBlock(0).getStartTime() / 1000 );
      
      // cCal.set(Calendar.DATE, 147);
      // cCal.set(Calendar.HOUR, 20);
      cCal.set(Calendar.MINUTE, 36);
      cCal.set(Calendar.SECOND, 0);
      // cCal.set(Calendar.MILLISECOND, 0);
      // System.out.println( "start: " + sdf.format( cCal.getTime() ) );
      long start = cCal.getTime().getTime() * 1000L;
      
      cCal.set(Calendar.MINUTE, 41);
      // System.out.println( "end: " + sdf.format( cCal.getTime() ) );
      long end = cCal.getTime().getTime() * 1000L;
      
      ds.trimAll(start, end);
      
      RandomizedExperiment rCal = (RandomizedExperiment)
          ExperimentFactory.createExperiment(ExperimentEnum.RANDM);
      
      rCal.setLowFreq(lowFreq);
      
      assertTrue( rCal.hasEnoughData(ds) );
      rCal.setData(ds);
      
      int width = 1280;
      int height = 960;
      
      // TODO: simplify most of this with calls to ReportingUtils functions
      
      List<XYSeriesCollection> xysc = rCal.getData();
      JFreeChart[] jfcl = new JFreeChart[xysc.size()];
      String[] yAxisTitles = new String[]{"Resp(f), dB", "Angle / TAU"};
      
      String xAxisTitle = "Frequency (Hz)";
      NumberAxis xAxis = new LogarithmicAxis(xAxisTitle);
      Font bold = xAxis.getLabelFont().deriveFont(Font.BOLD);
      xAxis.setLabelFont(bold);
      
      StringBuilder sb = new StringBuilder();
      sb.append( RandomizedPanel.getInsetString(rCal) );
      sb.append('\n');
      sb.append( RandomizedPanel.getTimeStampString(rCal) );
      sb.append('\n');
      sb.append("Response file used:");
      sb.append('\n');
      sb.append( rCal.getResponseName() );
      sb.append('\n');
      
      sb.append("BELOW RESULTS FOR EXPECTED BEST FIT (YELLOW CURVE)");
      sb.append('\n');
      
      // expected best fit params, for debugging
      double[] expectedParams = new double[]{-3.580104E+1, +7.122400E+1};
      ir = RandomizedExperiment.polesToResp(expectedParams, ir, lowFreq);
      ir.setName("Best-fit params");
      ds.setResponse(1, ir);
      rCal.setData(ds);

      // add initial curve from expected fit params to report
      XYSeries expectedInitialCurve = rCal.getData().get(0).getSeries(0);
      xysc.get(0).addSeries(expectedInitialCurve);
      XYSeries expectedInitialAngle = rCal.getData().get(1).getSeries(0);
      xysc.get(1).addSeries(expectedInitialAngle);
      
      for (int i = 0; i < xysc.size(); ++i) {
        
        jfcl[i] = ChartFactory.createXYLineChart(
            ExperimentEnum.RANDM.getName(),
            xAxisTitle,
            yAxisTitles[i],
            xysc.get(i),
            PlotOrientation.VERTICAL,
            true,
            false,
            false);
        
        XYPlot xyp = jfcl[i].getXYPlot();
                    
        //xyp.clearAnnotations();
        //xyp.addAnnotation(xyt);

        xyp.setDomainAxis( xAxis );
      }
      
      sb.append( RandomizedPanel.getInsetString(rCal) );
      
      PDDocument pdf = new PDDocument();
      ReportingUtils.chartsToPDFPage(width, height, pdf, jfcl);
      ReportingUtils.textToPDFPage( sb.toString(), pdf );
      
      String testResultFolder = currentDir + "/testResultImages/";
      File dir = new File(testResultFolder);
      if ( !dir.exists() ) {
        dir.mkdir();
      }
      
      String testResult = testResultFolder + "Random-Calib-Test.pdf";
      pdf.save( new File(testResult) );
      pdf.close();
      System.out.println("Output result has been written");
      
    } catch (IOException e) {
      e.printStackTrace();
      fail();
    }
  }
  
  @Test
  public void ResponseSetCorrectlyHighFreq() {
    
    String fname = "responses/TST5_response.txt";
    
    InstrumentResponse ir;

      try {
        ir = new InstrumentResponse(fname);
        boolean lowFreq = false;
        
        List<Complex> poles = new ArrayList<Complex>( ir.getPoles() );
        List<Complex> replacements = new ArrayList<Complex>();
        
        int start = 2;
        if ( poles.get(0).getImaginary() == 0 ) {
          start = 1;
        }
        
        for (int i = start; i < poles.size(); ++i) {
          if ( poles.get(i).getImaginary() == 0 ) {
            Complex c = new Complex(i, 0);
            replacements.add(c);
          } else {
            Complex c = new Complex(i, i * i);
            replacements.add(c);
            ++i;
          }
        }
        
        double[] newPoles = new double[replacements.size() * 2];
        for (int i = 0; i < newPoles.length; i += 2) {
          int poleIdx = i / 2;
          Complex c = replacements.get(poleIdx);
          newPoles[i] = c.getReal();
          newPoles[i + 1] = c.getImaginary();
        }
        
        InstrumentResponse ir2 = 
            RandomizedExperiment.polesToResp(newPoles, ir, lowFreq);
        
        List<Complex> testList = ir2.getPoles();
        int offsetIdx = 0;
        for (int i = 0; i < poles.size(); ++i) {
          if (i < start) {
            assertTrue( poles.get(i).equals( testList.get(i) ) );
          } else {
            Complex c = replacements.get(offsetIdx);
            assertTrue( testList.get(i).equals(c) );
            assertFalse( poles.get(i).equals(c) );
            if ( poles.get(i).getImaginary() != 0 ) {
              ++i;
              assertTrue( testList.get(i).equals( c.conjugate() ) );
              assertFalse( poles.get(i).equals( c.conjugate() ) );
            }
            ++offsetIdx;
          }
        }
        
      } catch (IOException e) {
        e.printStackTrace();
      }
    
  }
  
  @Test
  public void ResponseSetCorrectlyLowFreq() {
    
    String fname = "responses/TST5_response.txt";
    
      InstrumentResponse ir;
      try {
        ir = new InstrumentResponse(fname);
        boolean lowFreq = true;
        List<Complex> poles = new ArrayList<Complex>( ir.getPoles() );
        double[] newPoles = new double[2];
        newPoles[0] = 0.;
        newPoles[1] = 1.;
        
        Complex c = new Complex( newPoles[0], newPoles[1] );
        
        
        InstrumentResponse ir2 = 
            RandomizedExperiment.polesToResp(newPoles, ir, lowFreq);
        List<Complex> poles2 = ir2.getPoles();
        List<Complex> testList = new ArrayList<Complex>(poles);
        testList.set(0, c);
        testList.set( 1, c.conjugate() );
        
        // System.out.println(testList);
        
        for (int i = 0; i < poles.size(); ++i) {
          if (i < 2) {
            assertFalse( poles.get(i).equals( poles2.get(i) ) );
            assertTrue( poles2.get(i).equals( testList.get(i) ) );
          }
        }
        
        
      } catch (IOException e) {
        fail();
        e.printStackTrace();
      }
    
  }
  
  @Test
  public void ResponseCorrectlyConvertedToVectorLowFreq() {
    String fname = "responses/TST5_response.txt";
    
    InstrumentResponse ir;
    try {
      
      ir = new InstrumentResponse(fname);
      List<Complex> poles = new ArrayList<Complex>( ir.getPoles() );
      RealVector low = RandomizedExperiment.lowFreqPolesToVector(poles);
      
      // only test lower two poles
      
      assertEquals( low.getEntry(0), poles.get(0).getReal(), 0.0 );
      assertEquals( low.getEntry(1), poles.get(0).getImaginary(), 0.0 );
      
      assertEquals( low.getEntry(0), poles.get(1).getReal(), 0.0 );
      assertEquals( low.getEntry(1), -poles.get(1).getImaginary(), 0.0 );
      
    } catch (IOException e) {
      fail();
      e.printStackTrace();
    }
  }
  
  @Test
  public void ResponseCorrectConvertedToVectorHighFreq() {
    String fname = "responses/TST5_response.txt";
    
    InstrumentResponse ir;
    try {
      
      ir = new InstrumentResponse(fname);
      List<Complex> poles = new ArrayList<Complex>( ir.getPoles() );
      RealVector high = RandomizedExperiment.highFreqPolesToVector(poles);
      
      int complexIndex = 2; // start at second pole
      int vectorIndex = 0;
      
      while ( vectorIndex < high.getDimension() ) {
        // return current index
        double real = high.getEntry(vectorIndex++);
        double imag = high.getEntry(vectorIndex++);
        
        double poleImag = poles.get(complexIndex).getImaginary();
        
        assertEquals( real, poles.get(complexIndex).getReal(), 0.0 );
        assertEquals( imag, poleImag, 0.0 );
        
        if (poleImag != 0) {
          // complex conjugate case
          ++complexIndex;
          assertEquals( real, poles.get(complexIndex).getReal(), 0.0 );
          assertEquals( imag, -poles.get(complexIndex).getImaginary(), 0.0 );
        }
        
        ++complexIndex;
        
      }
      
    } catch (IOException e) {
      fail();
      e.printStackTrace();
    }
  }
}
