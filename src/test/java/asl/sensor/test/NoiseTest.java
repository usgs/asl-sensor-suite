package asl.sensor.test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Calendar;

import org.jfree.data.xy.XYDataItem;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.junit.Test;

import asl.sensor.experiment.NoiseExperiment;
import asl.sensor.input.DataBlock;
import asl.sensor.input.DataStore;
import asl.sensor.input.InstrumentResponse;
import asl.sensor.utils.FFTResult;
import asl.sensor.utils.TimeSeriesUtils;

public class NoiseTest {

  XYSeriesCollection xysc;
  
  //@Test
  public void testResp() {
    String resp = "T-compact_Q330HR_BH_40";
    try {
      InstrumentResponse ir = InstrumentResponse.loadEmbeddedResponse(resp);
      System.out.println(ir.toString());
    } catch (IOException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
      fail();
    }
    
  }
  
  public XYSeriesCollection setUpTest1() throws FileNotFoundException {
    String folder = "data/noise_1/";
    String[] data = new String[3];
    data[0] = "00_BH0.512.seed";
    data[1] = "10_BH0.512.seed";
    data[2] = "TST6." + data[0];
    String resp = "T-compact_Q330HR_BH_40";
    DataStore ds = new DataStore();
    for (int i = 0; i < data.length; ++i) {
      DataBlock db = TimeSeriesUtils.getFirstTimeSeries(folder + data[i]);
      ds.setBlock(i, db);
      ds.setEmbedResponse(i, resp);
    }
    Calendar startCal = ds.getBlock(0).getStartCalendar();
    startCal.set(Calendar.HOUR_OF_DAY, 0);
    startCal.set(Calendar.MINUTE, 59);
    startCal.set(Calendar.SECOND, 59);
    startCal.set(Calendar.MILLISECOND, 994);
    Calendar endCal = (Calendar) startCal.clone();
    endCal.set(Calendar.HOUR_OF_DAY, 7);
    endCal.set(Calendar.MINUTE, 0);
    endCal.set(Calendar.SECOND, 0);
    endCal.set(Calendar.MILLISECOND, 25);
    ds.trim(startCal, endCal);
    NoiseExperiment ne = new NoiseExperiment();
    ne.setFreqSpace(false); // use period units (s)
    ne.runExperimentOnData(ds);
    XYSeriesCollection xysc = ne.getData().get(0);
    return xysc;
  }
  
  public XYSeriesCollection setUpTest2() throws FileNotFoundException {
    String folder = "data/noise_1Hz/";
    String[] data = new String[3];
    data[0] = "00_LH0.512.seed";
    data[1] = "10_LH0.512.seed";
    data[2] = "TST6." + data[0];
    String resp = "T-compact_Q330HR_BH_40";
    DataStore ds = new DataStore();
    for (int i = 0; i < data.length; ++i) {
      DataBlock db = TimeSeriesUtils.getFirstTimeSeries(folder + data[i]);
      ds.setBlock(i, db);
      ds.setEmbedResponse(i, resp);
    }
    Calendar startCal = ds.getBlock(0).getStartCalendar();
    startCal.set(Calendar.HOUR_OF_DAY, 0);
    startCal.set(Calendar.MINUTE, 59);
    startCal.set(Calendar.SECOND, 59);
    startCal.set(Calendar.MILLISECOND, 994);
    Calendar endCal = (Calendar) startCal.clone();
    endCal.set(Calendar.HOUR_OF_DAY, 7);
    endCal.set(Calendar.MINUTE, 0);
    endCal.set(Calendar.SECOND, 0);
    endCal.set(Calendar.MILLISECOND, 25);
    ds.trim(startCal, endCal);
    NoiseExperiment ne = new NoiseExperiment();
    ne.setFreqSpace(false); // use period units (s)
    ne.runExperimentOnData(ds);
    XYSeriesCollection xysc = ne.getData().get(0);
    return xysc;
  }
  
  @Test
  public void testResultsData1PSD1() {
    int idx = 0;
    double psdCheck = -154.63;
    double noiseCheck = -155.52;
    // everything below here same for every test
    try{
      XYSeriesCollection xysc = setUpTest1();
      // first 3 data, PSDs of each input
      // second 3 data, self-noise of each input
      // want data from 30 to 100s
      double low = 30.;
      double high = 100.;
      double psdResults = 0.;
      double noiseResults = 0.;
      XYSeries psd = xysc.getSeries(idx);
      XYSeries noise = xysc.getSeries(idx + 3);
      int psdPoints = 0;
      int noisePoints = 0;
      for (int j = 0; j < psd.getItemCount(); ++j) {
        XYDataItem psdxy = psd.getDataItem(j);
        double x = psdxy.getX().doubleValue();
        if (x >= low && x <= high) {
          psdResults += psdxy.getY().doubleValue();
          ++psdPoints;
        }
        XYDataItem noisxy = noise.getDataItem(j);
        x = noisxy.getX().doubleValue();
        if (x >= low && x <= high) {
          noiseResults += noisxy.getY().doubleValue();
          ++noisePoints;
        }
      }
      psdResults /= psdPoints;
      noiseResults /= noisePoints;
      /*
      System.out.println(psdResults + "," + noiseResults);
      System.out.println(psdCheck + "," + noiseCheck);
      System.out.println("PSD DIFF: " + Math.abs(psdResults - psdCheck));
      System.out.println("NOISE DIFF: " + Math.abs(noiseResults - noiseCheck));
      */
      assertEquals(noiseCheck, noiseResults, 1E-2);
      assertEquals(psdCheck, psdResults, 1E-2);
      
    } catch (FileNotFoundException e) {
      e.printStackTrace();
      fail();
    }
  }
  
  @Test
  public void testResultsData2PSD1() {
    int idx = 0;
    double psdCheck = -155.67;
    double noiseCheck = -156.62;
    // everything below here same for every test
    try{
      XYSeriesCollection xysc = setUpTest2();
      // first 3 data, PSDs of each input
      // second 3 data, self-noise of each input
      // want data from 30 to 100s
      double low = 30.;
      double high = 100.;
      double psdResults = 0.;
      double noiseResults = 0.;
      XYSeries psd = xysc.getSeries(idx);
      XYSeries noise = xysc.getSeries(idx + 3);
      int psdPoints = 0;
      int noisePoints = 0;
      for (int j = 0; j < psd.getItemCount(); ++j) {
        XYDataItem psdxy = psd.getDataItem(j);
        double x = psdxy.getX().doubleValue();
        if (x >= low && x <= high) {
          psdResults += psdxy.getY().doubleValue();
          ++psdPoints;
        }
        XYDataItem noisxy = noise.getDataItem(j);
        x = noisxy.getX().doubleValue();
        if (x >= low && x <= high) {
          noiseResults += noisxy.getY().doubleValue();
          ++noisePoints;
        }
      }
      psdResults /= psdPoints;
      noiseResults /= noisePoints;
      assertEquals(noiseCheck, noiseResults, 1E-2);
      assertEquals(psdCheck, psdResults, 1E-2);
      
    } catch (FileNotFoundException e) {
      e.printStackTrace();
      fail();
    }
  }
  
  @Test
  public void testResultsData1PSD2() {
    int idx = 1;
    double psdCheck = -156.04;
    double noiseCheck = -157.78;
    // everything below here same for every test
    try{
      XYSeriesCollection xysc = setUpTest1();
      // first 3 data, PSDs of each input
      // second 3 data, self-noise of each input
      // want data from 30 to 100s
      double low = 30.;
      double high = 100.;
      double psdResults = 0.;
      double noiseResults = 0.;
      XYSeries psd = xysc.getSeries(idx);
      XYSeries noise = xysc.getSeries(idx + 3);
      int psdPoints = 0;
      int noisePoints = 0;
      for (int j = 0; j < psd.getItemCount(); ++j) {
        XYDataItem psdxy = psd.getDataItem(j);
        double x = psdxy.getX().doubleValue();
        if (x >= low && x <= high) {
          psdResults += psdxy.getY().doubleValue();
          ++psdPoints;
        }
        XYDataItem noisxy = noise.getDataItem(j);
        x = noisxy.getX().doubleValue();
        if (x >= low && x <= high) {
          noiseResults += noisxy.getY().doubleValue();
          ++noisePoints;
        }
      }
      psdResults /= psdPoints;
      noiseResults /= noisePoints;
      assertEquals(noiseCheck, noiseResults, 1E-2);
      assertEquals(psdCheck, psdResults, 1E-2);
      
    } catch (FileNotFoundException e) {
      e.printStackTrace();
      fail();
    }
  }
 
  @Test
  public void testResultsData1PSD3() {
    int idx = 2;
    double psdCheck = -152.55;
    double noiseCheck = -153.32;
    // everything below here same for every test
    try{
      XYSeriesCollection xysc = setUpTest1();
      // first 3 data, PSDs of each input
      // second 3 data, self-noise of each input
      // want data from 30 to 100s
      double low = 30.;
      double high = 100.;
      double psdResults = 0.;
      double noiseResults = 0.;
      XYSeries psd = xysc.getSeries(idx);
      XYSeries noise = xysc.getSeries(idx + 3);
      int psdPoints = 0;
      int noisePoints = 0;
      for (int j = 0; j < psd.getItemCount(); ++j) {
        XYDataItem psdxy = psd.getDataItem(j);
        double x = psdxy.getX().doubleValue();
        if (x >= low && x <= high) {
          psdResults += psdxy.getY().doubleValue();
          ++psdPoints;
        }
        XYDataItem noisxy = noise.getDataItem(j);
        x = noisxy.getX().doubleValue();
        if (x >= low && x <= high) {
          noiseResults += noisxy.getY().doubleValue();
          ++noisePoints;
        }
      }
      psdResults /= psdPoints;
      noiseResults /= noisePoints;
      assertEquals(noiseCheck, noiseResults, 1E-2);
      assertEquals(psdCheck, psdResults, 1E-2);
      
    } catch (FileNotFoundException e) {
      e.printStackTrace();
      fail();
    }
  }
  
  //@Test
  public void anotherDamnPrintFunction() {
    // still beats feeding the inputs in by hand though
    String folder = "data/noise_1/";
    String[] data = new String[3];
    data[0] = "00_BH0.512.seed";
    data[1] = "10_BH0.512.seed";
    data[2] = "TST6." + data[0];
    String resp = "T-compact_Q330HR_BH_40";
    try {
      InstrumentResponse ir = InstrumentResponse.loadEmbeddedResponse(resp);
      for (int i = 0; i < data.length; ++i) {
        String name = folder + data[i];
        DataBlock db = TimeSeriesUtils.getFirstTimeSeries(name);
        /*
        Calendar startCal = db.getStartCalendar();
        startCal.set(Calendar.HOUR_OF_DAY, 0);
        startCal.set(Calendar.MINUTE, 59);
        startCal.set(Calendar.SECOND, 59);
        startCal.set(Calendar.MILLISECOND, 994);
        Calendar endCal = (Calendar) startCal.clone();
        endCal.set(Calendar.HOUR_OF_DAY, 7);
        endCal.set(Calendar.MINUTE, 0);
        endCal.set(Calendar.SECOND, 0);
        endCal.set(Calendar.MILLISECOND, 25);
        db.trim(startCal, endCal);
        */
        FFTResult fftr1 = FFTResult.spectralCalc(db, db);
        FFTResult fftr2 = FFTResult.crossPower(db, db, ir, ir);
        double[] freqs = fftr1.getFreqs();
        StringBuilder csv = new StringBuilder();
        for (int j = 0; j < freqs.length; ++j) {
          csv.append(freqs[j]);
          csv.append(", "); 
          csv.append(fftr1.getFFT(j).getReal());
          csv.append(", ");
          csv.append(fftr2.getFFT(j).getReal());
          if (j+1 < freqs.length) {
            csv.append("\n");
          }
        }
        String outPath = "testResultImages/PSD_OF_"+data[i]+".txt";
        PrintWriter out = new PrintWriter(outPath);
        out.write(csv.toString());
        out.close();
      }
    } catch (IOException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
      fail();
    }

  }
  
}
