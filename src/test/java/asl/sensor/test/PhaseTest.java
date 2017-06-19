package asl.sensor.test;

import static org.junit.Assert.*;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.List;
import java.util.TimeZone;

import org.apache.commons.math3.complex.Complex;
import org.junit.Test;

import asl.sensor.gui.InputPanel;
import asl.sensor.input.DataBlock;
import asl.sensor.input.DataStore;
import asl.sensor.input.InstrumentResponse;
import asl.sensor.utils.FFTResult;
import asl.sensor.utils.TimeSeriesUtils;

public class PhaseTest {

  public DataStore getFromList(List<String> setUpFilenames) throws IOException {

    String respName = setUpFilenames.get(0);
    String calName = setUpFilenames.get(1);
    String sensOutName = setUpFilenames.get(2);
    String metaName;

    System.out.println(respName);
    System.out.println(calName);
    System.out.println(sensOutName);

    metaName = TimeSeriesUtils.getMplexNameList(calName).get(0);
    DataBlock calib = TimeSeriesUtils.getTimeSeries(calName, metaName);

    InstrumentResponse ir = new InstrumentResponse(respName);

    metaName = TimeSeriesUtils.getMplexNameList(sensOutName).get(0);
    DataBlock sensor = TimeSeriesUtils.getTimeSeries(sensOutName, metaName);

    DataStore ds = new DataStore();
    ds.setData(0, calib);
    ds.setData(1, sensor);
    ds.setResponse(1, ir);

    return ds;

  }

  public Calendar getStartCalendar(DataStore ds) {
    SimpleDateFormat sdf = InputPanel.SDF;
    sdf.setTimeZone( TimeZone.getTimeZone("UTC") );
    Calendar cCal = Calendar.getInstance( sdf.getTimeZone() );

    cCal.setTimeInMillis( ds.getBlock(0).getStartTime() / 1000 );
    return cCal;
  }

  public DataStore setUpTest1() throws IOException {

    List<String> fileList = new ArrayList<String>();
    String respName = "responses/RESP.XX.NS088..BHZ.STS1.360.2400";
    String dataFolderName = "data/random_cal/"; 
    String calName =  dataFolderName + "_EC0.512.seed";
    String sensOutName = dataFolderName + "00_EHZ.512.seed";

    fileList.add(respName);
    fileList.add(calName);
    fileList.add(sensOutName);

    DataStore ds = getFromList(fileList);   
    Calendar cCal = getStartCalendar(ds);

    cCal.set(Calendar.MINUTE, 36);
    cCal.set(Calendar.SECOND, 0);
    long start = cCal.getTime().getTime() * 1000L;

    cCal.set(Calendar.MINUTE, 41);
    // System.out.println( "end: " + sdf.format( cCal.getTime() ) );
    long end = cCal.getTime().getTime() * 1000L;

    ds.trimAll(start, end);

    return ds;
  }

  @Test
  public void testFFTPhaseCalc() {

    DataStore ds = new DataStore();

    try {
      ds = setUpTest1();
    } catch (IOException e) {
      fail();
      e.printStackTrace();
      return;
    }
    
    DataBlock cal = ds.getBlock(0);
    DataBlock out = ds.getBlock(1);
    
    FFTResult calFFT = FFTResult.singleSidedFFT(cal, false);
    FFTResult outFFT = FFTResult.singleSidedFFT(out, false);
    
    Complex[] result = new Complex[outFFT.size()];
    double[] freqs = calFFT.getFreqs();
    
    StringBuilder sb = new StringBuilder();
    
    for (int i = 0; i < result.length; ++i) {
      Complex numer = outFFT.getFFT(i);
      Complex denom = calFFT.getFFT(i);
      Complex conj = calFFT.getFFT(i).conjugate();
      numer = numer.multiply(conj);
      denom = denom.multiply(conj);
      
      Complex entry = numer.divide(denom);
      result[i] = entry;
      sb.append(entry.getReal());
      sb.append(", ");
      sb.append(entry.getImaginary());
      sb.append(", ");
      sb.append(freqs[i]);
      sb.append("\n");
    }
    
    String folderName = "testResultImages";
    File folder = new File(folderName);
    if ( !folder.exists() ) {
      System.out.println("Writing directory " + folderName);
      folder.mkdirs();
    }

    String textName = folderName + "/outputData.txt";
    PrintWriter write;
    try {
      write = new PrintWriter(textName);
      write.println( sb.toString() );
      write.close();
    } catch (FileNotFoundException e) {
      fail();
      e.printStackTrace();
    }

    
  }

}
