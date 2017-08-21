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
    ds.setBlock(0, calib);
    ds.setBlock(1, sensor);
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
    long end = cCal.getTimeInMillis() * 1000L;

    ds.trim(start, end);

    return ds;
  }
  
  public DataStore setUpTest3() throws IOException {
    List<String> fileList = new ArrayList<String>();
    String respName = "responses/RESP.XX.NS088..BHZ.STS1.360.2400";
    String dataFolderName = "data/random_cal_3/"; 
    String calName =  dataFolderName + "BC0.512.seed";
    String sensOutName = dataFolderName + "00_BHZ.512.seed";
    
    fileList.add(respName);
    fileList.add(calName);
    fileList.add(sensOutName);
    
    DataStore ds = getFromList(fileList);
    
    // response we want is embedded
    InstrumentResponse ir;
    ir = InstrumentResponse.loadEmbeddedResponse("KS54000_Q330HR_BH_40");
    ds.setResponse(1, ir);
    
    Calendar cCal = getStartCalendar(ds);
    
    cCal.set(Calendar.HOUR_OF_DAY, 21);
    cCal.set(Calendar.MINUTE, 24);
    cCal.set(Calendar.SECOND, 0);
    long start = cCal.getTime().getTime() * 1000L;
    
    // commented out -- calibration ends when the data does
    //int hour = cCal.get(Calendar.HOUR);
    /*
    cCal.set(Calendar.DAY_OF_YEAR, 4);
    cCal.set(Calendar.HOUR_OF_DAY, 0);
    cCal.set(Calendar.MINUTE, 0);
    */
    
    // System.out.println( "end: " + sdf.format( cCal.getTime() ) );
    long end = ds.getBlock(0).getEndTime();
    
    ds.trim(start, end);
    
    return ds;
  }
  
public DataStore setUpTest4() throws IOException {
    
    List<String> fileList = new ArrayList<String>();
    String respName = "responses/RESP.XX.NS088..BHZ.STS1.360.2400";
    String dataFolderName = "data/random_cal_4/"; 
    String calName =  dataFolderName + "CB_BC0.512.seed";
    String sensOutName = dataFolderName + "00_EHZ.512.seed";
    
    fileList.add(respName);
    fileList.add(calName);
    fileList.add(sensOutName);
    
    DataStore ds = getFromList(fileList);
    
    // response we want is embedded
    InstrumentResponse ir;
    ir = InstrumentResponse.loadEmbeddedResponse("KS54000_Q330HR_BH_40");
    ds.setResponse(1, ir);
    
    Calendar cCal = getStartCalendar(ds);
    
    cCal.set(Calendar.HOUR_OF_DAY, 20);
    cCal.set(Calendar.MINUTE, 16);
    cCal.set(Calendar.SECOND, 0);
    long start = cCal.getTime().getTime() * 1000L;
    
    cCal.set(Calendar.MINUTE, 26);
    // System.out.println( "end: " + sdf.format( cCal.getTime() ) );
    long end = cCal.getTime().getTime() * 1000L;
    
    ds.trim(start, end);
    
    return ds;
  }

  @Test
  public void testFFTDriver() {
    int[] tests = new int[]{1,3,4};
    // yes, folks, that's right, the abhorrent for-case paradigm
    for (int test : tests) {
      testFFTPhaseCalc(test);
    }
  }
  
  public void testFFTPhaseCalc(int testIteration) {

    DataStore ds = new DataStore();

    int test = testIteration;
    
    try {
      switch (test) {
      case 1:
        ds = setUpTest1();
        break;
      case 3:
        ds = setUpTest3();
        break;
      case 4:
        ds = setUpTest4();
        break;
        default:
          fail();
          return;
      }

    } catch (IOException e) {
      fail();
      e.printStackTrace();
      return;
    }
    
    long start = ds.getBlock(0).getStartTime();
    long interval = ds.getBlock(0).getInterval();
    
    DataBlock cal = ds.getBlock(0);
    DataBlock out = ds.getBlock(1);
    
    if ( interval != out.getInterval() ) {
      fail();
    }
    
    double[] calArr = new double[cal.size()];
    double[] outArr = new double[cal.size()];
    
    StringBuilder sbCalData = new StringBuilder();
    StringBuilder sbOutData = new StringBuilder();
    
    long nowTime = start;
    
    for (int i = 0; i < calArr.length; ++i) {
      calArr[i] = cal.getData().get(i).doubleValue();
      outArr[i] = out.getData().get(i).doubleValue();
      sbCalData.append(calArr[i]);
      sbCalData.append(", ");
      sbCalData.append(nowTime);
      sbOutData.append(outArr[i]);
      sbOutData.append(", ");
      sbOutData.append(nowTime);
      nowTime += interval;
      if ( i < (calArr.length - 1) ) {
        sbCalData.append("\n");
        sbOutData.append("\n");
      }
    }
    
    //double[] fq = FFTResult.singleSidedFFT(cal, false).getFreqs();
    //System.out.println( fq[fq.length - 1] );
    
    Complex[] calFFTCpx = FFTResult.simpleFFT(calArr);
    Complex[] outFFTCpx = FFTResult.simpleFFT(outArr);
    
    double[] freqs = new double[calFFTCpx.length];
    int numIdx = freqs.length / 2 + 1;
    double deltaFrq = cal.getSampleRate() / (numIdx * 2);
    for (int i = 1; i <= freqs.length / 2; ++i) {
      freqs[i] = i * deltaFrq;
      freqs[freqs.length - i] = - i * deltaFrq;
    }
    freqs[0] = 0;
    freqs[numIdx] = cal.getSampleRate() / 2;
    
    
    StringBuilder sbNumer = new StringBuilder();
    StringBuilder sbDenom = new StringBuilder();
    
    StringBuilder sbCalFFT = new StringBuilder();
    StringBuilder sbOutFFT = new StringBuilder();
    
    for (int i = 0; i < calFFTCpx.length; ++i) {
      Complex numer = outFFTCpx[i];
      Complex denom = calFFTCpx[i];
      Complex conj = calFFTCpx[i].conjugate();
      numer = numer.multiply(conj);
      denom = denom.multiply(conj);
      
      Complex entry = numer.divide(denom);
      
      // reset numerator in order to output it, without mult by conj
      // numer = outFFT.getFFT(i);
      // denom = calFFT.getFFT(i);
      
      sbCalFFT.append(calFFTCpx[i].getReal());
      sbCalFFT.append(", ");
      sbCalFFT.append(calFFTCpx[i].getImaginary());
      sbCalFFT.append(", ");
      sbCalFFT.append(freqs[i]);
      
      sbOutFFT.append(outFFTCpx[i].getReal());
      sbOutFFT.append(", ");
      sbOutFFT.append(outFFTCpx[i].getImaginary());
      sbOutFFT.append(", ");
      sbOutFFT.append(freqs[i]);
      
      sbNumer.append(numer.getReal());
      sbNumer.append(", ");
      sbNumer.append(numer.getImaginary());
      sbNumer.append(", ");
      sbNumer.append(freqs[i]);
      
      sbDenom.append(denom.getReal());
      sbDenom.append(", ");
      sbDenom.append(denom.getImaginary());
      sbDenom.append(", ");
      sbDenom.append(freqs[i]);
      
      if ( i < (calFFTCpx.length - 1) ) {
        // don't append final newline on last entry
        sbOutFFT.append("\n");
        sbCalFFT.append("\n");
        sbNumer.append("\n");
        sbDenom.append("\n");
      }

    }
    
    String folderName = "testResultImages";
    File folder = new File(folderName);
    if ( !folder.exists() ) {
      System.out.println("Writing directory " + folderName);
      folder.mkdirs();
    }

    String textNameNumer = folderName + "/outputData"+test+"_Numer.txt";
    String textNameDenom = folderName + "/outputData"+test+"_Denom.txt";
    String textNameCalIn = folderName + "/outputData"+test+"_CalSignal.txt";
    String textNameOutIn = folderName + "/outputData"+test+"_OutSignal.txt";
    String textNameCalFFT = folderName + "/outputData"+test+"_CalFFT.txt";
    String textNameOutFFT = folderName + "/outputData"+test+"_OutFFT.txt";
    
    PrintWriter write;
    try {
      write = new PrintWriter(textNameNumer);
      write.println( sbNumer.toString() );
      write.close();
      write = new PrintWriter(textNameDenom);
      write.println( sbDenom.toString() );
      write.close();
      write = new PrintWriter(textNameCalIn);
      write.println( sbCalData.toString() );
      write.close();
      write = new PrintWriter(textNameOutIn);
      write.println( sbOutData.toString() );
      write.close();
      write = new PrintWriter(textNameCalFFT);
      write.println( sbCalFFT.toString() );
      write.close();
      write = new PrintWriter(textNameOutFFT);
      write.println( sbOutFFT.toString() );
      write.close();
    } catch (FileNotFoundException e) {
      fail();
      e.printStackTrace();
    }

    
  }
  
  @Test
  public void quickNumericTest() {
    // get a number as integer
    Integer x = new Integer(123456);

    // get a number as float
    Float y = new Float(9876f);
    
    // get a number as double
    Double u = new Double(12304566983.045);

    // print their value as double
    System.out.println("x as integer: " + x
       + ", x as double: " + x.doubleValue());
    System.out.println("y as float: " + y
       + ", y as double: " + y.doubleValue());
    System.out.println("u as double: " + u
        + ", u as double: " + u.doubleValue());
  }

}
