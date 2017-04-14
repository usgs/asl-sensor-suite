package asl.sensor.test;

import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.List;
import java.util.TimeZone;

import javax.imageio.ImageIO;

import org.junit.Test;

import asl.sensor.experiment.ExperimentEnum;
import asl.sensor.gui.ExperimentPanel;
import asl.sensor.gui.ExperimentPanelFactory;
import asl.sensor.gui.InputPanel;
import asl.sensor.input.DataStore;
import asl.sensor.utils.TimeSeriesUtils;

public class NoiseNineTest {

  String currentDir = System.getProperty("user.dir");

  String folder = currentDir + "/data/noisenine/";
  String[] types = new String[]{"00","10","60"};
  String freqName = "_BH";
  String[] components = new String[]{"1","2","Z"};
  String ending = ".512.seed";
  
  String respName = "responses/RESP.XX.NS088..BHZ.STS1.360.2400";
  
  @Test
  public void canRunAndPlotTest(){
    
    DataStore ds = new DataStore();
    
    for (int i = 0; i < types.length; ++i) {
      for (int j = 0; j < components.length; ++j) {
        int indexInStore = i * types.length + j;
        
        String fName = folder + types[i] + freqName + components[j] + ending;
        
        System.out.println(fName);
        
        List<String> dataNames = new ArrayList<String>();
        try {
          dataNames = new ArrayList<String>( 
              TimeSeriesUtils.getMplexNameSet(fName) );
        } catch (FileNotFoundException e) {
          // TODO Auto-generated catch block
          fail();
          e.printStackTrace();
        }
        ds.setData(indexInStore, fName, dataNames.get(0) );
        ds.setResponse(indexInStore, respName);
      }
    }
    
    SimpleDateFormat sdf = InputPanel.SDF;
    sdf.setTimeZone( TimeZone.getTimeZone("UTC") );
    // sdf.setLenient(false);
    
    Calendar cCal = Calendar.getInstance( sdf.getTimeZone() );
    cCal.setTimeInMillis( ds.getBlock(0).getStartTime() / 1000 );
    cCal.set(Calendar.HOUR, 7);
    System.out.println( "start: " + sdf.format( cCal.getTime() ) );
    long start = cCal.getTime().getTime() * 1000L;
    cCal.set(Calendar.HOUR, 8);
    cCal.set(Calendar.MINUTE, 30);
    System.out.println( "end: " + sdf.format( cCal.getTime() ) );
    long end = cCal.getTime().getTime() * 1000L;
    
    ds.trimAll(start, end);
    
    ExperimentPanel nnp = 
        ExperimentPanelFactory.createPanel(ExperimentEnum.NOIS9);
    assertTrue( nnp.hasEnoughData(ds) );
    nnp.updateData(ds);
    
    BufferedImage bi = nnp.getAsImage(640,480 * 3);
    try {
      String testResultFolder = currentDir + "/testResultImages/";
      File dir = new File(testResultFolder);
      if ( !dir.exists() ) {
        dir.mkdir();
      }
      String testResult = testResultFolder + "Nine-Noise-Test.png";
      ImageIO.write(bi,"png", new File(testResult) );
      System.out.println("Image has been written");
    } catch (IOException e) {
      // TODO Auto-generated catch block
      fail();
      e.printStackTrace();
    }
    
  }
  
}
