package asl.sensor.test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.io.BufferedInputStream;
import java.io.DataInput;
import java.io.DataInputStream;
import java.io.EOFException;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math3.util.Pair;
import org.junit.Test;

import asl.sensor.input.DataBlock;
import asl.sensor.utils.TimeSeriesUtils;
import edu.iris.dmc.seedcodec.B1000Types;
import edu.iris.dmc.seedcodec.CodecException;
import edu.iris.dmc.seedcodec.DecompressedData;
import edu.iris.dmc.seedcodec.UnsupportedCompressionType;
import edu.sc.seis.seisFile.mseed.DataHeader;
import edu.sc.seis.seisFile.mseed.DataRecord;
import edu.sc.seis.seisFile.mseed.SeedFormatException;
import edu.sc.seis.seisFile.mseed.SeedRecord;

public class TimeSeriesUtilsTest {

  public String station = "TST5";
  public String location = "00";
  public String channel = "BH0";
  
  public String fileID = station+"_"+location+"_"+channel;
  
  public String filename1 = "./data/"+fileID+".512.seed";
  
  @Test
  public void canGetFile() {
    try{
      FileInputStream fis = new FileInputStream(filename1);
      fis.close();
    } catch (Exception e) {
      assertNull(e);
    }
  }
  
  
  @Test
  public void canGetMultiplexDataNames() {
    String filename2 = "./data/cat.seed";
    Set<String> names;
    
    try {
      names = TimeSeriesUtils.getMplexNameSet(filename2);
      assertTrue(names.contains("IU_ANMO_00_LH1"));
      assertTrue(names.contains("IU_ANMO_00_LH2"));
      assertTrue(names.contains("IU_ANMO_00_LHZ"));
      assertEquals( names.size(), 3 );
    } catch (FileNotFoundException e) {
      fail();
    }
  }
  
  @Test
  public void decimationTest() {
    
    long interval40Hz = (TimeSeriesUtils.ONE_HZ_INTERVAL / 40);
    long interval = TimeSeriesUtils.ONE_HZ_INTERVAL;
    
    List<Number> timeSeries = new ArrayList<Number>();
    
    for (int i = 0; i < 160; ++i) {
      timeSeries.add(i);
    }
    
    // System.out.println(timeSeries);
    
    timeSeries = TimeSeriesUtils.decimate(timeSeries, interval40Hz, interval);
    
    // System.out.println(timeSeries);
    
    assertEquals(timeSeries.size(), 4);
    for (int i = 0; i < timeSeries.size(); ++i) {
      assertEquals(timeSeries.get(i).doubleValue(), 40. * i, 0.5);
    }
  }
  
  @Test
  public void inputFileReaderCreatesXYSeries() {
    DataInput dis;
    DataBlock db = null;
    List<Number> data = new ArrayList<Number>();
    
    try {
      dis = new DataInputStream( 
            new BufferedInputStream( 
            new FileInputStream(filename1) ) );
      while ( true ) {

        try {
          long interval = 0L;
          SeedRecord sr = SeedRecord.read(dis,4096);
          if(sr instanceof DataRecord) {
            DataRecord dr = (DataRecord)sr;
            DataHeader dh = dr.getHeader();
            if (db == null){
              StringBuilder fileID = new StringBuilder();
              fileID.append(dh.getStationIdentifier() + "_");
              fileID.append(dh.getLocationIdentifier() + "_");
              fileID.append(dh.getChannelIdentifier());
              db = new DataBlock(data, interval, fileID.toString(), -1);
            }


            byte af = dh.getActivityFlags();
            byte correctionFlag = 0b00000010; // is there a time correction?
            int correction = 0;
            if ( (af & correctionFlag) == 0 ) {
              correction = dh.getTimeCorrection();
            }

            long start = dh.getStartBtime()
                .convertToCalendar()
                .getTimeInMillis() + correction;

            if(db.getStartTime() < 0) {
              db.setStartTime(start);
            }

            int fact = dh.getSampleRateFactor();
            int mult = dh.getSampleRateMultiplier();

            final long ONE_HZ_INTERVAL = TimeSeriesUtils.ONE_HZ_INTERVAL;
            
            if( fact > 0 && mult > 0) {
              interval = ONE_HZ_INTERVAL / (fact * mult);
            } else if (fact > 0 && mult < 0) {
              interval = Math.abs( (ONE_HZ_INTERVAL * mult) / fact);
            } else if (fact < 0 && mult > 0) {
              interval = Math.abs( (ONE_HZ_INTERVAL * fact) / mult);
            } else {
              interval = ONE_HZ_INTERVAL * fact * mult;
            }

            db.setInterval(interval);

            DecompressedData decomp = dr.decompress();

            // get the original datatype of the series (loads data faster)
            // otherwise the decompressed data gets converted (cloned) as
            // the other type instead
            int dataType = decomp.getType();

            // This is probably the best way to do this since
            // we have to add each point individually anyway
            // and converting between types for 



            switch (dataType) {
            case B1000Types.INTEGER:
              int[] decomArrayInt = decomp.getAsInt();
              for (int dataPoint : decomArrayInt ) {
                data.add(dataPoint);
              }
              break;
            case B1000Types.FLOAT:
              float[] decomArrayFlt = decomp.getAsFloat();
              for (float dataPoint : decomArrayFlt ) {
                data.add(dataPoint);
              }
              break;
            case B1000Types.SHORT:
              short[] decomArrayShr = decomp.getAsShort();
              for (short dataPoint : decomArrayShr ) {
                data.add(dataPoint);
              }
              break;
            default:
              double[] decomArrayDbl = decomp.getAsDouble();
              for (double dataPoint : decomArrayDbl ) {
                data.add(dataPoint);
              }
              break;
            }
          }
        } catch(EOFException e) {
          break;
        }
        
      }
      
      // quickly get the one name in the list
      Set<String> names = TimeSeriesUtils.getMplexNameSet(filename1);
      List<String> nameList = new ArrayList<String>(names);
      
      DataBlock testAgainst = 
          TimeSeriesUtils.getTimeSeries(filename1, nameList.get(0) );
      assertEquals( db.getData().size(), testAgainst.getData().size() );
      
    } catch (FileNotFoundException e) {
      assertNull(e);
    } catch (SeedFormatException e) {
      assertNull(e);
    } catch (IOException e) {
      assertNull(e);
    } catch (UnsupportedCompressionType e) {
      // TODO Auto-generated catch block
      assertNull(e);
    } catch (CodecException e) {
      // TODO Auto-generated catch block
      assertNull(e);
    }
  }
  
  @Test
  public void seisFileCanParseFile() {
    
    try {
      DataInput dis = new DataInputStream( new BufferedInputStream( 
          new FileInputStream(filename1) ) ); 
      try{
        while(true) {
          SeedRecord sr = SeedRecord.read(dis,4096);
          if (sr instanceof DataRecord) {
            DataRecord dr = (DataRecord)sr;

            String loc = dr.getHeader().getLocationIdentifier();
            assertTrue( loc.equals(location) );
            String stat = dr.getHeader().getStationIdentifier().trim();
            assertTrue( stat.equals(station) );

            String chan = dr.getHeader().getChannelIdentifier();
            assertTrue( chan.equals(channel) );
          }
        }
      } catch (EOFException e) {
        assertNotNull(e); // I haaates it! I haaaaaaaaaates it!
      } catch (SeedFormatException e) {
        assertNull(e);
      } catch (IOException e) {
        assertNull(e);
      }
      
    } catch (FileNotFoundException e) {
      // TODO Auto-generated catch block
      assertNull(e);
    }
    
  }
  
  @Test
  public void seisFileGivesCorrectSampleRateAndInterval() {
    DataInput dis;
    try {
      while (true) {
        
        dis = new DataInputStream( new BufferedInputStream( 
            new FileInputStream(filename1) ) );
        SeedRecord sr = SeedRecord.read(dis,4096);
        if(sr instanceof DataRecord) {
          DataRecord dr = (DataRecord)sr;
          
          int fact = dr.getHeader().getSampleRateFactor();
          int mult = dr.getHeader().getSampleRateMultiplier();
          
          //System.out.println(fact+","+mult);
          
          double rate = dr.getHeader().getSampleRate();
          assertTrue((double)fact/mult == rate);
          
          // checking the correct values for the intervals
          
          double multOf1Hz = rate/TimeSeriesUtils.ONE_HZ;
          long inverse = TimeSeriesUtils.ONE_HZ_INTERVAL/(long)multOf1Hz;
          
          long interval = TimeSeriesUtils.ONE_HZ_INTERVAL*mult/fact;
          
          assertEquals( inverse, interval);
          // System.out.println(interval);
          
          break;
          
        }
      }
    } catch (FileNotFoundException e) {
      assertNull(e); // only reading one record;
    } catch (SeedFormatException e) {
      assertNull(e);
    } catch (IOException e) {
      assertNull(e);
    }
  }
  
  @Test
  public void testInputParsing() {
    String dataFolderName = "data/random_cal_4/"; 
    String extension = "CB_BC0.512.seed";;
    doInputParseTest(dataFolderName, extension);
    extension = "00_EHZ.512.seed";
    doInputParseTest(dataFolderName, extension);
  }
  

  public void doInputParseTest(String dataFolderName, String extension) {

    String fileName =  dataFolderName + extension;
    String metaName;
    try {
      metaName = TimeSeriesUtils.getMplexNameList(fileName).get(0);
      Pair<Long, Map<Long, Number>> data = 
          TimeSeriesUtils.getTimeSeriesMap(fileName, metaName);
      Map<Long, Number> map = data.getSecond();
      
      List<Long> times = new ArrayList<Long>( map.keySet() );
      Collections.sort(times);
      long lastTime = times.get( times.size() - 1 );
      
      StringBuilder sb = new StringBuilder();
      for (Long time : times) {
        sb.append(time);
        sb.append(", ");
        sb.append( map.get(time) );
        if (time < lastTime) {
          sb.append("\n");
        }
      }
      
      String folderName = "testResultImages";
      File folder = new File(folderName);
      if ( !folder.exists() ) {
        System.out.println("Writing directory " + folderName);
        folder.mkdirs();
      }
      
      String outputFilename = folderName + "/"+extension+"_timeDataMap.txt";
      PrintWriter write;
      write = new PrintWriter(outputFilename);
      write.println( sb.toString() );
      write.close();
      
    } catch (FileNotFoundException e) {
      fail();
      e.printStackTrace();
    }
  }
  
}
