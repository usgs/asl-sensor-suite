package asl.sensor.test;

import static org.junit.Assert.*;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.net.HttpURLConnection;
import java.net.MalformedURLException;
import java.net.URL;

import org.apache.commons.io.FileUtils;
import org.junit.Test;

public class TestUtils {
  
  private static String URL_HEADER = 
      "https://code.usgs.gov/asl/sprockets/raw/master/";
  private static int LOGIN_PAGE_BYTE_SIZE = 7875;
  
  public static void 
  downloadTestData(String urlLoc, String localLoc, String fName) 
  throws IOException {
    String fullPath = URL_HEADER + urlLoc + fName;
    String localPath = "./test-data/sprockets/" + localLoc + fName;

    File target = new File(localPath);
    File fullDir = new File( target.getParent() );
    if ( !fullDir.exists() ) {
      fullDir.mkdirs();
    }
    if ( !target.exists() ) {
      URL website = new URL(fullPath);
      FileUtils.copyURLToFile(website, target);
      if (target.length() == LOGIN_PAGE_BYTE_SIZE) {
        // did we just download a login page?
        BufferedReader br = new BufferedReader( new FileReader(target) );
        String line = br.readLine();
        while ( line != null ) {
          line = line.toLowerCase();
          if (line.contains("sign in")) {
            br.close();
            target.delete();
            throw new IOException("Could not access file (repo not public?)");
          }
          line = br.readLine();
        }
        br.close();
      }
    }
    

  }
  
  // re-comment when repo made public
  //@Test
  public void canGetTestData() {
    String fullPath = URL_HEADER;
    try {
      URL web = new URL(fullPath);
      HttpURLConnection connection = (HttpURLConnection) web.openConnection();
      connection.connect();
      assertEquals( connection.getResponseCode(), 200 );
      String loc = "PSD_calculation/SyntheticData/";
      String file = "XX_KAS.00_BHZ.seed";
      downloadTestData(loc, "PSD/", file);
    } catch (MalformedURLException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    } catch (IOException e) {
      e.printStackTrace();
      fail();
    }
    
  }
  
  public static void makeTestDataDirectory() {
    File file = new File("./test-data/");
    if ( !file.exists() ) {
      file.mkdir();
    }
  }
}
