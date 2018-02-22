package asl.sensor.test;

import static org.junit.Assert.fail;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.net.URL;
import java.time.Instant;
import java.time.OffsetDateTime;
import java.time.ZoneOffset;
import org.apache.commons.io.FileUtils;
import org.junit.Test;
import asl.sensor.input.DataStore;

public class TestUtils {

  // may need to change this to deal with eventual migration to main usgs github
  private static String URL_HEADER =
      "https://github.com/kschramm-usgs/sprockets/raw/master/";
  private static int LOGIN_PAGE_BYTE_SIZE = 7875;

  // inName and fName are separated to make it possible to rename output file
  // to prevent collisions between data with the same filename from different inputs
  public static void
  downloadTestData(String urlLoc, String inName, String localLoc, String fName)
  throws IOException {
    // System.out.println("Acquiring data from " + urlLoc);
    String fullPath = URL_HEADER + urlLoc + inName;
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

  @Test
  public void canGetTestData() {

    String loc = "PSD_calculation/SyntheticData/";
    String file = "XX_KAS.00_BHZ.seed";
    try {
      downloadTestData(loc, file, "PSD/", file);
    } catch (IOException e1) {
      e1.printStackTrace();
      fail();
    }

  }

  public static OffsetDateTime getStartCalendar(DataStore ds) {
    Instant time = Instant.ofEpochMilli( ds.getBlock(0).getStartTime() );
    OffsetDateTime dt = OffsetDateTime.ofInstant(time, ZoneOffset.UTC);
    return dt;
  }

  public static OffsetDateTime getEndCalendar(DataStore ds) {
    Instant time = Instant.ofEpochMilli( ds.getBlock(0).getEndTime() );
    OffsetDateTime dt = OffsetDateTime.ofInstant(time, ZoneOffset.UTC);
    return dt;
  }

  public static void makeTestDataDirectory() {
    File file = new File("./test-data/");
    if ( !file.exists() ) {
      file.mkdir();
    }
  }
}
