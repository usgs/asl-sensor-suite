package asl.sensor.test;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.net.MalformedURLException;
import java.net.URL;
import java.nio.file.FileSystems;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;

public class TestUtils {
  
  private static String URL_HEADER = 
      "https://code.usgs.gov/asl/sprockets/raw/master/";
  
  public static void downloadTestData(String locationName) {
    try {
      String fullPath = URL_HEADER + locationName;
      String localName = "./test-data/" + locationName;
      URL web = new URL(fullPath);
      File target = new File(localName);
      if ( !target.exists() ) {
        try (InputStream in = web.openStream() ) {
          Path path = FileSystems.getDefault().getPath(localName);
          Files.copy(in, path, StandardCopyOption.REPLACE_EXISTING);
        } catch (IOException e) {
          // TODO Auto-generated catch block
          e.printStackTrace();
        }
      }
    } catch (MalformedURLException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
  }
  
  public static void makeTestDataDirectory() {
    File file = new File("./test-data/");
    if ( !file.exists() ) {
      file.mkdir();
    }
  }
}
