package asl.sensor.input;

import static org.junit.Assert.fail;

import edu.iris.dmc.seedcodec.CodecException;
import edu.sc.seis.seisFile.mseed.SeedFormatException;
import java.io.IOException;

public class DataStoreUtils {

  public static DataStore createFromNames(String respName, String calibrationInName,
      String sensorOutputName) {
    try {
      InstrumentResponse ir = new InstrumentResponse(respName);

      DataStore ds = new DataStore();

      ds.setBlock(0, calibrationInName);
      ds.setBlock(1, sensorOutputName);

      ds.setResponse(1, ir);

      return ds;
    } catch (IOException | SeedFormatException | CodecException e) {
      e.printStackTrace();
      fail();
    }
    //Never gets here
    return new DataStore();
  }
}
