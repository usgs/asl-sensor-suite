package asl.sensor.input;

import static org.junit.Assert.fail;

import edu.iris.dmc.seedcodec.CodecException;
import edu.sc.seis.seisFile.mseed.SeedFormatException;
import java.io.IOException;

public class DataStoreUtils {

  public static DataStore createFromNames(String respName, String calibrationInName,
      String sensorOutputName) {
    try {

      DataStore ds = new DataStore();

      if (calibrationInName != null) {
        ds.setBlock(0, calibrationInName);
      }
      if (sensorOutputName != null) {
        ds.setBlock(1, sensorOutputName);
      }

      if (respName != null) {
        InstrumentResponse ir = new InstrumentResponse(respName);
        ds.setResponse(1, ir);
      }

      return ds;
    } catch (IOException | SeedFormatException | CodecException e) {
      e.printStackTrace();
      fail();
    }
    //If it gets here, it has already failed.
    return new DataStore();
  }
}
