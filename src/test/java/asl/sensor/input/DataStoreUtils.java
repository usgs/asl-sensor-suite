package asl.sensor.input;

import static org.junit.Assert.fail;

import asl.utils.input.InstrumentResponse;
import edu.iris.dmc.seedcodec.CodecException;
import edu.sc.seis.seisFile.mseed.SeedFormatException;
import java.io.IOException;

public class DataStoreUtils {

  public static DataStore createFromNames(String respName, String calibrationInName,
      String sensorOutputName) {
    try {
      DataStore dataStore = new DataStore();

      if (calibrationInName != null) {
        dataStore.setBlock(0, calibrationInName);
      }
      if (sensorOutputName != null) {
        dataStore.setBlock(1, sensorOutputName);
      }

      if (respName != null) {
        InstrumentResponse ir = new InstrumentResponse(respName);
        dataStore.setResponse(1, ir);
      }

      return dataStore;
    } catch (IOException | SeedFormatException | CodecException e) {
      e.printStackTrace();
      fail(e.getMessage());
    }
    //If it gets here, it has already failed.
    return new DataStore();
  }

  public static DataStore createFromNamesEmbedResp(String respName, String calibrationInName,
      String sensorOutputName) {
    try {
      DataStore dataStore = new DataStore();

      if (calibrationInName != null) {
        dataStore.setBlock(0, calibrationInName);
      }
      if (sensorOutputName != null) {
        dataStore.setBlock(1, sensorOutputName);
      }

      if (respName != null) {
        InstrumentResponse ir = InstrumentResponse.loadEmbeddedResponse(respName);
        dataStore.setResponse(1, ir);
      }

      return dataStore;
    } catch (IOException | SeedFormatException | CodecException e) {
      e.printStackTrace();
      fail(e.getMessage());
    }
    //If it gets here, it has already failed.
    return new DataStore();
  }

  public static DataStore appendFromNames(DataStore dataStore, String calibrationInName,
      String sensorOutputName) {
    try {
      if (calibrationInName != null) {
        dataStore.getBlock(0).appendTimeSeries(calibrationInName);
      }
      if (sensorOutputName != null) {
        dataStore.getBlock(1).appendTimeSeries(sensorOutputName);
      }
      return dataStore;
    } catch (IOException | SeedFormatException | CodecException e) {
      e.printStackTrace();
      fail(e.getMessage());
    }
    //If it gets here, it has already failed.
    return new DataStore();
  }
}
