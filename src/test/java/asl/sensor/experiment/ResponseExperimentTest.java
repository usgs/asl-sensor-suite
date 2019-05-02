package asl.sensor.experiment;

import static asl.sensor.test.TestUtils.RESP_LOCATION;
import static org.junit.Assert.assertEquals;

import asl.sensor.input.DataStore;
import asl.sensor.input.InstrumentResponse;
import java.io.IOException;
import java.time.Instant;
import java.time.LocalDateTime;
import java.time.ZoneOffset;
import java.time.format.DateTimeFormatter;
import org.jfree.data.xy.XYSeriesCollection;
import org.junit.Test;

public class ResponseExperimentTest {

  @Test
  public void testMultipleEpochsPlotted() throws IOException {
    DataStore dataStore = new DataStore();
    ResponseExperiment respExperiment = new ResponseExperiment();
    DateTimeFormatter formatter =
        DateTimeFormatter.ofPattern("uuuu,DDD,HH:mm:ss").withZone(ZoneOffset.UTC);
    Instant startFirstEpoch =
        LocalDateTime.parse("2010,041,18:35:00", formatter).toInstant(ZoneOffset.UTC);
    Instant startSecondEpoch =
        LocalDateTime.parse("2015,055,00:00:00", formatter).toInstant(ZoneOffset.UTC);
    String respName = "RESP.CU.BCIP.00.BHZ_2017_268";
    String respFile = RESP_LOCATION + respName;
    System.out.println(InstrumentResponse.getRespFileEpochs(respFile));
    InstrumentResponse firstResp = new InstrumentResponse(respFile, startFirstEpoch);
    dataStore.setResponse(0, firstResp);
    InstrumentResponse secondResp = new InstrumentResponse(respFile, startSecondEpoch);
    dataStore.setResponse(1, secondResp);
    respExperiment.runExperimentOnData(dataStore);
    XYSeriesCollection firstPlottedData = respExperiment.getData().get(0);
    assertEquals(2, firstPlottedData.getSeriesCount());
    String[] chartRespNames = {respName + " [2010.041]", respName + " [2015.055]"};
    for (int i = 0; i < chartRespNames.length; ++i) {
      assertEquals(chartRespNames[i], firstPlottedData.getSeriesKey(i));
    }
  }

}
