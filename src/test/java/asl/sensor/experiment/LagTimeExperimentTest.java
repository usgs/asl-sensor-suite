package asl.sensor.experiment;

import static asl.utils.response.ResponseParser.loadEmbeddedResponse;
import static org.junit.Assert.assertEquals;

import asl.sensor.input.DataStore;
import asl.utils.response.ResponseUnits.ResolutionType;
import asl.utils.response.ResponseUnits.SensorType;
import asl.utils.timeseries.DataBlock;
import edu.iris.dmc.seedcodec.CodecException;
import edu.sc.seis.seisFile.mseed.SeedFormatException;
import java.io.IOException;
import org.jfree.data.xy.XYSeries;
import org.junit.Test;

public class LagTimeExperimentTest {

  @Test
  public void testProducesValidResult() throws SeedFormatException, CodecException, IOException {
    DataStore ds = new DataStore();
    ds.setBlock(0,
        "src/test/resources/seismic-test-data/seed_data/IU_ANMO/2018/005/00_BHZ.512.seed");
    ds.setBlock(1,
        "src/test/resources/seismic-test-data/seed_data/IU_ANMO/2018/005/10_BHZ.512.seed");
    ds.setResponse(0, loadEmbeddedResponse(SensorType.STS2gen3, ResolutionType.HIGH));
    ds.setResponse(1, loadEmbeddedResponse(SensorType.STS2gen3, ResolutionType.HIGH));
    {
      DataBlock first = ds.getBlock(0);
      long start = first.getStartTime();
      long end = start + first.getInitialInterval() * 240;
      ds.trim(start, end);
    }
    LagTimeExperiment exp = new LagTimeExperiment();
    exp.runExperimentOnData(ds);
    assertEquals(0, exp.getLagTime());
    XYSeries series = exp.getData().get(0).getSeries(0);
    assertEquals(10000, exp.getData().get(0).getSeries(0).getItemCount());
  }
}
