package asl.sensor.experiment;

import static asl.sensor.experiment.LagTimeExperiment.deconvolveResponse;
import static asl.sensor.experiment.LagTimeExperiment.weightedAverageSlopesInterp;
import static asl.utils.response.ResponseParser.loadEmbeddedResponse;
import static asl.utils.timeseries.TimeSeriesUtils.ONE_HZ_INTERVAL;
import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;

import asl.sensor.input.DataStore;
import asl.sensor.test.TestUtils;
import asl.utils.response.ChannelMetadata;
import asl.utils.response.ChannelMetadata.ResponseStageException;
import asl.utils.response.ResponseParser;
import asl.utils.response.ResponseParser.EpochIdentifier;
import asl.utils.response.ResponseUnits.ResolutionType;
import asl.utils.response.ResponseUnits.SensorType;
import asl.utils.timeseries.DataBlock;
import asl.utils.timeseries.TimeSeriesUtils;
import edu.iris.dmc.seedcodec.CodecException;
import edu.sc.seis.seisFile.mseed.SeedFormatException;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import org.jfree.data.xy.XYSeries;
import org.junit.Test;

public class LagTimeExperimentTest {

  private static final String folder = TestUtils.TEST_DATA_LOCATION + TestUtils.SUBPAGE;
  private static final long THOUSAND_SPS_INTERVAL = ONE_HZ_INTERVAL / 1000;

  @Test
  public void testProducesValidResult() throws SeedFormatException, CodecException, IOException {
    DataStore ds = new DataStore();
    ds.setBlock(0,
        "src/test/resources/seismic-test-data/seed_data/IU_ANMO/2018/005/00_BHZ.512.seed");
    ds.setBlock(1,
        "src/test/resources/seismic-test-data/seed_data/IU_ANMO/2018/005/00_BHZ.512.seed");
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
    XYSeries series = exp.getData().get(0).getSeries(0);
    int expectLength = 23901;
    assertEquals(expectLength, series.getItemCount());
    assertEquals(0, exp.getLagTime());
    int xValue = expectLength/2;
    assertEquals(series.getMaxY(), series.getY(xValue));
  }

  @Test
  public void testsCorrectLagIdentification_1Second()
      throws SeedFormatException, CodecException, IOException {
    String seedFile = folder + "lag-tests/Test_data_1s.seed"; // IU ANMO data
    List<String> traces = TimeSeriesUtils.getMplexNameList(seedFile);
    // should be in order 00_BHZ then 10_BHZ
    assertEquals("IU_ANMO_00_BHZ", traces.get(0));
    DataStore ds = new DataStore();
    for (int i = 0; i < traces.size(); ++i) {
      ds.setBlock(i, TimeSeriesUtils.getTimeSeries(seedFile, traces.get(i)));
    }
    ds.setResponse(0, folder + "lag-tests/RESP.IU.ANMO.00.BHZ");
    ds.setResponse(1, folder + "lag-tests/RESP.IU.ANMO.10.BHZ");
    LagTimeExperiment exp = new LagTimeExperiment();
    exp.runExperimentOnData(ds);
    // we'll do the comparison in seconds rather than ms
    assertEquals(-1., exp.getLagTime()/1000., 2E-3);
  }

  @Test
  public void testsCorrectLagIdentification_5Seconds()
      throws SeedFormatException, CodecException, IOException {
    String seedFile = folder + "lag-tests/Test_data_5s.seed"; // IU ANMO data
    List<String> traces = TimeSeriesUtils.getMplexNameList(seedFile);
    // should be in order 00_BHZ then 10_BHZ
    assertEquals("IU_ANMO_00_BHZ", traces.get(0));
    DataStore ds = new DataStore();
    for (int i = 0; i < traces.size(); ++i) {
      ds.setBlock(i, TimeSeriesUtils.getTimeSeries(seedFile, traces.get(i)));
    }
    ds.setResponse(0, folder + "lag-tests/RESP.IU.ANMO.00.BHZ");
    ds.setResponse(1, folder + "lag-tests/RESP.IU.ANMO.10.BHZ");
    LagTimeExperiment exp = new LagTimeExperiment();
    exp.runExperimentOnData(ds);
    // we'll do the comparison in seconds rather than ms
    assertEquals(-5., exp.getLagTime()/1000., 2E-3);
  }

  @Test
  public void testsCorrectLagIdentification_under1Second()
      throws SeedFormatException, CodecException, IOException {
    String seedFile = folder + "lag-tests/Test_data_XXs.seed"; // IU ANMO data
    List<String> traces = TimeSeriesUtils.getMplexNameList(seedFile);
    // should be in order 00_BHZ then 10_BHZ
    assertEquals("IU_ANMO_00_BHZ", traces.get(0));
    DataStore ds = new DataStore();
    for (int i = 0; i < traces.size(); ++i) {
      ds.setBlock(i, TimeSeriesUtils.getTimeSeries(seedFile, traces.get(i)));
    }
    ds.setResponse(0, folder + "lag-tests/RESP.IU.ANMO.00.BHZ");
    ds.setResponse(1, folder + "lag-tests/RESP.IU.ANMO.10.BHZ");
    LagTimeExperiment exp = new LagTimeExperiment();
    exp.runExperimentOnData(ds);
    // we'll do the comparison in seconds rather than ms
    assertEquals(-0.100, exp.getLagTime()/1000., 1E-2);
  }

  @Test
  public final void testWeightedAverageSlopesInterpQuadruple()
      throws SeedFormatException, CodecException, IOException {
    String filename = "asl-java-utils/src/test/resources/fft-full/00_LHZ.512.seed";
    // some ANMO LHZ data
    DataBlock db = TimeSeriesUtils.getFirstTimeSeries(filename);
    double[] data = Arrays.copyOfRange(db.getData(), 0, 300);
    long resampleRate = db.getInitialInterval() / 4;
    double[] interpolated =
        weightedAverageSlopesInterp(data, db.getInitialInterval(), resampleRate);
    assertEquals(1197, interpolated.length);
  }

  @Test
  public void testDifferentiation() {
    double[] data = new double[100];
    for (int i = 0; i < data.length; ++i) {
      data[i] = Math.pow(i, 3);
    }
    double[] diff = LagTimeExperiment.differentiate(data, 1);
    double[] expected = {1, 7, 19, 37, 61, 91, 127, 169, 217, 271, 331, 397, 469, 547, 631, 721,
        817, 919, 1027, 1141, 1261, 1387, 1519, 1657, 1801, 1951, 2107, 2269, 2437, 2611, 2791,
        2977, 3169, 3367, 3571, 3781, 3997, 4219, 4447, 4681, 4921, 5167, 5419, 5677, 5941, 6211,
        6487, 6769, 7057, 7351, 7651, 7957, 8269, 8587, 8911, 9241, 9577, 9919, 10267, 10621, 10981,
        11347, 11719, 12097, 12481, 12871, 13267, 13669, 14077, 14491, 14911, 15337, 15769, 16207,
        16651, 17101, 17557, 18019, 18487, 18961, 19441, 19927, 20419, 20917, 21421, 21931, 22447,
        22969, 23497, 24031, 24571, 25117, 25669, 26227, 26791, 27361, 27937, 28519, 29107};
    assertArrayEquals(expected, diff, 0.);
  }

  @Test
  public void testCorrelation() {
    double[] first = {1, 2, 5, 6, 8, 9, 12, 13};
    double[] second = {8, 7, 6, 9, 8, 7, 2, 1, 2, 1, 2};
    double[] expected = {104, 187, 234, 316, 370, 398, 341, 270, 224, 168, 136, 87, 54,
        36, 22, 14, 5, 2};
    double[] result = LagTimeExperiment.getCorrelation(first, second);

    assertEquals(18, result.length);
    assertArrayEquals(expected, result, 0.);
  }

  @Test
  public void testResponseDeconvolution()
      throws SeedFormatException, IOException, CodecException, ResponseStageException {
    String seedFile = folder + "lag-tests/Test_data_1s.seed"; // IU ANMO data
    List<String> traces = TimeSeriesUtils.getMplexNameList(seedFile);
    // should be in order 00_BHZ then 10_BHZ
    assertEquals("IU_ANMO_00_BHZ", traces.get(0));
    DataBlock db = TimeSeriesUtils.getTimeSeries(seedFile, traces.get(0));
    String respFile = folder + "lag-tests/RESP.IU.ANMO.00.BHZ";
    EpochIdentifier epoch =  ResponseParser.getRespFileClosestEpoch(respFile, db.getStartTime(),
        db.getEndTime());
    ChannelMetadata resp = ResponseParser.parseResponse(respFile, epoch.filePointer);
    double[] data = db.getData();
    double sampleRate = db.getSampleRate();
    double[] deconvolved = deconvolveResponse(data, resp, sampleRate);
    // these values match obspy within 4 significant figures
    double[] referenceData = {5.670340176262498E-8, 6.12698644469828E-8, 5.432928106370917E-8,
        4.764814385177327E-8, 4.214522525619337E-8, 3.652414727623224E-8, 3.068923397265336E-8,
        2.7142722266613733E-8, 2.786476614032104E-8, 2.6870732011867994E-8, 1.843696868286195E-8,
        8.046415930251987E-9, 3.3936456994508232E-9, 1.9865677144096453E-9, -2.0947023160831E-9,
        -6.772540148534291E-9, -1.1926652308888597E-8, -1.909082587492656E-8,
        -2.3325641021675717E-8, -2.660409245167799E-8, -2.8987782490221795E-8,
        -3.311331703941894E-8, -3.935194254986788E-8, -4.48440559485019E-8, -5.491579336671025E-8};
    for (int i = 0; i < referenceData.length; ++i) {
      assertEquals(referenceData[i], deconvolved[i], 0.);
    }
  }
}
