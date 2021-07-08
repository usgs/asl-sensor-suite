package asl.sensor.experiment;

import static asl.sensor.experiment.LagTimeExperiment.deconvolveResponse;
import static asl.sensor.experiment.LagTimeExperiment.weightedAverageSlopesInterp;
import static asl.utils.NumericUtils.demeanInPlace;
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
  public final void testWeightedAverageSlopesInterp()
      throws SeedFormatException, CodecException, IOException {
    String filename =  "asl-java-utils/src/test/resources/fft-full/00_LHZ.512.seed";
      // some ANMO LHZ data
    DataBlock db = TimeSeriesUtils.getFirstTimeSeries(filename);
    double[] data = Arrays.copyOfRange(db.getData(), 0, 300);
    long resampleRate = db.getInitialInterval() / 2;
    double[] interpolated =
        weightedAverageSlopesInterp(data, db.getInitialInterval(), resampleRate);
    assertEquals(599, interpolated.length);
    // these values taken from running obspy data
    double[] expected = {-514397., -514614.84973226, -514819., -515027.40026774, -515144.,
        -515124.49315068, -515080., -514980.42351598, -514852., -514733.44270833, -514648.,
        -514611.640625, -514596., -514732.34602273, -514913., -514996.65397727, -515036.,
        -514912.71280603, -514725., -514582.78719397, -514505., -514564.93012422, -514691.,
        -514929.77161491, -515149., -515230.79826087, -515266., -515062.12021858, -514746.,
        -514492.37978142, -514351., -514450.57450331, -514625., -514803.22314123, -514955.,
        -515047.70235546, -515092., -515046.49379771, -514949., -514732.50620229, -514568.,
        -514585.33766234, -514628., -514783.60918346, -514953., -515024.05315421, -515056.,
        -515002.15412186, -514890., -514706.76412631, -514498., -514309.58175182, -514205.,
        -514307.51077586, -514520., -514903.56108537, -515249., -515374.42813877, -515428.,
        -515290.77884615, -515073., -514909.63725362, -514778., -514683.14122507, -514632.,
        -514625.87510759, -514621., -514609.04673423, -514595., -514581.02083333, -514573.,
        -514592.20568928, -514640., -514795.90683972, -515030., -515319.387471, -515502.,
        -515303.2079489, -514935., -514429.2920511, -514093., -514187.67776584, -514382.,
        -514724.08821881, -515024., -515122.73401535, -515164., -515000.53157895, -514758.,
        -514588.96842105, -514499., -514756.06276446, -515069., -515166.43723554, -515208.,
        -514824.41114983, -514336., -514145.58885017, -514060., -514360.47845497, -514817.,
        -515160.52154503, -515348., -515284.85039894, -515175., -515066.32829234, -514972.,
        -514926.58164936, -514877., -514739.27936763, -514561., -514371.96029173, -514260.,
        -514263.59839357, -514274., -514445.25280643, -514758., -515170.48463584, -515524.,
        -515685.16416416, -515757., -515524.95618051, -515106., -514614.10138681, -514228.};
    for (int i = 0; i < expected.length; ++i) {
      assertEquals(expected[i], interpolated[i], 1E-6);
    }
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
  public final void testWeightedAverageSlopesInterp_bhzData()
      throws SeedFormatException, CodecException, IOException {
    String seedFile = folder + "lag-tests/Test_data_XXs.seed"; // IU ANMO data
    DataBlock db = TimeSeriesUtils.getFirstTimeSeries(seedFile);
    double[] data = db.getData();
    demeanInPlace(data);
    double[] interpolated =
        weightedAverageSlopesInterp(data, db.getInitialInterval(), THOUSAND_SPS_INTERVAL);
    assertEquals(99951, interpolated.length);

    double[] expectedInterpolated = {246.897, 246.42801863, 245.94162578, 245.43915729,
        244.92160302, 244.38993841, 243.8455209, 243.28932914, 242.7224669, 242.14589944,
        241.56100562, 240.96875285, 240.37010154, 239.76644107, 239.1587332, 238.54793582,
        237.93544412, 237.32221656, 236.70935691, 236.0978234, 235.48901158, 234.8838792,
        234.28338754, 233.68892687, 233.10145719, 232.52194516, 231.95177106, 231.39189949,
        230.84330489, 230.30735278, 229.78501474, 229.27739454, 228.78548064, 228.31061156,
        227.85377185, 227.41596328, 226.99849857, 226.60237477, 226.22860932, 225.8784844,
        225.55301193, 225.25329569, 224.98037767, 224.73549196, 224.51967401, 224.33398701,
        224.17962347, 224.05763924, 223.96912106, 223.9152151, 223.897, 223.897, 223.897,
        223.897, 223.897, 223.897, 223.897, 223.897, 223.897, 223.897, 223.897, 223.897, 223.897,
        223.897, 223.897, 223.897, 223.897, 223.897, 223.897, 223.897, 223.897, 223.897, 223.897,
        223.897, 223.897, 223.897, 223.897, 223.897, 223.897, 223.897, 223.897, 223.897, 223.897,
        223.897, 223.897, 223.897, 223.897, 223.897, 223.897, 223.897, 223.897, 223.897, 223.897,
        223.897, 223.897, 223.897, 223.897, 223.897, 223.897, 223.897};

    for (int i = 0; i < expectedInterpolated.length; ++i) {
      assertEquals("Array mismatch at index " + i,
          expectedInterpolated[i], interpolated[i], 1E-4);
    }
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
