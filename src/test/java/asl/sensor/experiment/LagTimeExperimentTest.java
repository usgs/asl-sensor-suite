package asl.sensor.experiment;

import static asl.sensor.experiment.LagTimeExperiment.weightedAverageSlopesInterp;
import static asl.utils.response.ResponseParser.loadEmbeddedResponse;
import static asl.utils.timeseries.TimeSeriesUtils.ONE_HZ_INTERVAL;
import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;

import asl.sensor.input.DataStore;
import asl.utils.response.ResponseUnits.ResolutionType;
import asl.utils.response.ResponseUnits.SensorType;
import asl.utils.timeseries.DataBlock;
import asl.utils.timeseries.TimeSeriesUtils;
import edu.iris.dmc.seedcodec.CodecException;
import edu.sc.seis.seisFile.mseed.SeedFormatException;
import java.io.IOException;
import java.util.Arrays;
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
    XYSeries series = exp.getData().get(0).getSeries(0);
    int expectLength = 23901;
    assertEquals(expectLength, series.getItemCount());
    assertEquals(0, exp.getLagTime());
    int xValue = expectLength/2;
    assertEquals(series.getMaxY(), series.getY(xValue));
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
  public void testDifferentiation() {
    double[] data = new double[100];
    for (int i = 0; i < data.length; ++i) {
      data[i] = Math.pow(i, 3);
    }
    double[] diff = LagTimeExperiment.differentiate(data, ONE_HZ_INTERVAL);
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
    double[] expected = {0, 0, 0, 104, 187, 234, 316, 370, 398, 341, 270, 224, 168, 136, 87, 54,
        36, 22, 14, 5, 2};
    assertArrayEquals(expected, LagTimeExperiment.getCorrelation(first, second), 0.);
  }
}
