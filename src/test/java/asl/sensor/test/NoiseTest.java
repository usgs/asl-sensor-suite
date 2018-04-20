package asl.sensor.test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;

import asl.sensor.experiment.NoiseExperiment;
import asl.sensor.input.DataStore;
import asl.sensor.utils.TimeSeriesUtils;
import edu.iris.dmc.seedcodec.CodecException;
import edu.sc.seis.seisFile.mseed.SeedFormatException;
import java.io.IOException;
import java.time.OffsetDateTime;
import java.time.ZoneOffset;
import org.jfree.data.xy.XYDataItem;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.junit.Test;

public class NoiseTest {

  public static String folder = TestUtils.TEST_DATA_LOCATION + TestUtils.SUBPAGE;

  private XYSeriesCollection setUpTest1() {
    String testFolder = folder + "noise-neg159db/";
    String[] data = new String[3];
    data[0] = "00_BH0.512.seed";
    data[1] = "10_BH0.512.seed";
    data[2] = "TST6." + data[0];
    String resp = "T-compact_Q330HR_BH_40";
    DataStore ds = new DataStore();
    for (int i = 0; i < data.length; ++i) {
      try {
        ds.setBlock(i, testFolder + data[i]);
      } catch (IOException | SeedFormatException | CodecException e) {
        e.printStackTrace();
        fail();
      }
      ds.setEmbedResponse(i, resp);
    }

    OffsetDateTime startCal =
        OffsetDateTime.ofInstant(ds.getBlock(0).getStartInstant(), ZoneOffset.UTC);
    startCal = startCal.withHour(0);
    startCal = startCal.withMinute(59);
    startCal = startCal.withSecond(59);
    startCal = startCal.withNano(994 * TimeSeriesUtils.TO_MILLI_FACTOR);
    OffsetDateTime endCal = startCal.withHour(7);
    endCal = endCal.withMinute(0);
    endCal = endCal.withSecond(0);
    endCal = endCal.withNano(25 * TimeSeriesUtils.TO_MILLI_FACTOR);
    ds.trim(startCal.toInstant(), endCal.toInstant());
    NoiseExperiment ne = new NoiseExperiment();
    ne.setFreqSpace(false); // use period units (s)
    ne.runExperimentOnData(ds);
    return ne.getData().get(0);
  }

  private XYSeriesCollection setUpTest2() {
    String testFolder = folder + "noise-neg160db/";
    String[] data = new String[3];
    data[0] = "00_LH0.512.seed";
    data[1] = "10_LH0.512.seed";
    data[2] = "TST6." + data[0];
    String resp = "T-compact_Q330HR_BH_40";
    DataStore ds = new DataStore();
    for (int i = 0; i < data.length; ++i) {
      try {
        ds.setBlock(i, testFolder + data[i]);
      } catch (IOException | SeedFormatException | CodecException e) {
        e.printStackTrace();
        fail();
      }
      ds.setEmbedResponse(i, resp);
    }
    OffsetDateTime startCal =
        OffsetDateTime.ofInstant( ds.getBlock(0).getStartInstant(), ZoneOffset.UTC);
    startCal = startCal.withHour(0);
    startCal = startCal.withMinute(59);
    startCal = startCal.withSecond(59);
    startCal = startCal.withNano(994 * TimeSeriesUtils.TO_MILLI_FACTOR);
    OffsetDateTime endCal = startCal.withHour(7);
    endCal = endCal.withMinute(0);
    endCal = endCal.withSecond(0);
    endCal = endCal.withNano(25 * TimeSeriesUtils.TO_MILLI_FACTOR);
    ds.trim( startCal.toInstant(), endCal.toInstant() );

    NoiseExperiment ne = new NoiseExperiment();
    ne.setFreqSpace(false); // use period units (s)
    ne.runExperimentOnData(ds);
    return ne.getData().get(0);
  }

  @Test
  public void testResultsData1PSD1() {
    int idx = 0;
    double psdCheck = -159.73;
    double noiseCheck = -161.17;
    // everything below here same for every test
    XYSeriesCollection xysc = setUpTest1();
    // first 3 data, PSDs of each input
    // second 3 data, self-noise of each input
    // want data from 30 to 100s
    double low = 30.;
    double high = 100.;
    double psdResults = 0.;
    double noiseResults = 0.;
    XYSeries psd = xysc.getSeries(idx);
    XYSeries noise = xysc.getSeries(idx + 3);
    int psdPoints = 0;
    int noisePoints = 0;
    for (int j = 0; j < psd.getItemCount(); ++j) {
      XYDataItem psdxy = psd.getDataItem(j);
      double x = psdxy.getX().doubleValue();
      if (x >= low && x <= high) {
        psdResults += psdxy.getY().doubleValue();
        ++psdPoints;
      }
      XYDataItem noisxy = noise.getDataItem(j);
      x = noisxy.getX().doubleValue();
      if (x >= low && x <= high) {
        noiseResults += noisxy.getY().doubleValue();
        ++noisePoints;
      }
    }
    psdResults /= psdPoints;
    noiseResults /= noisePoints;

    System.out.println(psdResults + "," + noiseResults);
    System.out.println(psdCheck + "," + noiseCheck);
    System.out.println("PSD DIFF: " + Math.abs(psdResults - psdCheck));
    System.out.println("NOISE DIFF: " + Math.abs(noiseResults - noiseCheck));

    assertEquals(noiseCheck, noiseResults, 1E-2);
    assertEquals(psdCheck, psdResults, 1E-2);

  }

  @Test
  public void testResultsData1PSD2() {
    int idx = 1;
    double psdCheck = -161.19;
    double noiseCheck = -162.41;
    // everything below here same for every test
    XYSeriesCollection xysc = setUpTest1();
    // first 3 data, PSDs of each input
    // second 3 data, self-noise of each input
    // want data from 30 to 100s
    double low = 30.;
    double high = 100.;
    double psdResults = 0.;
    double noiseResults = 0.;
    XYSeries psd = xysc.getSeries(idx);
    XYSeries noise = xysc.getSeries(idx + 3);
    int psdPoints = 0;
    int noisePoints = 0;
    for (int j = 0; j < psd.getItemCount(); ++j) {
      XYDataItem psdxy = psd.getDataItem(j);
      double x = psdxy.getX().doubleValue();
      if (x >= low && x <= high) {
        psdResults += psdxy.getY().doubleValue();
        ++psdPoints;
      }
      XYDataItem noisxy = noise.getDataItem(j);
      x = noisxy.getX().doubleValue();
      if (x >= low && x <= high) {
        noiseResults += noisxy.getY().doubleValue();
        ++noisePoints;
      }
    }
    psdResults /= psdPoints;
    noiseResults /= noisePoints;
    assertEquals(noiseCheck, noiseResults, 1E-2);
    assertEquals(psdCheck, psdResults, 1E-2);

  }

  @Test
  public void testResultsData1PSD3() {
    int idx = 2;
    double psdCheck = -157.76;
    double noiseCheck = -158.66;
    // everything below here same for every test
    XYSeriesCollection xysc = setUpTest1();
    // first 3 data, PSDs of each input
    // second 3 data, self-noise of each input
    // want data from 30 to 100s
    double low = 30.;
    double high = 100.;
    double psdResults = 0.;
    double noiseResults = 0.;
    XYSeries psd = xysc.getSeries(idx);
    XYSeries noise = xysc.getSeries(idx + 3);
    int psdPoints = 0;
    int noisePoints = 0;
    for (int j = 0; j < psd.getItemCount(); ++j) {
      XYDataItem psdxy = psd.getDataItem(j);
      double x = psdxy.getX().doubleValue();
      if (x >= low && x <= high) {
        psdResults += psdxy.getY().doubleValue();
        ++psdPoints;
      }
      XYDataItem noisxy = noise.getDataItem(j);
      x = noisxy.getX().doubleValue();
      if (x >= low && x <= high) {
        noiseResults += noisxy.getY().doubleValue();
        ++noisePoints;
      }
    }
    psdResults /= psdPoints;
    noiseResults /= noisePoints;
    assertEquals(noiseCheck, noiseResults, 1E-2);
    assertEquals(psdCheck, psdResults, 1E-2);
  }

  @Test
  public void testResultsData2PSD1() {
    int idx = 0;
    double psdCheck = -159.77;
    double noiseCheck = -161.16;
    // everything below here same for every test
    XYSeriesCollection xysc = setUpTest2();
    // first 3 data, PSDs of each input
    // second 3 data, self-noise of each input
    // want data from 30 to 100s
    double low = 30.;
    double high = 100.;
    double psdResults = 0.;
    double noiseResults = 0.;
    XYSeries psd = xysc.getSeries(idx);
    XYSeries noise = xysc.getSeries(idx + 3);
    int psdPoints = 0;
    int noisePoints = 0;
    for (int j = 0; j < psd.getItemCount(); ++j) {
      XYDataItem psdxy = psd.getDataItem(j);
      double x = psdxy.getX().doubleValue();
      if (x >= low && x <= high) {
        psdResults += psdxy.getY().doubleValue();
        ++psdPoints;
      }
      XYDataItem noisxy = noise.getDataItem(j);
      x = noisxy.getX().doubleValue();
      if (x >= low && x <= high) {
        noiseResults += noisxy.getY().doubleValue();
        ++noisePoints;
      }
    }
    psdResults /= psdPoints;
    noiseResults /= noisePoints;
    assertEquals(noiseCheck, noiseResults, 1E-2);
    assertEquals(psdCheck, psdResults, 1E-2);

  }

}
