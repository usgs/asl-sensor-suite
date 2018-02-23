package asl.sensor.test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.time.OffsetDateTime;
import java.time.ZoneOffset;
import org.jfree.data.xy.XYDataItem;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.junit.Test;
import asl.sensor.experiment.NoiseExperiment;
import asl.sensor.input.DataStore;
import asl.sensor.input.InstrumentResponse;
import asl.sensor.utils.TimeSeriesUtils;
import edu.iris.dmc.seedcodec.CodecException;
import edu.sc.seis.seisFile.mseed.SeedFormatException;

public class NoiseTest {

  XYSeriesCollection xysc;

  //@Test
  public void testResp() {
    String resp = "T-compact_Q330HR_BH_40";
    try {
      InstrumentResponse ir = InstrumentResponse.loadEmbeddedResponse(resp);
      System.out.println(ir.toString());
    } catch (IOException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
      fail();
    }

  }

  public XYSeriesCollection setUpTest1() throws FileNotFoundException {
    String folder = "test-data/noise-neg159db/";
    String[] data = new String[3];
    data[0] = "00_BH0.512.seed";
    data[1] = "10_BH0.512.seed";
    data[2] = "TST6." + data[0];
    String resp = "T-compact_Q330HR_BH_40";
    DataStore ds = new DataStore();
    for (int i = 0; i < data.length; ++i) {
      try {
        ds.setBlock(i, folder + data[i]);
      } catch (SeedFormatException | CodecException e) {
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
    XYSeriesCollection xysc = ne.getData().get(0);
    return xysc;
  }

  public XYSeriesCollection setUpTest2() throws FileNotFoundException {
    String folder = "test-data/noise-neg160db/";
    String[] data = new String[3];
    data[0] = "00_LH0.512.seed";
    data[1] = "10_LH0.512.seed";
    data[2] = "TST6." + data[0];
    String resp = "T-compact_Q330HR_BH_40";
    DataStore ds = new DataStore();
    for (int i = 0; i < data.length; ++i) {
      try {
        ds.setBlock(i, folder + data[i]);
      } catch (SeedFormatException | CodecException e) {
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
    XYSeriesCollection xysc = ne.getData().get(0);
    return xysc;
  }

  @Test
  public void testResultsData1PSD1() {
    int idx = 0;
    double psdCheck = -157.64;
    double noiseCheck = -158.53;
    // everything below here same for every test
    try{
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

    } catch (FileNotFoundException e) {
      e.printStackTrace();
      fail();
    }
  }

  @Test
  public void testResultsData2PSD1() {
    int idx = 0;
    double psdCheck = -158.68;
    double noiseCheck = -159.63;
    // everything below here same for every test
    try{
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

    } catch (FileNotFoundException e) {
      e.printStackTrace();
      fail();
    }
  }

  @Test
  public void testResultsData1PSD2() {
    int idx = 1;
    double psdCheck = -159.05;
    double noiseCheck = -160.79;
    // everything below here same for every test
    try{
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

    } catch (FileNotFoundException e) {
      e.printStackTrace();
      fail();
    }
  }

  @Test
  public void testResultsData1PSD3() {
    int idx = 2;
    double psdCheck = -155.56;
    double noiseCheck = -156.33;
    // everything below here same for every test
    try{
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

    } catch (FileNotFoundException e) {
      e.printStackTrace();
      fail();
    }
  }

}
