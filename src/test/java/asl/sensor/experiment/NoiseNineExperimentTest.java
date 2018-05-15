package asl.sensor.experiment;

import static org.junit.Assert.assertTrue;

import asl.sensor.ExperimentFactory;
import asl.sensor.gui.ExperimentPanel;
import asl.sensor.gui.NoiseNinePanel;
import asl.sensor.input.DataStore;
import asl.sensor.test.TestUtils;
import asl.sensor.utils.ReportingUtils;
import java.awt.Font;
import java.io.File;
import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.util.List;
import org.apache.pdfbox.pdmodel.PDDocument;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.LogarithmicAxis;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.data.xy.XYSeriesCollection;
import org.junit.Test;

public class NoiseNineExperimentTest {

  public static String folder = TestUtils.TEST_DATA_LOCATION + TestUtils.SUBPAGE;

  private String currentDir = System.getProperty("user.dir");

  @Test
  public void canRunAndPlotTest1() throws Exception {

    String testFolder = folder + "noisenine/";
    String[] types = new String[]{"00", "10", "60"};
    String freqName = "_BH";
    String[] components = new String[]{"1", "2", "Z"};
    String ending = ".512.seed";
    String respName = TestUtils.RESP_LOCATION + "STS-1_Q330HR_BH_20";

    DataStore ds =
        setUpTest(testFolder, types, freqName, components, ending, respName, false);

    SimpleDateFormat sdf = ExperimentPanel.DATE_TIME_FORMAT.get();

    Calendar cCal = Calendar.getInstance(sdf.getTimeZone());
    cCal.setTimeInMillis(ds.getBlock(0).getStartTime());
    cCal.set(Calendar.HOUR, 7);
    long start = cCal.getTime().getTime();
    cCal.set(Calendar.HOUR, 8);
    cCal.set(Calendar.MINUTE, 0);
    long end = cCal.getTime().getTime();

    ds.trim(start, end);

    NoiseNineExperiment nne = new NoiseNineExperiment();
    assertTrue(nne.hasEnoughData(ds));
    nne.runExperimentOnData(ds);

    List<XYSeriesCollection> xysc = nne.getData();
    JFreeChart[] jfcl = new JFreeChart[xysc.size()];

    String xAxisTitle = "Frequency (Hz)";
    NumberAxis xAxis = new LogarithmicAxis(xAxisTitle);
    Font bold = xAxis.getLabelFont().deriveFont(Font.BOLD);
    xAxis.setLabelFont(bold);
    String yAxisTitle = "Power (rel. 1 (m/s^2)^2/Hz)";
    String[] orientations = new String[]{" (North)", " (East)", " (Vertical)"};

    for (int i = 0; i < xysc.size(); ++i) {
      jfcl[i] = ChartFactory.createXYLineChart(
          ExperimentFactory.RANDOMCAL.getName() + orientations[i],
          xAxisTitle,
          yAxisTitle,
          xysc.get(i),
          PlotOrientation.VERTICAL,
          true,
          false,
          false);

      XYPlot xyp = jfcl[i].getXYPlot();

      xyp.setDomainAxis(xAxis);
    }

    StringBuilder sb = new StringBuilder();
    String insets = NoiseNinePanel.getInsetString(nne);
    sb.append('\n');
    sb.append(NoiseNinePanel.getTimeStampString(nne));
    sb.append('\n');
    sb.append("INPUTTED FILES:");
    sb.append('\n');

    List<String> names = nne.getInputNames();

    for (String name : names) {
      sb.append(name);
      sb.append('\n');
    }

    int width = 1280;
    int height = 960;

    PDDocument pdf = new PDDocument();
    ReportingUtils.chartsToPDFPage(width, height, pdf, jfcl);
    ReportingUtils.textListToPDFPages(pdf, insets, sb.toString());

    String testResultFolder = currentDir + "/testResultImages/";
    File dir = new File(testResultFolder);
    if (!dir.exists()) {
      dir.mkdir();
    }

    String testResult = testResultFolder + "Nine-Noise-Test.pdf";
    pdf.save(new File(testResult));
    pdf.close();
  }

  @Test
  public void canRunAndPlotTest2() throws Exception {
    String testFolder = folder + "noisenine2/";
    String[] types = new String[]{"00", "10", "30"};
    String freqName = "_BH";
    String[] components = new String[]{"1", "2", "Z"};
    String ending = ".512.seed";
    String respName = folder + "noisenine2/RESP.XX.MOFO.00.BHZ";

    DataStore ds =
        setUpTest(testFolder, types, freqName, components, ending, respName, false);

    SimpleDateFormat sdf = ExperimentPanel.DATE_TIME_FORMAT.get();

    Calendar cCal = Calendar.getInstance(sdf.getTimeZone());
    cCal.setTimeInMillis(ds.getBlock(0).getStartTime());
    cCal.set(Calendar.HOUR_OF_DAY, 12);
    cCal.set(Calendar.MINUTE, 0);
    cCal.set(Calendar.SECOND, 0);
    System.out.println("start: " + sdf.format(cCal.getTime()));
    long start = cCal.getTime().getTime();
    cCal.set(Calendar.HOUR_OF_DAY, 13);
    cCal.set(Calendar.MINUTE, 0);
    cCal.set(Calendar.SECOND, 0);
    System.out.println("end: " + sdf.format(cCal.getTime()));
    long end = cCal.getTime().getTime();

    ds.trim(start, end, DataStore.FILE_COUNT);

    NoiseNineExperiment nne = new NoiseNineExperiment();
    assertTrue(nne.hasEnoughData(ds));
    nne.runExperimentOnData(ds);

    List<XYSeriesCollection> xysc = nne.getData();
    JFreeChart[] jfcl = new JFreeChart[xysc.size()];

    String xAxisTitle = "Frequency (Hz)";
    NumberAxis xAxis = new LogarithmicAxis(xAxisTitle);
    Font bold = xAxis.getLabelFont().deriveFont(Font.BOLD);
    xAxis.setLabelFont(bold);
    String yAxisTitle = "Power (rel. 1 (m/s^2)^2/Hz)";
    String[] orientations = new String[]{" (North)", " (East)", " (Vertical)"};

    for (int i = 0; i < xysc.size(); ++i) {

      // make sure each series is populated with data from correct orientation
      // i.e., BH1s for north, BH2s for south, etc.
      String check = freqName + components[i];
      XYSeriesCollection coll = xysc.get(i);
      // has 3 components, so do 3 checks
      // nlnm, nhnm should be 4th and 5th entries in data
      for (int j = 0; j < 3; ++j) {
        String key = (String) coll.getSeriesKey(j);
        assertTrue(key.contains(check));
      }

      jfcl[i] = ChartFactory.createXYLineChart(
          ExperimentFactory.RANDOMCAL.getName() + orientations[i],
          xAxisTitle,
          yAxisTitle,
          coll,
          PlotOrientation.VERTICAL,
          true,
          false,
          false);

      XYPlot xyp = jfcl[i].getXYPlot();
      xyp.setDomainAxis(xAxis);
    }

    String insets = NoiseNinePanel.getInsetString(nne);
    StringBuilder sb = new StringBuilder();
    sb.append(NoiseNinePanel.getTimeStampString(nne));
    sb.append('\n');
    sb.append("INPUTTED FILES:");
    sb.append('\n');

    List<String> names = nne.getInputNames();

    for (String name : names) {
      sb.append(name);
      sb.append('\n');
    }

    int width = 1280;
    int height = 960;

    PDDocument pdf = new PDDocument();
    ReportingUtils.chartsToPDFPage(width, height, pdf, jfcl);
    ReportingUtils.textListToPDFPages(pdf, insets, sb.toString());

    String testResultFolder = currentDir + "/testResultImages/";
    File dir = new File(testResultFolder);
    if (!dir.exists()) {
      dir.mkdir();
    }

    String testResult = testResultFolder + "Nine-Noise-Test-2.pdf";
    pdf.save(new File(testResult));
    pdf.close();
  }

  private DataStore setUpTest(String folder, String[] types, String freqName,
      String[] components, String ending, String respName, boolean isEmbed) throws Exception {

    DataStore ds = new DataStore();

    for (int i = 0; i < types.length; ++i) {
      for (int j = 0; j < components.length; ++j) {
        int indexInStore = i * types.length + j;

        String fName = folder + types[i] + freqName + components[j] + ending;
        ds.setBlock(indexInStore, fName);
        if (isEmbed) {
          ds.setEmbedResponse(indexInStore, respName);
        } else {
          ds.setResponse(indexInStore, respName);
        }
      }
    }
    return ds;
  }
}
