package asl.sensor.gui;

import static java.nio.file.StandardCopyOption.REPLACE_EXISTING;

import asl.sensor.ExperimentFactory;
import asl.sensor.experiment.ResponseExperiment;
import asl.sensor.input.DataStore;
import asl.sensor.input.InstrumentResponse;
import asl.sensor.utils.ReportingUtils;
import java.awt.Color;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Set;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JDialog;
import javax.swing.JOptionPane;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.LogarithmicAxis;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.data.xy.XYSeriesCollection;

/**
 * Panel used to display response curves. Includes plots of response magnitude
 * and argument (rotation in complex space). Because this includes two charts
 * and also does not require any input timeseries, there are more overridden
 * methods in this class than in most other panels, including the image output.
 * The chart to be displayed at any moment is selected with a combobox.
 *
 * @author akearns
 */
public class ResponsePanel extends ExperimentPanel {

  private static final long serialVersionUID = 1L;

  private final ValueAxis frequencyAxis, degreeAxis;

  private final JCheckBox freqSpaceBox;
  private final JComboBox<String> plotSelection;

  private final JButton copyEmbeddedResp;

  private JFreeChart magnitudeChart, argumentChart;

  public ResponsePanel(ExperimentFactory experiment) {
    super(experiment);

    for (int i = 0; i < 3; ++i) {
      channelType[i] = "Response data (SEED data not used)";
    }

    String xAxisTitle = "Period (s)";
    String freqAxisTitle = "Frequency (Hz)";
    String yAxisTitle = "10 * log10( RESP(f) )";
    String degreeAxisTitle = "phi(RESP(f))";

    xAxis = new LogarithmicAxis(xAxisTitle);
    frequencyAxis = new LogarithmicAxis(freqAxisTitle);
    xAxis.setAutoRange(true);
    frequencyAxis.setAutoRange(true);

    yAxis = new NumberAxis(yAxisTitle);
    yAxis.setAutoRange(true);

    degreeAxis = new NumberAxis(degreeAxisTitle);
    degreeAxis.setAutoRange(true);
    ((NumberAxis) degreeAxis).setAutoRangeIncludesZero(false);

    ((NumberAxis) yAxis).setAutoRangeIncludesZero(false);
    Font bold = xAxis.getLabelFont();
    bold = bold.deriveFont(Font.BOLD, bold.getSize() + 2);
    xAxis.setLabelFont(bold);
    yAxis.setLabelFont(bold);
    frequencyAxis.setLabelFont(bold);
    degreeAxis.setLabelFont(bold);

    freqSpaceBox = new JCheckBox("Use Hz units (requires regen)");
    freqSpaceBox.setSelected(true);

    copyEmbeddedResp = new JButton("Extract an embedded response for editing");
    copyEmbeddedResp.addActionListener(this);

    plotSelection = new JComboBox<>();
    plotSelection.addItem(ResponseExperiment.MAGNITUDE);
    plotSelection.addItem(ResponseExperiment.ARGUMENT);
    plotSelection.addActionListener(this);

    applyAxesToChart(); // now that we've got axes defined

    // set the GUI components
    this.setLayout(new GridBagLayout());
    GridBagConstraints constraints = new GridBagConstraints();

    constraints.fill = GridBagConstraints.BOTH;
    constraints.gridx = 0;
    constraints.gridy = 0;
    constraints.weightx = 1.0;
    constraints.weighty = 1.0;
    constraints.gridwidth = 3;
    constraints.anchor = GridBagConstraints.CENTER;
    this.add(chartPanel, constraints);

    constraints.gridy += 1;
    constraints.weighty = 0;
    this.add(copyEmbeddedResp, constraints);

    // place the other UI elements in a single row below the chart
    constraints.gridwidth = 1;
    constraints.weighty = 0.0;
    constraints.weightx = 0.0;
    constraints.anchor = GridBagConstraints.WEST;
    constraints.fill = GridBagConstraints.NONE;
    constraints.gridy += 1;
    constraints.gridx = 0;
    this.add(freqSpaceBox, constraints);

    constraints.gridx += 1;
    constraints.weightx = 1.0;
    constraints.fill = GridBagConstraints.NONE;
    constraints.anchor = GridBagConstraints.CENTER;
    this.add(save, constraints);

    // add an empty panel as a spacer to keep the save button in the center
    constraints.fill = GridBagConstraints.HORIZONTAL;
    constraints.gridx += 1;
    constraints.weightx = 0;
    constraints.anchor = GridBagConstraints.WEST;
    this.add(plotSelection, constraints);

  }

  @Override
  public void actionPerformed(ActionEvent event) {

    super.actionPerformed(event);

    if (event.getSource() == plotSelection) {
      if (!set) {
        return;
      }

      JFreeChart[] charts = new JFreeChart[]{magnitudeChart, argumentChart};
      int index = plotSelection.getSelectedIndex();
      chartPanel.setChart(charts[index]);

      return;
    }

    if (event.getSource() == copyEmbeddedResp) {
      Set<String> respFilenames = InstrumentResponse.parseInstrumentList();

      List<String> names = new ArrayList<>(respFilenames);
      Collections.sort(names);

      JDialog dialog = new JDialog();
      Object result = JOptionPane.showInputDialog(
          dialog,
          "Select a response to copy:",
          "RESP File Selection",
          JOptionPane.PLAIN_MESSAGE,
          null, names.toArray(),
          0);

      // did user cancel operation?
      if (result == null) {
        return; // nothing left to do here, so let's close this out
      }

      String resultStr = (String) result;

      try {
        // copy response file out of embedded set and into responses folder
        File respDir = new File("responses/");
        if (!respDir.exists()) {
          //noinspection ResultOfMethodCallIgnored
          respDir.mkdir();
        }

        ClassLoader cl = ResponsePanel.class.getClassLoader();
        InputStream respStream =
            cl.getResourceAsStream(InstrumentResponse.RESP_DIRECTORY + resultStr);
        Path destinationPath = Paths.get(respDir.getCanonicalPath(), resultStr);
        Files.copy(respStream, destinationPath, REPLACE_EXISTING);

      } catch (IOException e) {
        e.printStackTrace();
      }
    }
  }

  @Override
  protected void drawCharts() {
    plotSelection.setSelectedIndex(0);
    chart = magnitudeChart;
    chartPanel.setChart(chart);
    chartPanel.setMouseZoomable(true);
  }

  @Override
  public String[] getAdditionalReportPages() {
    ResponseExperiment respExp = (ResponseExperiment) expResult;
    InstrumentResponse[] responses = respExp.getResponses();
    String[] pages = new String[responses.length];
    for (int i = 0; i < pages.length; ++i) {
      pages[i] = responses[i].toString();
    }
    return pages;
  }

  @Override
  public JFreeChart[] getCharts() {
    return new JFreeChart[]{magnitudeChart, argumentChart};
  }

  /**
   * Produce the filename of the report generated from this experiment.
   * Since response data is not directly associated with data at a given
   * time, rather than a sensor as a whole, we merely use the current date
   * and the first response used in the experiment.
   *
   * @return String that will be default filename of PDF generated from data
   */
  @Override
  public String getPDFFilename() {

    SimpleDateFormat dateFormatter = ExperimentPanel.DATE_FORMAT.get();
    Calendar calendar = Calendar.getInstance(dateFormatter.getTimeZone());
    // experiment has no time metadata to be associated with it, get time now
    String date = dateFormatter.format(calendar.getTime());

    String test = expType.getName().replace(' ', '_');

    int idx = getIndexOfMainData(); // first resp in list
    String name = expResult.getInputNames().get(idx);

    return test + '_' + name + '_' + date + ".pdf";
  }

  @Override
  public ValueAxis getXAxis() {
    // true if using Hz units
    if (freqSpaceBox.isSelected()) {
      return frequencyAxis;
    }

    return xAxis;
  }

  @Override
  public ValueAxis getYAxis() {
    ValueAxis[] axes = new ValueAxis[]{yAxis, degreeAxis};
    if (null == plotSelection) {
      return yAxis;
    }

    return axes[plotSelection.getSelectedIndex()];
  }

  @Override
  public int panelsNeeded() {
    return 3;
  }

  @Override
  public int plotsToShow() {
    return 0;
  }

  @Override
  protected void updateData(DataStore dataStore) {

    set = true;

    seriesColorMap = new HashMap<>();

    boolean freqSpace = freqSpaceBox.isSelected();
    ResponseExperiment respExp = (ResponseExperiment) expResult;
    respExp.setFreqSpace(freqSpace);
    expResult.runExperimentOnData(dataStore);

    List<XYSeriesCollection> timeSeries = expResult.getData();
    XYSeriesCollection magSeries = timeSeries.get(0);
    XYSeriesCollection argSeries = timeSeries.get(1);

    for (int i = 0; i < magSeries.getSeriesCount(); ++i) {
      Color toColor = ReportingUtils.COLORS[i % ReportingUtils.COLORS.length];
      String magName = (String) magSeries.getSeriesKey(i);
      String argName = (String) argSeries.getSeriesKey(i);
      seriesColorMap.put(magName, toColor);
      seriesColorMap.put(argName, toColor);
    }

    argumentChart = buildChart(argSeries, getXAxis(), degreeAxis);
    argumentChart.getXYPlot().getRangeAxis().setAutoRange(true);
    magnitudeChart = buildChart(magSeries, getXAxis(), yAxis);
    magnitudeChart.getXYPlot().getRangeAxis().setAutoRange(true);

    appendChartTitle(argumentChart, " " + ResponseExperiment.ARGUMENT);
    appendChartTitle(magnitudeChart, " " + ResponseExperiment.MAGNITUDE);
  }

}
