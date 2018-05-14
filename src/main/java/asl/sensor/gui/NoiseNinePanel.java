package asl.sensor.gui;

import asl.sensor.ExperimentFactory;
import asl.sensor.experiment.NoiseNineExperiment;
import asl.sensor.input.DataStore;
import java.awt.Color;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import javax.swing.JComboBox;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.annotations.XYTitleAnnotation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.title.TextTitle;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.RectangleAnchor;

/**
 * Panel for 9-input self noise. Similar to 3-input self noise (NoisePanel)
 * but includes multiple plots, one for each linear axis in 3D space
 * (north-south, east-west, up-down) and a combo box to select them
 *
 * @author akearns - KBRWyle
 */
public class NoiseNinePanel extends NoisePanel {

  private static final long serialVersionUID = -8049021432657749975L;

  /**
   * Get string of data to used when building PDFs in specific circumstances
   *
   * @param experiment Experiment to extract data from
   * @return String representing experiment data (rotation angles)
   */
  public static String getInsetString(NoiseNineExperiment experiment) {
    double[] angles = experiment.getNorthAngles();
    StringBuilder sb = new StringBuilder();
    sb.append("Angle of rotation of north sensor 2 (deg): ");
    sb.append(DECIMAL_FORMAT.get().format(Math.toDegrees(angles[0])));
    sb.append("\nAngle of rotation of north sensor 3 (deg): ");
    sb.append(DECIMAL_FORMAT.get().format(Math.toDegrees(angles[1])));
    sb.append("\n");
    angles = experiment.getEastAngles();
    sb.append("Angle of rotation of east sensor 2 (deg): ");
    sb.append(DECIMAL_FORMAT.get().format(Math.toDegrees(angles[0])));
    sb.append("\nAngle of rotation of east sensor 3 (deg): ");
    sb.append(DECIMAL_FORMAT.get().format(Math.toDegrees(angles[1])));
    return sb.toString();
  }

  private final JComboBox<String> plotSelection;

  private JFreeChart northChart, eastChart, verticalChart;

  /**
   * Construct panel and lay out its components
   *
   * @param experiment Enum to get relevant experiment backend from factory (NoiseNineExperiment)
   */
  NoiseNinePanel(ExperimentFactory experiment) {
    super(experiment);

    expResult = experiment.createExperiment();

    set = false;

    for (int i = 0; i < 3; ++i) {
      int num = i + 1;
      channelType[3 * i] = "North sensor " + num + " (RESP required)";
      channelType[(3 * i) + 1] = "East sensor " + num + " (RESP required)";
      channelType[(3 * i) + 2] = "Vertical sensor " + num + " (RESP required)";
    }

    this.setLayout(new GridBagLayout());
    GridBagConstraints constraints = new GridBagConstraints();

    String xTitle = getXAxis().getLabel();
    String yTitle = getYAxis().getLabel();

    northChart =
        ChartFactory.createXYLineChart(expType.getName() + " (North)",
            xTitle, yTitle, null);
    eastChart =
        ChartFactory.createXYLineChart(expType.getName() + " (East)",
            xTitle, yTitle, null);
    verticalChart =
        ChartFactory.createXYLineChart(expType.getName() + " (Vertical)",
            xTitle, yTitle, null);

    chart = northChart;
    chartPanel.setChart(chart);

    removeAll(); // get rid of blank spacer JPanel from super
    // (everything else will get redrawn)

    constraints.fill = GridBagConstraints.BOTH;
    constraints.gridx = 0;
    constraints.gridy = 0;
    constraints.weightx = 1.0;
    constraints.weighty = 1.0;
    constraints.gridwidth = 3;
    constraints.anchor = GridBagConstraints.CENTER;
    add(chartPanel, constraints);

    // place the other UI elements in a single row below the chart
    constraints.gridwidth = 1;
    constraints.weighty = 0.0;
    constraints.weightx = 0.0;
    constraints.anchor = GridBagConstraints.WEST;
    constraints.fill = GridBagConstraints.NONE;
    constraints.gridy += 1;
    constraints.gridx = 0;
    add(freqSpaceBox, constraints);

    constraints.gridx += 1;
    constraints.weightx = 1.0;
    constraints.fill = GridBagConstraints.NONE;
    constraints.anchor = GridBagConstraints.CENTER;
    add(save, constraints);

    // combo box to select items
    constraints.fill = GridBagConstraints.HORIZONTAL;
    constraints.gridx += 1;
    constraints.weightx = 0;
    constraints.anchor = GridBagConstraints.WEST;
    plotSelection = new JComboBox<>();
    plotSelection.addItem("North component plot");
    plotSelection.addItem("East component plot");
    plotSelection.addItem("Vertical component plot");
    plotSelection.addActionListener(this);
    add(plotSelection, constraints);

    revalidate();

  }

  @Override
  public void actionPerformed(ActionEvent event) {

    super.actionPerformed(event);

    if (event.getSource() == plotSelection) {

      int index = plotSelection.getSelectedIndex();
      switch (index) {
        case 0:
          chart = northChart;
          break;
        case 1:
          chart = eastChart;
          break;
        default:
          chart = verticalChart;
          break;
      }

      chartPanel.setChart(chart);
      chartPanel.setMouseZoomable(true);
    }

  }

  @Override
  protected void drawCharts() {

    int index = plotSelection.getSelectedIndex();

    switch (index) {
      case 0:
        chart = northChart;
        break;
      case 1:
        chart = eastChart;
        break;
      default:
        chart = verticalChart;
        break;
    }

    chartPanel.setChart(chart);
    chartPanel.setMouseZoomable(true);

  }

  @Override
  public JFreeChart[] getCharts() {
    return new JFreeChart[]{northChart, eastChart, verticalChart};
  }

  /**
   * Get text representation of angles used to rotate data
   *
   * @return String displaying angles of rotation for the 2nd, 3rd east sensors
   */
  private String getEastChartString() {
    NoiseNineExperiment experiment = (NoiseNineExperiment) expResult;
    double[] angles = experiment.getEastAngles();
    return "Angle of rotation of east sensor 2 (deg): "
        + DECIMAL_FORMAT.get().format(Math.toDegrees(angles[0]))
        + "\nAngle of rotation of east sensor 3 (deg): "
        + DECIMAL_FORMAT.get().format(Math.toDegrees(angles[1]));
  }

  @Override
  public String getInsetStrings() {
    return getNorthChartString() + "\n"
        + getEastChartString() + "\n";
  }

  /**
   * Get text representation of angles used to rotate data
   *
   * @return String displaying angles of rotation for the 2nd, 3rd north sensors
   */
  private String getNorthChartString() {
    NoiseNineExperiment nne = (NoiseNineExperiment) expResult;
    double[] angles = nne.getNorthAngles();
    return "Angle of rotation of north sensor 2 (deg): "
        + DECIMAL_FORMAT.get().format(Math.toDegrees(angles[0]))
        + "\nAngle of rotation of north sensor 3 (deg): "
        + DECIMAL_FORMAT.get().format(Math.toDegrees(angles[1]));
  }

  @Override
  public int panelsNeeded() {
    return 9;
  }

  @Override
  public void updateData(final DataStore dataStore) {

    set = true;

    final boolean freqSpaceImmutable = freqSpaceBox.isSelected();

    NoiseNineExperiment noiseExperiment = (NoiseNineExperiment) expResult;
    noiseExperiment.setFreqSpace(freqSpaceImmutable);

    expResult.runExperimentOnData(dataStore);

    for (int j = 0; j < 3; ++j) {
      XYSeriesCollection timeseries = expResult.getData().get(j);

      for (int i = 0; i < NOISE_PLOT_COUNT; ++i) {
        String name = (String) timeseries.getSeriesKey(i);
        System.out.println(name);
        Color plotColor = COLORS[i % 3];
        seriesColorMap.put(name, plotColor);
        System.out.println(name + "," + plotColor);
        if (i >= 3) {
          System.out.println(name + "," + i);
          seriesDashedSet.add(name);
        }

      }
    }

    set = true;

    northChart = buildChart(expResult.getData().get(0));
    northChart.setTitle("Self-noise (NORTH)");
    XYPlot plot = northChart.getXYPlot();
    TextTitle angle = new TextTitle();

    angle.setText(getNorthChartString());
    angle.setBackgroundPaint(Color.white);
    XYTitleAnnotation title = new XYTitleAnnotation(0.98, 0.98, angle,
        RectangleAnchor.TOP_RIGHT);
    plot.clearAnnotations();
    plot.addAnnotation(title);

    eastChart = buildChart(expResult.getData().get(1));
    eastChart.setTitle("Self-noise (EAST)");
    plot = eastChart.getXYPlot();
    angle = new TextTitle();
    angle.setText(getEastChartString());
    angle.setBackgroundPaint(Color.white);
    title = new XYTitleAnnotation(0.98, 0.98, angle, RectangleAnchor.TOP_RIGHT);
    plot.clearAnnotations();
    plot.addAnnotation(title);

    verticalChart = buildChart(expResult.getData().get(2));
    verticalChart.setTitle("Self-noise (VERTICAL)");

  }

}
