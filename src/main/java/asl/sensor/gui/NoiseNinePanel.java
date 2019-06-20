package asl.sensor.gui;

import asl.sensor.ExperimentFactory;
import asl.sensor.experiment.NoiseNineExperiment;
import asl.sensor.input.DataStore;
import asl.utils.ReportingUtils;
import java.awt.Color;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import javax.swing.BoxLayout;
import javax.swing.ButtonGroup;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.XYPlot;
import org.jfree.data.xy.XYSeriesCollection;

/**
 * Panel for 9-input self noise. Similar to 3-input self noise (NoisePanel)
 * but includes multiple plots, one for each linear axis in 3D space
 * (north-south, east-west, up-down) and a combo box to select them
 *
 * @author akearns - KBRWyle
 */
public class NoiseNinePanel extends NoisePanel {

  private static final long serialVersionUID = -8049021432657749975L;
  private final JComboBox<String> plotSelection;
  private JFreeChart northChart, eastChart, verticalChart;
  private JRadioButton firstRadioButton, secondRadioButton;

  /**
   * Construct panel and lay out its components
   *
   * @param experiment Enum to get relevant experiment backend from factory (NoiseNineExperiment)
   */
  public NoiseNinePanel(ExperimentFactory experiment) {
    super(experiment);

    expResult = experiment.createExperiment();

    set = false;

    for (int i = 0; i < 3; ++i) {
      int num = i + 1;
      channelType[3 * i] = "North sensor " + num + " (RESP required)";
      channelType[(3 * i) + 1] = "East sensor " + num + " (RESP required)";
      channelType[(3 * i) + 2] = "Vertical sensor " + num + " (RESP required)";
    }

    firstRadioButton = new JRadioButton("first ");
    firstRadioButton.setSelected(true);
    secondRadioButton = new JRadioButton("second ");
    // we don't keep thisRadioButton as a global var since it's only selected if others aren't
    JRadioButton thirdRadioButton = new JRadioButton("third ");


    ButtonGroup group = new ButtonGroup();
    group.add(firstRadioButton);
    group.add(secondRadioButton);
    group.add(thirdRadioButton);

    JPanel angleRefSelection = new JPanel();
    angleRefSelection.setLayout(new BoxLayout(angleRefSelection, BoxLayout.X_AXIS));
    angleRefSelection.add(new JLabel("Use the "));
    angleRefSelection.add(firstRadioButton);
    angleRefSelection.add(secondRadioButton);
    angleRefSelection.add(thirdRadioButton);
    angleRefSelection.add(new JLabel("input set for angle reference (requires regen)"));
    angleRefSelection.setMaximumSize(angleRefSelection.getMinimumSize());
    angleRefSelection.setPreferredSize(angleRefSelection.getMinimumSize());

    this.setLayout(new GridBagLayout());
    GridBagConstraints constraints = new GridBagConstraints();

    northChart =
        ChartFactory.createXYLineChart(expType.getName() + " (North)",
            "", "", null);
    eastChart =
        ChartFactory.createXYLineChart(expType.getName() + " (East)",
            "", "", null);
    verticalChart =
        ChartFactory.createXYLineChart(expType.getName() + " (Vertical)",
            "", "", null);
    for (JFreeChart chart : getCharts()) {
      chart.getXYPlot().setDomainAxis(getXAxis());
      chart.getXYPlot().setRangeAxis(getYAxis());
    }

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

    constraints.gridy += 1;
    constraints.weighty = 0;
    constraints.fill = GridBagConstraints.NONE;
    this.add(angleRefSelection, constraints);

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

  /**
   * Get string of data to used when building PDFs in specific circumstances
   *
   * @param experiment Experiment to extract data from
   * @return String representing experiment data (rotation angles)
   */
  public static String getInsetString(NoiseNineExperiment experiment) {

    String[] output = experiment.getInsetStrings();
    StringBuilder sb = new StringBuilder();
    for (String s : output) {
      sb.append(s).append('\n');
    }
    return sb.toString();
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

    if (firstRadioButton.isSelected())  {
      noiseExperiment.setFirstDataAsAngleReference();
    } else if (secondRadioButton.isSelected()) {
      noiseExperiment.setSecondDataAsAngleReference();
    } else {
      noiseExperiment.setThirdDataAsAngleReference();
    }

    expResult.runExperimentOnData(dataStore);

    for (int j = 0; j < 3; ++j) {
      XYSeriesCollection timeseries = expResult.getData().get(j);

      for (int i = 0; i < NOISE_PLOT_COUNT; ++i) {
        String name = (String) timeseries.getSeriesKey(i);
        Color plotColor = ReportingUtils.COLORS[i % 3];
        seriesColorMap.put(name, plotColor);
        if (i >= 3) {
          seriesDashedSet.add(name);
        }

      }
    }

    set = true;

    NoiseNineExperiment nineExperiment = (NoiseNineExperiment) expResult;
    String[] insetStrings = nineExperiment.getInsetStrings();
    XYPlot plot;

    northChart = buildChart(expResult.getData().get(0));
    northChart.setTitle("Self-noise (NORTH)");
    plot = northChart.getXYPlot();
    setTitle(plot, insetStrings[0]);

    eastChart = buildChart(expResult.getData().get(1));
    eastChart.setTitle("Self-noise (EAST)");
    plot = eastChart.getXYPlot();
    setTitle(plot, insetStrings[1]);

    verticalChart = buildChart(expResult.getData().get(2));
    verticalChart.setTitle("Self-noise (VERTICAL)");
    plot = verticalChart.getXYPlot();
    setTitle(plot, insetStrings[2]);

  }

}
