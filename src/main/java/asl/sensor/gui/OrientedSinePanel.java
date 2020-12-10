package asl.sensor.gui;

import asl.sensor.ExperimentFactory;
import asl.sensor.experiment.OrientedSineExperiment;
import asl.sensor.input.DataStore;
import java.awt.Color;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import javax.swing.BoxLayout;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.event.ChangeEvent;
import org.jfree.chart.annotations.XYTitleAnnotation;
import org.jfree.chart.axis.DateAxis;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.title.TextTitle;
import org.jfree.chart.ui.RectangleAnchor;
import org.jfree.data.xy.XYSeriesCollection;

public class OrientedSinePanel extends ExperimentPanel {

  private final JCheckBox requiresRotation;
  private final JComboBox<String> rotationTypeSelector;
  private static final String TRILLIUM = "Trillium";

  /**
   * Construct a new panel, using a backend defined by the passed-in enum
   *
   * @param experiment Experiment enum with corresponding backend for factory instantiation
   */
  public OrientedSinePanel(ExperimentFactory experiment) {
    super(experiment);
    channelType[0] = "North calibration input";
    channelType[1] = "East calibration input";
    channelType[2] = "Vertical calibration input";

    xAxis = new DateAxis("Time (UTC)");
    xAxis.setAutoRange(true);
    ((DateAxis) xAxis).setDateFormatOverride(ExperimentPanel.DATE_TIME_FORMAT.get());
    yAxis = new NumberAxis("Mean-removed sine wave signal from data");
    yAxis.setAutoRange(true);
    Font bold = xAxis.getLabelFont();
    bold = bold.deriveFont(Font.BOLD, bold.getSize() + 2);
    xAxis.setLabelFont(bold);
    yAxis.setLabelFont(bold);
    applyAxesToChart();

    JPanel rotationSubpanel = new JPanel();
    rotationSubpanel.setLayout(new BoxLayout(rotationSubpanel, BoxLayout.X_AXIS));
    requiresRotation = new JCheckBox();
    requiresRotation.setEnabled(true);
    requiresRotation.setSelected(true);
    requiresRotation.addChangeListener(this);
    rotationSubpanel.add(requiresRotation);
    rotationTypeSelector = new JComboBox<>(new String[]{"STS", TRILLIUM});
    rotationSubpanel.setEnabled(true);
    rotationSubpanel.add(new JLabel("Enable rotation and use rotation model"));
    rotationSubpanel.add(rotationTypeSelector);


    this.setLayout(new GridBagLayout());
    GridBagConstraints constraints = new GridBagConstraints();
    constraints.weighty = 1;
    constraints.weightx = 1;
    constraints.gridx = 0;
    constraints.gridy = 0;
    constraints.anchor = GridBagConstraints.CENTER;
    constraints.fill = GridBagConstraints.BOTH;
    this.add(chartPanel, constraints);
    constraints.fill = GridBagConstraints.NONE;
    constraints.weightx = 0;
    constraints.weighty = 0;
    ++constraints.gridy;
    this.add(rotationSubpanel, constraints);
    ++constraints.gridy;
    this.add(save, constraints);

  }

  public void stateChanged(ChangeEvent e) {
    if (e.getSource() == requiresRotation) {
      rotationTypeSelector.setEnabled(requiresRotation.isSelected());
    }
  }

  @Override
  protected void drawCharts() {
    OrientedSineExperiment ose = (OrientedSineExperiment) expResult;
    chart = buildChart(expResult.getData().get(0));
    XYPlot plot = (XYPlot) chart.getPlot();

    TextTitle result = getDefaultTextTitle();
    if (ose.hasError()) {
      result.setBackgroundPaint(Color.red);
      result.setPaint(Color.white);
    }
    result.setText(expResult.getReportString());
    XYTitleAnnotation title = new XYTitleAnnotation(0.98, 0.98, result,
        RectangleAnchor.TOP_RIGHT);
    plot.clearAnnotations();
    plot.addAnnotation(title);

    chartPanel.setChart(chart);
    chartPanel.setMouseZoomable(true);
  }

  @Override
  public int panelsNeeded() {
    return expResult.blocksNeeded();
  }

  @Override
  protected void updateData(DataStore dataStore) {
    OrientedSineExperiment ose = (OrientedSineExperiment) expResult;
    boolean rotate = requiresRotation.isSelected();
    ose.setDoRotation(rotate);
    if (rotate) {
      String selection = (String) rotationTypeSelector.getSelectedItem();
      assert (selection != null); // why have a null in the menu when the checkbox is there?
      ose.setRotationTrillium(selection.equals(TRILLIUM));
    }

    ose.runExperimentOnData(dataStore);

    XYSeriesCollection timeseries = ose.getData().get(0);
    for (int i = 0; i < timeseries.getSeriesCount(); ++i) {
      Color toColor = getColor(i);
      String curve = (String) timeseries.getSeriesKey(i);
      seriesColorMap.put(curve, toColor);
    }
    set = true;
  }
}
