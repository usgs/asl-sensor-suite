package asl.sensor.gui;

import static org.junit.Assert.assertNotNull;

import asl.sensor.experiment.ExperimentEnum;
import asl.sensor.input.DataStore;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.junit.Test;

public class ExperimentPanelTest {

  @Test
  public void buildChart_nullxyDataset_DefaultAxis() {
    ExperimentPanel panel = new MockPanel(ExperimentEnum.AZMTH);
    //Since nothing else is initialized, this should just return without Exception.
    JFreeChart chart = panel.buildChart(null);
    assertNotNull(chart);
  }

  @Test
  public void buildChart_nullxyDataset_AxisDefined() {
    ExperimentPanel panel = new MockPanel(ExperimentEnum.AZMTH);
    //Since nothing else is initialized, this should just return without an Exception.
    JFreeChart chart = panel.buildChart(
        null, new NumberAxis("x Custom"), new NumberAxis("y Custom"));
    assertNotNull(chart);
  }


  class MockPanel extends ExperimentPanel {

    MockPanel(ExperimentEnum exp) {
      super(exp);
      xAxis = new NumberAxis("x Axis");
      yAxis = new NumberAxis("y Axis");
    }

    @Override
    protected void drawCharts() {

    }

    @Override
    public int panelsNeeded() {
      return 0;
    }

    @Override
    protected void updateData(DataStore ds) {

    }
  }
}