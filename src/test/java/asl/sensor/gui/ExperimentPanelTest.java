package asl.sensor.gui;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import asl.sensor.ExperimentFactory;
import asl.sensor.input.DataStore;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.junit.Before;
import org.junit.Test;

public class ExperimentPanelTest {

  ExperimentPanel panel;

  @Before
  public void setUp() {
    this.panel = new MockPanel(ExperimentFactory.AZIMUTH);
  }

  @Test
  public void buildChart_nullxyDataset_DefaultAxis() {
    //Since nothing else is initialized, this should just return without Exception.
    JFreeChart chart = panel.buildChart(null);
    assertNotNull(chart);
  }

  @Test
  public void buildChart_nullxyDataset_AxisDefined() {
    //Since nothing else is initialized, this should just return without an Exception.
    JFreeChart chart = panel.buildChart(
        null, new NumberAxis("x Custom"), new NumberAxis("y Custom"));
    assertNotNull(chart);
  }

  @Test
  public void setChart_nullxyDataset() {
    //Since nothing else is initialized, this should just return without Exception.
    panel.setChart(null);
    assertNotNull(panel.chart);
  }

  @Test
  public void hasRun_defaultValueIsFalse() {
    assertFalse(panel.hasRun());
  }

  @Test
  public void clearChart_setsHasRunToFalse() {
    panel.updateData(null);
    assertTrue(panel.hasRun());
    panel.clearChart();
    assertFalse(panel.hasRun());
  }

  @Test
  public void getXAxis_returnsDefault() {
    assertEquals("x Axis", panel.getXAxis().getLabel());
  }

  @Test
  public void getYAxis_returnsDefault() {
    assertEquals("y Axis", panel.getYAxis().getLabel());
  }

  @Test
  public void getPDFFilename() {
    assertTrue(panel.getPDFFilename().startsWith("Azimuth__"));
    assertTrue(panel.getPDFFilename().endsWith(".pdf"));
  }

  @Test
  public void plotsToShow_DefaultToPanelsNeeded() {
    //This since plotsToShow isn't overridden, it should match PanelsNeeded().
    assertEquals(panel.panelsNeeded(), panel.plotsToShow());
  }


  class MockPanel extends ExperimentPanel {

    MockPanel(ExperimentFactory exp) {
      super(exp);
      xAxis = new NumberAxis("x Axis");
      yAxis = new NumberAxis("y Axis");
    }

    @Override
    protected void drawCharts() {

    }

    @Override
    public int panelsNeeded() {
      return 61;
    }

    @Override
    protected void updateData(DataStore dataStore) {
      set = true;
    }
  }
}