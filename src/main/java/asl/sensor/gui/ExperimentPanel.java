package asl.sensor.gui;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics2D;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JPanel;
import javax.swing.filechooser.FileNameExtensionFilter;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.annotations.XYTitleAnnotation;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.chart.title.TextTitle;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.RectangleAnchor;

import asl.sensor.experiment.Experiment;
import asl.sensor.experiment.ExperimentEnum;
import asl.sensor.experiment.ExperimentFactory;
import asl.sensor.input.DataStore;

/**
 * Panel used to display the data produced from a specified sensor test.
 * This primarily exists as a chartpanel, plus a filechooser and button used
 * in saving the chart held in the panel to a file.
 * 
 * A default construction of the GUI components exists in this class, but
 * implementing methods are suggested to override this in their constructor
 * (see existing classes such as GainPanel and NoisePanel for examples).
 * 
 * @author akearns
 *
 */
public abstract class ExperimentPanel extends JPanel implements ActionListener {
 
  private static final long serialVersionUID = -5591522915365766604L;

  protected JButton save;
  
  protected JFreeChart chart; // replace with plot object
  protected ChartPanel chartPanel;
  
  protected JFileChooser fc; // save image when image save button clicked
  
  public final ExperimentEnum expType; 
          // used to define experiment of each plot object
  
  protected Experiment expResult;
          // used to get the actual data from loaded-in files
  
  // axes and titles must be instantiated in implementing functino
  protected String xAxisTitle, yAxisTitle;
  protected ValueAxis xAxis, yAxis;
  
  public String[] channelType;
  
  protected String[] plotTheseInBold; // given in the implementing function
  // this is a String because bolded names are intended to be fixed
  
  protected Map<String, Color> seriesColorMap;
  protected Set<String> seriesDashedSet;
  // these are map/set because they are based on the data read in, not fixed
  
  public ExperimentPanel(ExperimentEnum exp) {
    
    channelType = new String[DataStore.FILE_COUNT];
    
    for (int i = 0; i < channelType.length; ++i) {
      channelType[i] = "NOT USED";
    }
    
    seriesColorMap = new HashMap<String, Color>();
    seriesDashedSet = new HashSet<String>();
    plotTheseInBold = new String[]{};
    
    expType = exp;
    expResult = ExperimentFactory.createExperiment(exp);
    
    chart = ChartFactory.createXYLineChart( expType.getName(), 
        "", "", null);
    chartPanel = new ChartPanel(chart);
    // chartPanel.setMouseZoomable(false);
    
    fc = new JFileChooser();
    
    save = new JButton("Save plot (PNG)");
    save.addActionListener(this);
    
    // basic layout for components (recommended to override in concrete class)
    this.setLayout( new BoxLayout(this, BoxLayout.Y_AXIS) );
    this.add(chartPanel);
    this.add(save);
  }
  
  public String[] getChannelTypes() {
    return channelType;
  }
  
  /**
   * Handle's saving this plot's chart to file (PNG image) 
   * when the save button is clicked.
   */
  @Override
  public void actionPerformed(ActionEvent e) {
    
    if( e.getSource() == save ) {
      String ext = ".png";
      fc.addChoosableFileFilter(
          new FileNameExtensionFilter("PNG image (.png)",ext) );
      fc.setFileFilter(fc.getChoosableFileFilters()[1]);
      int returnVal = fc.showSaveDialog(save);
      if (returnVal == JFileChooser.APPROVE_OPTION) {
        File selFile = fc.getSelectedFile();
        if( !selFile.getName().endsWith( ext.toLowerCase() ) ) {
          selFile = new File( selFile.toString() + ext);
        }
        try {
          ChartUtilities.saveChartAsPNG(selFile,chart,640,480);
        } catch (IOException e1) {
          e1.printStackTrace();
        }
      }
    }
  }
  
  /**
   * Gets the axes to be used to plot the data 
   */
  protected void applyAxesToChart() {
    XYPlot xyp = chart.getXYPlot();
    xyp.setDomainAxis( getXAxis() );
    xyp.setRangeAxis( getYAxis() );
  }
  
  /**
   * Overlay an error message in the event of an exception or other issue
   * @param errMsg Text of the message to be displayed
   */
  public void displayErrorMessage(String errMsg) {
    XYPlot xyp = (XYPlot) chartPanel.getChart().getPlot();
    TextTitle result = new TextTitle();
    result.setText(errMsg);
    result.setBackgroundPaint(Color.red);
    result.setPaint(Color.white);
    XYTitleAnnotation xyt = new XYTitleAnnotation(0.5, 0.5, result,
        RectangleAnchor.CENTER);
    xyp.clearAnnotations();
    xyp.addAnnotation(xyt);
  }
  
  /**
   * Overlay informational text, such as extra results and statistics for plots
   * @param infoMsg
   */
  public void displayInfoMessage(String infoMsg) {
    XYPlot xyp = (XYPlot) chartPanel.getChart().getPlot();
    TextTitle result = new TextTitle();
    result.setText(infoMsg);
    XYTitleAnnotation xyt = new XYTitleAnnotation(0.5, 0.5, result,
        RectangleAnchor.CENTER);
    xyp.clearAnnotations();
    xyp.addAnnotation(xyt);
  }
  
  /**
   * Return image of this chart with specified dimensions
   * Used to compile PNG image of all currently-displayed charts
   * @param width Width of output image in pixels
   * @param height Height of output image in pixels
   * @return buffered image of this panel's chart
   */
  public BufferedImage getAsImage(int width, int height) {
    
    ChartPanel outPanel = new ChartPanel(chart);
    outPanel.setSize( new Dimension(width, height) );
    
    BufferedImage bi = new BufferedImage(
        width, 
        height, 
        BufferedImage.TYPE_INT_ARGB);

    Graphics2D g = bi.createGraphics();
    outPanel.printAll(g);
    g.dispose();

    return bi;
  }

  /**
   * Default x-axis return function.
   * Though the x-axis is a local variable, some panels may have multiple unit
   * types for the x-axis (i.e., for units of seconds vs. Hz); accessing
   * the x-axis object through this function allows for overrides allowing for
   * more flexibility.
   * @return
   */
  public ValueAxis getXAxis() {
    return xAxis;
  }
  

  /**
   * Default x-axis title return. Displays the string used for the x-axis, 
   * which is set when the panel's chart is constructed.
   * @return
   */
  public String getXTitle() {
    return xAxisTitle;
  }

  /**
   * Default y-axis return function. As with x-axis, designed to be overridden
   * for charts that may use multiple scales
   * @return
   */
  public ValueAxis getYAxis() {
    return yAxis;
  }


  public String getYTitle() {
    return yAxisTitle;
  }

  /**
   * Used to plot the results of a backend function from an experiment
   * using a collection of XYSeries mapped by strings.
   * If the color and dashed maps have been updated by the concrete class
   * extending this one, the corresponding XYSeries will have their data
   * displayed in the specified color or dashing, etc.
   * @param xyDataset collection of XYSeries to plot
   */
  protected void populateChart(XYDataset xyDataset) {
    
    JFreeChart chart = ChartFactory.createXYLineChart(
        expType.getName(),
        getXTitle(),
        getYTitle(),
        xyDataset,
        PlotOrientation.VERTICAL,
        true, // include legend
        false, 
        false);
    
    // apply effects to the components that require it (i.e., NLNM time series)
    XYPlot xyPlot = chart.getXYPlot();
    XYItemRenderer xyir = xyPlot.getRenderer();
    
    if (xyDataset instanceof XYSeriesCollection) {
      // force certain colors and whether or not a line should be dashed
      XYSeriesCollection xysc = (XYSeriesCollection) xyDataset;
      
      for ( String series : seriesColorMap.keySet() ) {
        int seriesIdx = xysc.getSeriesIndex(series);
        if (seriesIdx >= 0) {
          xyir.setSeriesPaint( seriesIdx, seriesColorMap.get(series) );
        }
        
        if ( seriesDashedSet.contains(series) ) {
          xyir.setSeriesPaint( seriesIdx, seriesColorMap.get(series).darker() );
          
          BasicStroke stroke = (BasicStroke) xyir.getSeriesStroke(seriesIdx);
          if (stroke == null) {
            stroke = (BasicStroke) xyir.getBaseStroke();
          }
          float width = stroke.getLineWidth();
          int join = stroke.getLineJoin();
          int cap = stroke.getEndCap();
          
          float[] dashing = new float[]{1,4};
          
          stroke = new BasicStroke(width, cap, join, 10f, dashing, 0f);
          xyir.setSeriesStroke(seriesIdx, stroke);
        }
        
      }

      if ( !(plotTheseInBold.length == 0) ) {
        for (String series : plotTheseInBold) {
          int seriesIdx = xysc.getSeriesIndex(series);
          
          BasicStroke stroke = (BasicStroke) xyir.getSeriesStroke(seriesIdx);
          if (stroke == null) {
            stroke = (BasicStroke) xyir.getBaseStroke();
          }
          stroke = new BasicStroke( stroke.getLineWidth()*2 );
          xyir.setSeriesStroke(seriesIdx, stroke);
          xyir.setSeriesPaint(seriesIdx, new Color(0,0,0) );
        }
      }
    }
    

    
    xyPlot.setDomainAxis( getXAxis() );
    xyPlot.setRangeAxis( getYAxis() );
    
    this.chart = chart;
  }
  
  /**
   * Function template for sending input to a backend fucntion and collecting
   * the corresponding data
   * @param ds DataStore object containing seed and resp files
   */
  public abstract void updateData(final DataStore ds);
  
  public abstract int panelsNeeded();
  
  public boolean hasEnoughData(final DataStore ds) {
    return expResult.hasEnoughData(ds);
  }
  
  // details of how to run updateData are left up to the implementing panel
  // however, it is advised to wrap the code inside a swingworker,
  // using doInBackground() to make calls to the experiment backend to
  // compute the data and then using the done() method to update chartPanel

}
