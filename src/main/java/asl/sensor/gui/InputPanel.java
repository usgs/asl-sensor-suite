package asl.sensor.gui;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.TimeZone;

import javax.imageio.ImageIO;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSlider;
import javax.swing.JTextField;
import javax.swing.ScrollPaneConstants;
import javax.swing.SwingConstants;
import javax.swing.SwingWorker;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.filechooser.FileNameExtensionFilter;
import javax.swing.text.JTextComponent;

import org.jfree.chart.ChartColor;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.annotations.XYTitleAnnotation;
import org.jfree.chart.axis.DateAxis;
import org.jfree.chart.plot.Marker;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.ValueMarker;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.title.TextTitle;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.RectangleAnchor;

import asl.sensor.input.DataBlock;
import asl.sensor.input.DataStore;
import asl.sensor.input.InstrumentResponse;
import asl.sensor.utils.ReportingUtils;
import asl.sensor.utils.TimeSeriesUtils;


/**
 * Panel used to hold the plots for the files taken in as input
 * Handles the UI for loading SEED and RESP files, plotting their data
 * and selecting regions of that data, as well as outputting the plots to PNG
 * The data exists in two forms: one is an otherwise-hidden DataStore object
 * that is used to hold in the loaded data directly, and the other one is also
 * a DataStore object but has its data passed to the panel's input plots and
 * to experiment calculations.
 * @author akearns
 *
 */
public class InputPanel 
extends JPanel 
implements ActionListener, ChangeListener {

  /**
   * auto-generated serialization UID
   */
  private static final long serialVersionUID = -7302813951637543526L;
  
  public static final SimpleDateFormat SDF = 
      new SimpleDateFormat("Y.DDD.HH:mm");
  
  /**
   * Default height of image produced by the save-as-image function
   * (each chart is 240 pixels tall)
   */
  public static final int IMAGE_HEIGHT = 240;
  public static final int IMAGE_WIDTH = 480;
  
  public static final int MAX_UNSCROLLED = 4;
  
  public static final int FILE_COUNT = DataStore.FILE_COUNT;
  
  private static final int MARGIN = 10; // min space of the two sliders
  public static final int SLIDER_MAX = 10000;
  
  /**
   * Gets the value of start or end time from slider value and DataBlock
   * @param db DataBlock corresponding to one of the plots
   * @param sliderValue Value of starting or ending time slider [0-SLIDER_MAX]
   * @return Long that represents start or end time matching slider's value
   */
  public static long getMarkerLocation(DataBlock db, int sliderValue) {
    long start = db.getStartTime();
    long len = (db.getInterval()) * db.size();
    return start + (sliderValue * len) / SLIDER_MAX;
  }
  
  private int activePlots = FILE_COUNT; // how much data is being displayed
  
  private DataStore ds;
  private DataStore zooms;
  private ChartPanel[] chartPanels = new ChartPanel[FILE_COUNT];
  private Color[] defaultColor = {
          ChartColor.LIGHT_RED, 
          ChartColor.LIGHT_BLUE, 
          ChartColor.LIGHT_GREEN };
  private JButton save;
  private JButton zoomIn;
  private JButton zoomOut;
  private JButton clearAll;
  private JFileChooser fc;
  private JPanel allCharts; // parent of the chartpanels, used for image saving
  private JSlider leftSlider;
  private JSlider rightSlider;
  private JScrollPane inputScrollPane;
      
  private JButton[] seedLoaders  = new JButton[FILE_COUNT];
  private JTextComponent[] seedFileNames = 
      new JTextComponent[FILE_COUNT];
  private JButton[] respLoaders  = new JButton[FILE_COUNT];
  private JTextComponent[] respFileNames = 
      new JTextComponent[FILE_COUNT];
  private JButton[] clearButton = new JButton[FILE_COUNT];
  
  private JLabel[] channelType = new JLabel[FILE_COUNT];
  
  private JPanel[] chartSubpanels = new JPanel[FILE_COUNT];
  
  // used to store current directory locations
  private String seedDirectory = "data";
  private String respDirectory = "responses";

  
  private String saveDirectory = System.getProperty("user.home");
  
  /**
   * Used to construct the panels for loading and displaying SEED data
   * (as well as the corresponding response file)
   * @param i Index of panel to be created, for getting references to the chart
   * panel and the appropriate actionlisteners for the loaders
   * @return composite panel of chart, loaders, and clear button
   */
  private JPanel makeChartSubpanel(int i) {
    
    JPanel chartSubpanel = new JPanel();
    chartSubpanel.setLayout( new GridBagLayout() );
    GridBagConstraints gbc = new GridBagConstraints();
    
    gbc.weightx = 1.0;
    gbc.weighty = 1.0;
    gbc.gridy = 0;
    gbc.anchor = GridBagConstraints.CENTER;
    
    instantiateChart(i);
    
    /*
    chartPanels[i].setMaximumSize(
        new Dimension(Integer.MAX_VALUE, Integer.MAX_VALUE) );
    */
    
    channelType[i] = new JLabel("");
    
    chartPanels[i].setMouseZoomable(true);
    
    seedLoaders[i] = new JButton( "Load SEED file " + (i+1) );
    seedLoaders[i].addActionListener(this);
    seedLoaders[i].setMaximumSize( seedLoaders[i].getMinimumSize() );
    
    JTextField text = new JTextField( "NO FILE LOADED" );
    text.setHorizontalAlignment(SwingConstants.CENTER);
    text.setMaximumSize( text.getPreferredSize() );
    seedFileNames[i] = text;
    seedFileNames[i].setEditable(false);
   
    respLoaders[i] = new JButton( "Load RESP file " + (i+1) );
    respLoaders[i].addActionListener(this);
    respLoaders[i].setMaximumSize( respLoaders[i].getMinimumSize() );
    
    text = new JTextField( "NO FILE LOADED" );
    text.setHorizontalAlignment(SwingConstants.CENTER);
    text.setMaximumSize( text.getPreferredSize() );
    respFileNames[i] = text;
    respFileNames[i].setEditable(false);
    
    clearButton[i] = new JButton( "Clear data " + (i+1) );
    clearButton[i].setMaximumSize( clearButton[i].getMinimumSize() );
    
    gbc.gridx = 0; gbc.gridy = 0;

    gbc.weightx = 0; gbc.weighty = 0;
    gbc.fill = GridBagConstraints.BOTH;
    
    gbc.fill = GridBagConstraints.NONE;
    chartSubpanel.add(channelType[i], gbc);
    
    gbc.weightx = 1; gbc.weighty = 1;
    gbc.fill = GridBagConstraints.BOTH;
    gbc.gridwidth = 1; gbc.gridheight = 5;
    gbc.gridy += 1;
    chartSubpanel.add(chartPanels[i], gbc);
    

    // Temoved a line to resize the chartpanels
    // This made sense before switching to gridbaglayout, but since that
    // tries to fill space with whatever panels it can, we can just get rid
    // of the code to do that. This also fixes the issue with the text boxes
    // resizing -- just let the charts fill as much space as they can instead
    
    gbc.fill = GridBagConstraints.BOTH;
    gbc.gridx = 1;
    gbc.gridwidth = 1; gbc.gridheight = 1;
    gbc.weightx = 0; gbc.weighty = 0.25;
    chartSubpanel.add(seedLoaders[i], gbc);
    
    gbc.fill = GridBagConstraints.BOTH;
    gbc.weighty = 1;
    gbc.gridy += 1;
    JScrollPane jsp = new JScrollPane();
    jsp.setMaximumSize( jsp.getMinimumSize() );
    jsp.setViewportView(seedFileNames[i]);
    jsp.setVerticalScrollBarPolicy(
        ScrollPaneConstants.VERTICAL_SCROLLBAR_NEVER);
    jsp.setHorizontalScrollBarPolicy(
        ScrollPaneConstants.HORIZONTAL_SCROLLBAR_ALWAYS);
    chartSubpanel.add(jsp, gbc);
    
    gbc.fill = GridBagConstraints.BOTH;
    gbc.weighty = 0.25;
    gbc.gridy += 1;
    chartSubpanel.add(respLoaders[i], gbc);

    
    gbc.fill = GridBagConstraints.BOTH;
    gbc.weighty = 1;
    gbc.gridy += 1;
    jsp = new JScrollPane();
    jsp.setMaximumSize( jsp.getMinimumSize() );
    jsp.setViewportView(respFileNames[i]);
    jsp.setVerticalScrollBarPolicy(
        ScrollPaneConstants.VERTICAL_SCROLLBAR_NEVER);
    jsp.setHorizontalScrollBarPolicy(
        ScrollPaneConstants.HORIZONTAL_SCROLLBAR_ALWAYS);
    chartSubpanel.add(jsp, gbc);
    
    gbc.fill = GridBagConstraints.HORIZONTAL;
    gbc.weighty = 0;
    gbc.gridy += 1;
    clearButton[i] = new JButton( "Clear data " + (i+1) );
    clearButton[i].setOpaque(true);
    clearButton[i].setBackground( Color.RED.darker() );
    clearButton[i].addActionListener(this);
    clearButton[i].setEnabled(false);
    chartSubpanel.add(clearButton[i], gbc);
    
    return chartSubpanel;
  }
  
  
  /**
   * Creates a new data panel -- instantiates each chart, to be populated with
   * data when a file is loaded in. Also creates a save button for writing all
   * the inputted data plots into a single PNG file.
   */
  public InputPanel() {
    
    this.setLayout( new GridBagLayout() );
    GridBagConstraints gbc = new GridBagConstraints();
   
    ds = new DataStore();
    zooms = new DataStore();
    
    gbc.weightx = 1.0;
    gbc.weighty = 1.0;
    gbc.gridy = 0;
    gbc.gridwidth = 8;
    gbc.anchor = GridBagConstraints.CENTER;
    gbc.fill = GridBagConstraints.BOTH;
    
    inputScrollPane = new JScrollPane();
    inputScrollPane.setVerticalScrollBarPolicy(
        ScrollPaneConstants.VERTICAL_SCROLLBAR_AS_NEEDED);
    inputScrollPane.setHorizontalScrollBarPolicy(
        ScrollPaneConstants.HORIZONTAL_SCROLLBAR_AS_NEEDED);
    
    // 
    VScrollPanel cont = new VScrollPanel();
    cont.setLayout( new GridBagLayout() );
    GridBagConstraints contConstraints = new GridBagConstraints();
    contConstraints.weightx = 1.0;
    contConstraints.weighty = 1.0;
    contConstraints.gridy = 0;
    contConstraints.anchor = GridBagConstraints.CENTER;
    contConstraints.fill = GridBagConstraints.BOTH;
    Dimension minDim = new Dimension(0, 0);
    
    for (int i = 0; i < FILE_COUNT; ++i) {
      Dimension d = cont.getPreferredSize();
      chartSubpanels[i] = makeChartSubpanel(i);
      Dimension d2 = chartSubpanels[i].getPreferredSize();
      minDim = chartSubpanels[i].getMinimumSize();
      
      d2.setSize( d2.getWidth(), d.getHeight() );
      
      chartSubpanels[i].setSize( d2 );
      cont.add(chartSubpanels[i], contConstraints);
      contConstraints.gridy += 1;
      // gbc.gridy += 1;
      
    }
    
    inputScrollPane.getViewport().setView(cont);
    inputScrollPane.setVisible(true);
    inputScrollPane.setMinimumSize(minDim);
    
    this.add(inputScrollPane, gbc);
    gbc.gridy += 1;
    
    // set size so that the result pane isn't distorted on window launch
    Dimension d = inputScrollPane.getPreferredSize();
    d.setSize( d.getWidth() + 5, d.getHeight() );
    this.setPreferredSize(d);
    
    leftSlider = new JSlider(0, SLIDER_MAX, 0);
    leftSlider.setEnabled(false);
    leftSlider.addChangeListener(this);
    
    rightSlider = new JSlider(0, SLIDER_MAX, SLIDER_MAX);
    rightSlider.setEnabled(false);
    rightSlider.addChangeListener(this);
    
    zoomIn = new JButton("Zoom in (on selection)");
    zoomIn.addActionListener(this);
    zoomIn.setEnabled(false);
    
    zoomOut = new JButton("Zoom out (show all)");
    zoomOut.addActionListener(this);
    zoomOut.setEnabled(false);
    
    save = new JButton("Save input (PNG)");
    save.addActionListener(this);
    save.setEnabled(false);
    
    clearAll = new JButton("Clear ALL data");
    clearAll.setOpaque(true);
    clearAll.setBackground( Color.RED.darker() );
    clearAll.addActionListener(this);
    clearAll.setEnabled(false);
    
    int yOfClear = gbc.gridy;
    
    //this.add(allCharts);
    gbc.fill = GridBagConstraints.HORIZONTAL;
    gbc.gridy += 2;
    gbc.weighty = 0;
    gbc.weightx = 1;
    gbc.gridwidth = 1;
    gbc.gridx = 0;

    this.add(leftSlider, gbc);
    
    gbc.gridx += 1;
    this.add(rightSlider, gbc);
    
    // now we can add the space between the last plot and the save button
    //this.add( Box.createVerticalStrut(5) );
    
    gbc.fill = GridBagConstraints.NONE;
    gbc.gridy += 1;
    gbc.weighty = 0; gbc.weightx = 0;
    
    gbc.gridx = 0;
    gbc.anchor = GridBagConstraints.EAST;
    this.add(zoomIn, gbc);
    
    gbc.gridx += 1;
    gbc.anchor = GridBagConstraints.WEST;
    this.add(zoomOut, gbc);

    
    gbc.gridwidth = 1;
    gbc.anchor = GridBagConstraints.CENTER;
    gbc.gridx = 7;
    gbc.gridy = yOfClear+1;
    gbc.gridheight = 2;
    gbc.fill = GridBagConstraints.BOTH;

    this.add(save, gbc);

    
    gbc.gridy += 2;
    gbc.gridheight = GridBagConstraints.REMAINDER;

    this.add(clearAll, gbc);

    
    //this.add(buttons);
    
    fc = new JFileChooser();
    
  }
  
  public InputPanel(DataStore ds, int panelsNeeded) {
    this();
    for (int i = 0; i < panelsNeeded; ++i) {
      if ( ds.blockIsSet(i) ) {
        seedFileNames[i].setText( ds.getBlock(i).getName() );
        XYSeries ts = ds.getPlotSeries(i);
        JFreeChart chart = ChartFactory.createXYLineChart(
            ts.getKey().toString(),
            "Time",
            "Counts",
            new XYSeriesCollection(ts),
            PlotOrientation.VERTICAL,
            false, false, false);
        chartPanels[i].setChart(chart);
      }
      if ( ds.responseIsSet(i) ) {
        respFileNames[i].setText( ds.getResponse(i).getName() );
      }
    }
    showDataNeeded(panelsNeeded);
  }
  
  /**
   * Dispatches commands when interface buttons are clicked.
   * When the save button is clicked, dispatches the command to save plots as
   * an image. When the zoom buttons are clicked, scales the plot to only
   * contain the data within a specific range. Resets plots and removes 
   * underlying data when the clear buttons are clicked. Prompts user to
   * load in a file for miniseed and resp data when the corresponding
   * buttons are clicked; because seed files can be large and contian a lot
   * of data to plot, runs seed-loading code backend in a separate thread.
   * 
   * When new data is loaded in, this also fires a change event to any listeners
   */
  @Override
  public void actionPerformed(ActionEvent e) {

    for (int i = 0; i < FILE_COUNT; ++i) {
      JButton clear = clearButton[i];
      JButton seed = seedLoaders[i];
      JButton resp = respLoaders[i];
      
      if( e.getSource() == clear ) {
        instantiateChart(i);
        ds.removeData(i);
        clear.setEnabled(false);
        seedFileNames[i].setText("NO FILE LOADED");
        respFileNames[i].setText("NO FILE LOADED");
        
        // plot all valid range of loaded-in data or else
        // disable clearAll button if there's no other data loaded in
        clearAll.setEnabled( ds.isAnythingSet() );
        
        zooms = new DataStore(ds);
        zooms.trimToCommonTime();
        
        showRegionForGeneration();
        
        fireStateChanged();
      }
      
      if ( e.getSource() == seed ) {
        loadData(i, seed);
      }
      
      if ( e.getSource() == resp ) {
        // don't need a new thread because resp loading is pretty prompt
        
        Set<String> respFilenames = InstrumentResponse.parseInstrumentList();
        
        List<String> names = new ArrayList<String>(respFilenames);
        Collections.sort(names);
        
        String custom = "Load custom response...";
        
        names.add(custom);
        
        JDialog dialog = new JDialog();
        Object result = JOptionPane.showInputDialog(
            dialog,
            "Select a response to load:",
            "RESP File Selection",
            JOptionPane.PLAIN_MESSAGE,
            null, names.toArray(),
            names.get( names.size() - 1 ) );
        
        String resultStr = (String) result;
        
        if (resultStr == null) {
          return;
        }
        
        if ( respFilenames.contains(resultStr) ) {
          // final used here in the event of thread weirdness
          final String fname = resultStr;
          try {
            InstrumentResponse ir = 
                InstrumentResponse.loadEmbeddedResponse(fname);
            ds.setResponse(i, ir);
            zooms.setResponse(i, ir);
            respFileNames[i].setText( ir.getName() );
            clear.setEnabled(true);
            clearAll.setEnabled(true);
            
            fireStateChanged();
          } catch (IOException e1) {
            e1.printStackTrace();
          }
        } else {
          fc.setCurrentDirectory( new File(respDirectory) );
          fc.resetChoosableFileFilters();
          fc.setDialogTitle("Load response file...");
          int returnVal = fc.showOpenDialog(resp);
          if (returnVal == JFileChooser.APPROVE_OPTION) {
            File file = fc.getSelectedFile();
            respDirectory = file.getParent();
            ds.setResponse(i, file.getAbsolutePath() );
            zooms.setResponse(i, file.getAbsolutePath() );
            respFileNames[i].setText( file.getName() );
            clear.setEnabled(true);
            clearAll.setEnabled(true);
            
            fireStateChanged();
          }
        }
        

        return;
      }
      
    }
    
    if ( e.getSource() == clearAll) {
      clearAllData();
      
      fireStateChanged();
      
      return;
    }
    
    
    if ( e.getSource() == save ) {
      String ext = ".png";
      fc = new JFileChooser();
      fc.setCurrentDirectory( new File(saveDirectory) );
      fc.addChoosableFileFilter(
          new FileNameExtensionFilter("PNG image (.png)",ext) );
      fc.setFileFilter(fc.getChoosableFileFilters()[1]);
      int returnVal = fc.showSaveDialog(save);
      if (returnVal == JFileChooser.APPROVE_OPTION) {
        File selFile = fc.getSelectedFile();
        saveDirectory = selFile.getParent();
        if( !selFile.getName().endsWith( ext.toLowerCase() ) ) {
          selFile = new File( selFile.toString() + ext);
        }
        try {
          int height = IMAGE_HEIGHT * activePlots;
          BufferedImage bi = getAsImage( 640, height, activePlots );
          
          ImageIO.write(bi,"png",selFile);
        } catch (IOException e1) {
          e1.printStackTrace();
        }
        
      }
      return;
    }
    
    if ( e.getSource() == zoomIn ) {
      
      showRegionForGeneration();
      return;
    }
    
    if ( e.getSource() == zoomOut ) {
      
      // restore original loaded datastore
      zooms = new DataStore(ds);
      zooms.trimToCommonTime();
      for (int i = 0; i < FILE_COUNT; ++i) {
        if ( !zooms.blockIsSet(i) ) {
          continue;
        }
        resetPlotZoom(i);
      }
      
      leftSlider.setValue(0); rightSlider.setValue(SLIDER_MAX);
      setVerticalBars();
      zoomOut.setEnabled(false);
      return;
    }
  }
  
  /**
   * Used to add objects to the list that will be informed when data is loaded
   * @param listener
   */
  public void addChangeListener(ChangeListener listener) {
    listenerList.add(ChangeListener.class, listener);
  }
  
  /**
   * Resets the data and blanks out all charts
   */
  public void clearAllData() {
    ds = new DataStore();
    zooms = new DataStore();
    
    zoomIn.setEnabled(false);
    zoomOut.setEnabled(false);
    
    leftSlider.setEnabled(false);
    rightSlider.setEnabled(false);
    
    clearAll.setEnabled(false);
    save.setEnabled(false);
    
    for (JTextComponent fn : seedFileNames) {
      fn.setText("NO FILE LOADED");
    }
    for (JTextComponent fn : respFileNames) {
      fn.setText("NO FILE LOADED");
    }
    
    for (int i = 0; i < chartPanels.length; ++i) {
      clearButton[i].setEnabled(false);
      instantiateChart(i);
    }
    
  }
  
  /**
   * Informs listening objects that the state of the inputs has changed
   * This is done when new seed or resp data has been loaded in, mainly
   * to tell whether enough data exists to run one of the experiments
   */
  protected void fireStateChanged() {
    ChangeListener[] lsners = listenerList.getListeners(ChangeListener.class);
    if (lsners != null && lsners.length > 0) {
      ChangeEvent evt = new ChangeEvent(this);
      for (ChangeListener lsnr : lsners) {
        lsnr.stateChanged(evt);
      }
    }
  }
  
  /**
   * Return this panel's charts as a single buffered image
   * @return Buffered image of the plots, writeable to file
   */
  public BufferedImage getAsImage(int plotsToShow) {
    // if all 3 plots are set, height of panel is height of image
    int height = getImageHeight(plotsToShow);
    // otherwise, we only use the height of images actually set
    return getAsImage( IMAGE_WIDTH, height, plotsToShow );
  }
  
  /**
   * Return this panel's charts as a single buffered image 
   * with specified dimensions
   * @param width Width of returned image
   * @param height Height of returned image
   * @return Buffered image of the plots, writeable to file
   */
  public BufferedImage getAsImage(int width, int height, int plotsToShow) {
    
    if (plotsToShow <= 0) {
      // should never be called like this but just in case
      // return an empty image
      return new BufferedImage(0, 0, BufferedImage.TYPE_INT_RGB);
    }
    
    // int shownHeight = allCharts.getHeight();
    
    // width = Math.min( width, chartPanels[0].getWidth() );
    // height = Math.max( height, shownHeight );
    
    int loaded = plotsToShow;
    // cheap way to make sure height is a multiple of the chart count
    height = (height*loaded)/loaded;
    int chartHeight = height/loaded;
    
    JFreeChart[] chartsToPrint = new JFreeChart[plotsToShow];
    for (int i = 0; i < plotsToShow; ++i) {
      chartsToPrint[i] = chartPanels[i].getChart();
    }
    
    return ReportingUtils.chartsToImage(width, chartHeight, chartsToPrint);

  }
  
  
  /**
   * Returns the selected region of underlying DataStore, to be fed 
   * into experiments for processing (the results of which will be plotted)
   * When this function is called, the graphs zoom to the currently active
   * range to display selected by the sliders, which is also the range
   * passed into an experiment
   * @return A DataStore object (contains arrays of DataBlocks & Responses)
   */
  public DataStore getData() {
    
    showRegionForGeneration();
    if ( zooms.numberOfBlocksSet() > 1 ) {
      zooms.matchIntervals(activePlots);
      zooms.trimToCommonTime(activePlots);
    }

    return zooms;
  }

  /**
   * Gets the height of resulting image of plots given default parameters,
   * so that it only needs to fit the plots that have data in them 
   * @return height of image to output
   */
  public int getImageHeight(int plotsToShow) {
    return IMAGE_HEIGHT * plotsToShow;
  }
  
  /**
   * Returns a default image width for writing plots to file
   * @return width of image to output
   */
  public int getImageWidth()  {
    return allCharts.getWidth();
  }
  
  /**
   * Instantiates the underlying chart of a chartpanel with default data
   * @param idx Index of the chartpanel to instantiate
   */
  private void instantiateChart(int idx) {
    JFreeChart chart = ChartFactory.createXYLineChart(
        "SEED input " + (idx + 1),
        "Time",
        "Counts",
        new XYSeriesCollection(),
        PlotOrientation.VERTICAL,
        false, false, false);
    
    if (chartPanels[idx] == null) {
      chartPanels[idx] = new ChartPanel(chart);
    } else {
      chartPanels[idx].setChart(chart);

    }
    chartPanels[idx].setMouseZoomable(true);
  }
  
  /**
   * Used to remove an object from the list of those informed when
   * data is loaded in or cleared out
   * @param listener
   */
  public void removeChangeListener(ChangeListener listener) {
    listenerList.remove(ChangeListener.class, listener);
  }
  
  /**
   * Does the work to reset the zoom of a chart when the zoom button is hit
   * @param idx Index of appropriate chart/panel
   */
  private void resetPlotZoom(int idx) {
    XYPlot xyp = chartPanels[idx].getChart().getXYPlot();
    XYSeriesCollection xys = new XYSeriesCollection();
    xys.addSeries( zooms.getBlock(idx).toXYSeries() );
    xyp.setDataset( xys );
    xyp.getRenderer().setSeriesPaint(0,
        defaultColor[idx % defaultColor.length]);
    if ( xyp.getSeriesCount() > 1 ) {
      throw new RuntimeException("TOO MUCH DATA");
    }
    chartPanels[idx].repaint();
  }
  
  /**
   * Parent function to load in a response file to this panel's
   * DataStore object
   * @param idx Index of the file/chart this response corresponds to
   * @param filepath Full address of the file to be loaded in
   */
  public void setResponse(int idx, String filepath) {
    
    ds.setResponse(idx, filepath);
    zooms.setResponse(idx, filepath);
  }
  
  /**
   * Displays the range set by the sliders using
   * vertical bars at the min and max values
   */
  public void setVerticalBars() {
    
    // zooms.trimToCommonTime();
    
    for (int i = 0; i < FILE_COUNT; ++i) {
      if ( !zooms.blockIsSet(i) ) {
        continue;
      }
      
      int leftValue = leftSlider.getValue();
      int rightValue = rightSlider.getValue();
      
      XYPlot xyp = (XYPlot) chartPanels[i].getChart().getPlot();
      xyp.clearDomainMarkers();
      
      DataBlock db = zooms.getBlock(i);
      
      long startMarkerLocation = getMarkerLocation(db, leftValue);
      long endMarkerLocation = getMarkerLocation(db, rightValue);
      
      // divide by 1000 here to get time value in ms
      Marker startMarker = new ValueMarker(startMarkerLocation/1000);
      startMarker.setStroke( new BasicStroke( (float) 1.5 ) );
      Marker endMarker = new ValueMarker(endMarkerLocation/1000);
      endMarker.setStroke( new BasicStroke( (float) 1.5 ) );
      
      xyp.addDomainMarker(startMarker);
      xyp.addDomainMarker(endMarker);
      
      chartPanels[i].repaint();
    }
    
    
  }

  /**
   * Zooms in on the current range of data, which will be passed into
   * backend functions for experiment calculations
   */
  public void showRegionForGeneration() {
    
    if ( zooms.numberOfBlocksSet() < 1 ) {
      return;
    }
    
    // get (any) loaded data block to map slider to domain boundary
    // all data should have the same range
    DataBlock db = zooms.getXthLoadedBlock(1);

    if ( leftSlider.getValue() != 0 || rightSlider.getValue() != SLIDER_MAX ) {
      long start = getMarkerLocation(db, leftSlider.getValue() );
      long end = getMarkerLocation(db, rightSlider.getValue() );
      zooms = new DataStore(ds, start, end, activePlots);
      leftSlider.setValue(0); rightSlider.setValue(SLIDER_MAX);
      zoomOut.setEnabled(true);
    }
    
    for (int i = 0; i < activePlots; ++i) {
      if ( !zooms.blockIsSet(i) ) {
        continue;
      }
      resetPlotZoom(i);
    }

    setVerticalBars();
    
  }

  @Override
  /**
   * Handles changes in value by the sliders below the charts
   */
  public void stateChanged(ChangeEvent e) {
    
    int leftSliderValue = leftSlider.getValue();
    int rightSliderValue = rightSlider.getValue();
    
    // probably can refactor this
    // the conditionals are effectively the same
    
    if ( e.getSource() == leftSlider ) {
      if (leftSliderValue > rightSliderValue || 
          leftSliderValue + MARGIN > rightSliderValue) {
        leftSliderValue = rightSliderValue - MARGIN;
        if (leftSliderValue < 0) {
          leftSliderValue = 0;
          rightSliderValue = MARGIN;
        }
      }
    } else if ( e.getSource() == rightSlider ) {
      if (rightSliderValue < leftSliderValue ||
          rightSliderValue - MARGIN < leftSliderValue) {
        rightSliderValue = leftSliderValue + MARGIN;
        if (rightSliderValue > SLIDER_MAX) {
          rightSliderValue = SLIDER_MAX;
          leftSliderValue = SLIDER_MAX - MARGIN;
        }
      }
    }
    
    leftSlider.setValue(leftSliderValue);
    rightSlider.setValue(rightSliderValue);
    
    setVerticalBars();
    
  }
  
  /**
   * Load in data for a specified SEED file, to be run in a specific thread.
   * Because loading can be a slow operation, this runs in a background thread
   * @param idx Index into datastore/plots this data should be loaded
   * @param seed The JButton to passed into the fileloader
   */
  private void loadData(final int idx, final JButton seed) {

    fc.setCurrentDirectory( new File(seedDirectory) );
    fc.resetChoosableFileFilters();
    fc.setDialogTitle("Load SEED file...");
    int returnVal = fc.showOpenDialog(seed);
    if (returnVal == JFileChooser.APPROVE_OPTION) {
      final File file = fc.getSelectedFile();
      seedDirectory = file.getParent();
      String oldName = seedFileNames[idx].getText();

      seedFileNames[idx].setText("LOADING: " + file.getName());
      final String filePath = file.getAbsolutePath();
      String filterName = "";
      try {
        Set<String> nameSet = TimeSeriesUtils.getMplexNameSet(filePath);

        if (nameSet.size() > 1) {
          // more than one series in the file? prompt user for it
          String[] names = nameSet.toArray(new String[0]);
          Arrays.sort(names);
          JDialog dialog = new JDialog();
          Object result = JOptionPane.showInputDialog(
              dialog,
              "Select the subseries to load:",
              "Multiplexed File Selection",
              JOptionPane.PLAIN_MESSAGE,
              null, names,
              names[0]);
          if (result instanceof String) {
            filterName = (String) result;
          } else {
            // if the user cnacelled selecting a subseries
            seedFileNames[idx].setText(oldName);
            return;
          }
        } else {
          // just get the first one; it's the only one in the list
          filterName = new ArrayList<String>(nameSet).get(0);
        }

      } catch (FileNotFoundException e1) {
        e1.printStackTrace();
        return;
      }

      final String immutableFilter = filterName;

      // create swingworker to load large files in the background
      SwingWorker<Integer, Void> worker = new SwingWorker<Integer, Void>(){

        JFreeChart chart;
        boolean caughtException = false;

        @Override
        public Integer doInBackground() {

          try {
            ds.setData(idx, filePath, immutableFilter);
          } catch (RuntimeException e) {
            e.printStackTrace();
            caughtException = true;
            return 1;
          }

          zooms = new DataStore(ds);
          zooms.matchIntervals(activePlots);
          zooms.trimToCommonTime(activePlots);

          XYSeries ts = zooms.getPlotSeries(idx);
          chart = ChartFactory.createXYLineChart(
              ts.getKey().toString(),
              "Time",
              "Counts",
              new XYSeriesCollection(ts),
              PlotOrientation.VERTICAL,
              false, false, false);

          XYPlot xyp = (XYPlot) chart.getPlot();
          DateAxis da = new DateAxis();

          SDF.setTimeZone( TimeZone.getTimeZone("UTC") );
          da.setLabel("UTC Time (Year.Day.Hour:Minute)");
          Font bold = da.getLabelFont();
          bold = bold.deriveFont(Font.BOLD);
          da.setLabelFont(bold);
          da.setDateFormatOverride(SDF);
          xyp.setDomainAxis(da);
          int colorIdx = idx % defaultColor.length;
          xyp.getRenderer().setSeriesPaint(0, defaultColor[colorIdx]);
          
          return 0;
          // setData(idx, filePath, immutableFilter);
          // return 0;
        }

        public void done() {

          if (caughtException) {
            instantiateChart(idx);
            XYPlot xyp = (XYPlot) chartPanels[idx].getChart().getPlot();
            TextTitle result = new TextTitle();
            String errMsg = "COULD NOT LOAD IN " + file.getName();
            errMsg += "\nTIME RANGE DOES NOT INTERSECT";
            result.setText(errMsg);
            result.setBackgroundPaint(Color.red);
            result.setPaint(Color.white);
            XYTitleAnnotation xyt = new XYTitleAnnotation(0.5, 0.5, result,
                RectangleAnchor.CENTER);
            xyp.clearAnnotations();
            xyp.addAnnotation(xyt);
            seedFileNames[idx].setText("NO FILE LOADED");
            clearButton[idx].setEnabled(true);
            return;
          }
          
          // seedFileNames[idx].setText("PLOTTING: " + file.getName());

          chartPanels[idx].setChart(chart);
          chartPanels[idx].setMouseZoomable(true);

          clearButton[idx].setEnabled(true);

          for (int i = 0; i < FILE_COUNT; ++i) {
            if ( !zooms.blockIsSet(i) ) {
              continue;
            }
            resetPlotZoom(i);
          }
          
          leftSlider.setValue(0); rightSlider.setValue(SLIDER_MAX);
          setVerticalBars();

          zoomOut.setEnabled(false);
          zoomIn.setEnabled(true);
          leftSlider.setEnabled(true);
          rightSlider.setEnabled(true);
          save.setEnabled(true);
          clearAll.setEnabled(true);

          leftSlider.setValue(0);
          rightSlider.setValue(SLIDER_MAX);
          setVerticalBars();

          seedFileNames[idx].setText( 
              file.getName() + ": " + immutableFilter);

          fireStateChanged();
        }
      };

      worker.execute();
      return;
    }
  }

  /**
   * Show the number of panels needed to load in data for a specific experiment
   * @param panelsNeeded Number of panels to show
   */
  public void showDataNeeded(int panelsNeeded) {
    
    VScrollPanel cont = new VScrollPanel();
    cont.setLayout( new GridBagLayout() );
    GridBagConstraints contConstraints = new GridBagConstraints();
    contConstraints.weightx = 1.0;
    contConstraints.weighty = 1.0;
    contConstraints.gridy = 0;
    contConstraints.anchor = GridBagConstraints.CENTER;
    contConstraints.fill = GridBagConstraints.BOTH;
    
    activePlots = panelsNeeded;
    
    zooms = new DataStore(ds, activePlots);
    zooms.trimToCommonTime(activePlots);
    int i;
    for (i = 0; i < activePlots; ++i) {
      if ( zooms.blockIsSet(i) ){
        resetPlotZoom(i);
      }
      
      cont.add(chartSubpanels[i], contConstraints);
      contConstraints.gridy += 1;
      // gbc.gridy += 1;
    }
    
    leftSlider.setValue(0); rightSlider.setValue(SLIDER_MAX);
    setVerticalBars();
    
    zoomIn.setEnabled( zooms.numberOfBlocksSet() > 0 );
    
    zoomOut.setEnabled(false);
    
    // using this test means the panel doesn't try to scroll when it's
    // only got a few inputs to deal with, when stuff is still pretty readable
    cont.setScrollableTracksViewportHeight(activePlots <= MAX_UNSCROLLED);
    
    inputScrollPane.getViewport().setView(cont);
    inputScrollPane.setPreferredSize( cont.getPreferredSize() );
  }
  
  /**
   * Used to get labels for each plot to idenitify what data they need to
   * contain in order for an experiment to have enough data to run
   * @param channels List of strings to be used as panel title
   */
  public void setChannelTypes(String[] channels) {
    
    int len = Math.min(channels.length, channelType.length);
    
    for (int i = 0; i < len; ++i) {
      channelType[i].setText(channels[i]);
      channelType[i].setHorizontalAlignment(SwingConstants.CENTER);
    }
  }
  
}

