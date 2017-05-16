package asl.sensor;

import java.awt.Component;
import java.awt.Dimension;
import java.awt.Graphics2D;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;

import javax.imageio.ImageIO;
import javax.swing.Box;
import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JSplitPane;
import javax.swing.JTabbedPane;
import javax.swing.SwingUtilities;
import javax.swing.SwingWorker;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.filechooser.FileNameExtensionFilter;

import org.apache.log4j.BasicConfigurator;
import org.apache.pdfbox.pdmodel.PDDocument;
import org.jfree.chart.JFreeChart;

import asl.sensor.experiment.ExperimentEnum;
import asl.sensor.gui.ExperimentPanel;
import asl.sensor.gui.ExperimentPanelFactory;
import asl.sensor.gui.InputPanel;
import asl.sensor.input.DataStore;
import asl.sensor.utils.ReportingUtils;

/**
 * Main window of the sensor test program and the program's launcher
 * Mainly used for handling the input (InputPanel) 
 * and output (ExperimentPanel)
 * GUI frames and making sure they fit togther and cooperate.
 * @author akearns
 *
 */
public class SensorSuite extends JPanel 
                         implements ActionListener, ChangeListener {

  /**
   * 
   */
  private static final long serialVersionUID = 2866426897343097822L;
 

  /**
   * Loads the main window for the program on launch
   */
  private static void createAndShowGUI() {
    JFrame frame = new JFrame("Sensor Tests");
    frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

    frame.add( new SensorSuite() );

    frame.pack();
    frame.setVisible(true);
  }
  /**
   * Create space with the given dimensions as a buffer between plots
   * when writing the plots to an image file
   * @param width The width of the spacer
   * @param height The height of the spacer
   * @return A blank BufferedImage that can be concatenated with other plots
   */
  public static BufferedImage getSpace(int width, int height) {
    BufferedImage space = new BufferedImage(
        width, 
        height, 
        BufferedImage.TYPE_INT_RGB);
    Graphics2D tmp = space.createGraphics();
    JPanel margin = new JPanel();
    margin.add( Box.createRigidArea( new Dimension(width,height) ) );
    margin.printAll(tmp);
    tmp.dispose();

    return space;
  }
  /**
   * Starts the program -- instantiate the top-level GUI
   * @param args (Any parameters fed in on command line are currently ignored)
   */
  public static void main(String[] args) {
    //Schedule a job for the event dispatch thread:
    //creating and showing this application's GUI.
    
    SwingUtilities.invokeLater(new Runnable() {
      public void run() {
        BasicConfigurator.configure();
        //Turn off metal's use of bold fonts
        // UIManager.put("swing.boldMetal", Boolean.FALSE); 
        createAndShowGUI();
      }
    });

  }

  /**
   * Plots the data from an output panel and its associated input
   * @param file Filename to write to
   * @param ep Experiment panel with data to be plotted
   * @param ip Input panel holding data associated with the experiment
   * @throws IOException If the file cannot be written to
   */
  public static void plotsToPDF(File file, ExperimentPanel ep, InputPanel ip)
      throws IOException{

    int inPlotCount = ep.plotsToShow();
    String[] responses = ip.getResponseStrings( ep.getResponseIndices() );
    // BufferedImage toFile = getCompiledImage();

    // START OF UNIQUE CODE FOR PDF CREATION HERE
    PDDocument pdf = new PDDocument();
    ep.savePDFResults( pdf );
    
    if (inPlotCount > 0) {

      int inHeight = ip.getImageHeight(inPlotCount) * 2;
      int width = 1280; // TODO: set as global static variable somewhere?

      BufferedImage[] toFile = 
          ip.getAsMultipleImages(width, inHeight, inPlotCount);
      
      ReportingUtils.imageListToPDFPages(pdf, toFile);

      // TODO: response string generation goes here
      ReportingUtils.textListToPDFPages(pdf, responses);
    }

    pdf.save(file);
    pdf.close();
  }
  
  private JFileChooser fc; // loads in files based on parameter

  private InputPanel inputPlots;

  private JTabbedPane tabbedPane; // holds set of experiment panels

  private JButton generate, savePDF; // run all calculations

  // used to store current directory locations
  private String saveDirectory = System.getProperty("user.home");

  private SwingWorker<Integer, Void> worker;
  
  /**
   * Creates the main window of the program when called
   * (Three main panels: the top panel for displaying the results
   * of sensor tests; the lower panel for displaying plots of raw data from
   * miniSEED files; the side panel for most file-IO operations
   */
  public SensorSuite() {

    super();
    
    // set up experiment panes in a tabbed pane
    tabbedPane = new JTabbedPane();
    
    for ( ExperimentEnum exp : ExperimentEnum.values() ) {
      JPanel tab = ExperimentPanelFactory.createPanel(exp);
      tabbedPane.addTab( exp.getName(), tab);
    }
    
    Dimension d = tabbedPane.getPreferredSize();
    d.setSize( d.getWidth() * 1.5, d.getHeight() );
    tabbedPane.setPreferredSize(d);
    
    tabbedPane.addChangeListener(this);
    
    inputPlots = new InputPanel();
    inputPlots.addChangeListener(this);
    
    // experiments on left, input on the right; split to allow resizing
    JSplitPane mainSplit = 
        new JSplitPane(JSplitPane.HORIZONTAL_SPLIT, true);
        // boolean allows redrawing immediately on resize
    mainSplit.setLeftComponent(tabbedPane);
    mainSplit.setRightComponent(inputPlots);
    // set the left-pane to resize more when window is horizontally stretched
    mainSplit.setResizeWeight(.5);
    mainSplit.setOneTouchExpandable(true);
    
    // we want to make sure the split pane fills the window
    this.setLayout( new GridBagLayout() );
    GridBagConstraints c = new GridBagConstraints();
    c.gridx = 0; c.gridy = 0;
    c.weightx = 1; c.weighty = 1;
    c.gridwidth = 2;
    c.fill = GridBagConstraints.BOTH;
    
    this.add(mainSplit, c);
    
    c.fill = GridBagConstraints.BOTH;
    c.gridx = 0; c.gridy = 1;
    c.weightx = 1.0; c.weighty = 0.0;
    c.gridwidth = 1;
    
    // now add the buttons
    savePDF = new JButton("Save input and output plots (PDF)");
    savePDF.setEnabled(false);
    savePDF.addActionListener(this);
    this.add(savePDF, c);
    c.gridx += 1;

    generate = new JButton("Generate test result");
    generate.setEnabled(false);
    generate.addActionListener(this);
    d = generate.getPreferredSize();
    d.setSize( d.getWidth(), d.getHeight() * 2 );
    generate.setPreferredSize(d);
    this.add(generate, c);
    
    fc = new JFileChooser();
    
    ExperimentPanel ep = (ExperimentPanel) tabbedPane.getSelectedComponent();
    inputPlots.showDataNeeded( ep.panelsNeeded() );

    inputPlots.setChannelTypes( ep.getChannelTypes() );
    
  }
  
  /**
   * Handles actions when the buttons are clicked -- either the 'save PDF'
   * button, which compiles the input and output plots into a single PDF, or
   * the 'generate result' button.
   * Because generating results of an experiment can be slow, the operation
   * is set to run in a separate thread.
   */
  @Override
  public void actionPerformed(ActionEvent e) {


    if ( e.getSource() == generate ) {
      
      if (worker != null) {
        savePDF.setEnabled(false);
        // clear out any old data in the chart
        // since we only have one worker thread for experiment calculations
        // we clear all charts -- if we run a new experiment while another one
        // is calculating from another panel, we should make it clear that
        // other panel was cancelled, and thus clear the chart / unset data
        worker.cancel(true); // doesn't matter if experiment already completed
        for (Component component : tabbedPane.getComponents()) {
          ExperimentPanel ep = (ExperimentPanel) component;
          if ( ep.hasRun() ) {
            ep.clearChart();
          }
        }
      }
      
      worker = new SwingWorker<Integer, Void>() {
        
        @Override
        public Integer doInBackground() {
          resetTabPlots();
          return 0;
        }
        
      };
      worker.execute();
      return;
    } else if ( e.getSource() == savePDF ) {

      String ext = ".pdf";
      fc.setCurrentDirectory( new File(saveDirectory) );
      fc.addChoosableFileFilter(
          new FileNameExtensionFilter("PDF file (.pdf)",ext) );
      fc.setFileFilter(fc.getChoosableFileFilters()[1]);
      fc.setDialogTitle("Save PDF report...");
      ExperimentPanel ep = (ExperimentPanel) tabbedPane.getSelectedComponent();
      String defaultName = ep.getPDFFilename();

      String text = ep.getAllTextData();
      JFreeChart[] charts = ep.getCharts();
      
      fc.setSelectedFile( new File(defaultName) );
      int returnVal = fc.showSaveDialog(savePDF);
      
      if (returnVal == JFileChooser.APPROVE_OPTION) {
        
        File selFile = fc.getSelectedFile();
        saveDirectory = selFile.getParent();
        if( !selFile.getName().endsWith( ext.toLowerCase() ) ) {
          selFile = new File( selFile.getName() + ext);
        }
        
        // start in the folder the pdf is saved, add data into a new
        // subfolder for calibration data
        StringBuilder folderName = new StringBuilder( selFile.getParent() );
        folderName.append("/test_results/");
        folderName.append( selFile.getName().replace(".pdf","") );
        
        File folder = new File( folderName.toString() );
        if ( !folder.exists() ) {
          System.out.println("Writing directory " + folderName);
          folder.mkdirs();
        }
        
        String textName = folderName + "/outputData.txt";
        
        try {
          PrintWriter out = new PrintWriter(textName);
          out.println(text);
          out.close();
        } catch (FileNotFoundException e2) {
          System.out.println("Can't write the text");
          e2.printStackTrace();
        }
        
        for (int i = 0; i < charts.length; ++i) {
          JFreeChart chart = charts[i];
          String plotName = folderName + "/chart" + (i + 1) + ".png";
          BufferedImage chartImage = 
              ReportingUtils.chartsToImage(1280, 960, chart);
          File plotPNG = new File(plotName);
          try {
            ImageIO.write(chartImage, "png", plotPNG);
          } catch (IOException e1) {
            e1.printStackTrace();
          }
        }
        
        try {
          plotsToPDF(selFile);
        } catch (IOException e1) {
          e1.printStackTrace();
        }
      }
      return;
    }


  }
  
  /**
   * Produces a buffered image of all active charts
   * @return BufferedImage that can be written to file
   */
  private BufferedImage getCompiledImage() {
    
    ExperimentPanel ep = (ExperimentPanel) tabbedPane.getSelectedComponent();
    int inPlotCount = ep.plotsToShow();
    
    int width = 1280;
    BufferedImage outPlot = ep.getAsImage(width, 960);
    
    width = outPlot.getWidth();
    
    int inHeight = inputPlots.getImageHeight(inPlotCount) * 2;

    int height = outPlot.getHeight();

    BufferedImage inPlot = null; // unfortunate, but can't make an empty BI
    if (inPlotCount > 0) {
      inPlot = inputPlots.getAsImage(width, inHeight, inPlotCount);
      height = inPlot.getHeight() + outPlot.getHeight();
    }

    // int width = Math.max( inPlot.getWidth(), outPlot.getWidth() );
    // 5px tall buffer used to separate result plot from inputs
    // BufferedImage space = getSpace(width, 0);


    // System.out.println(space.getHeight());

    BufferedImage returnedImage = new BufferedImage(width, height, 
        BufferedImage.TYPE_INT_ARGB);

    Graphics2D combined = returnedImage.createGraphics();
    combined.drawImage(outPlot, null, 0, 0);
    if (null != inPlot) {
      combined.drawImage( inPlot, null, 0, 
          outPlot.getHeight() );
    }

    combined.dispose();
    
    return returnedImage;
  }
  
  /**
   * Output the currently displayed plots as a PDF file.
   * @param file Filename to write to
   * @throws IOException If the file cannot be written
   */
  public void plotsToPDF(File file) throws IOException {

    ExperimentPanel ep = (ExperimentPanel) tabbedPane.getSelectedComponent();
    InputPanel ip = inputPlots;

    plotsToPDF(file, ep, ip);
  }


  /**
   * Handles function to create a PNG image with all currently-displayed plots
   * (active experiment and read-in time series data)
   * @param file File (PNG) that image will be saved to
   * @throws IOException
   * @deprecated Use PDF output (plotstoPDF) instead
   */
  @Deprecated
  public void plotsToPNG(File file) throws IOException {

    // just write the bufferedimage to file
    ImageIO.write( getCompiledImage(), "png", file );

  }

  /**
   * Resets plot data and gets the inputted data to send to experiments
   */
  private void resetTabPlots() {
    
    ExperimentPanel ep = (ExperimentPanel) tabbedPane.getSelectedComponent();
    
    // pass the inputted data to the panels that handle them
    DataStore ds = inputPlots.getData();
    // update the input plots to show the active region being calculated
    inputPlots.showRegionForGeneration();
    
    // now, update the data

    ep.updateData(ds);
    
    savePDF.setEnabled( ep.hasRun() );
  }

  /**
   * Checks when input panel gets new data or active experiment changes
   * to determine whether or not the experiment can be run yet
   */
  @Override
  public void stateChanged(ChangeEvent e) {
    if ( e.getSource() == inputPlots ){
      ExperimentPanel ep = (ExperimentPanel) tabbedPane.getSelectedComponent();
      DataStore ds = inputPlots.getData();
      boolean canGenerate = ep.hasEnoughData(ds);
      generate.setEnabled(canGenerate);
    } else if ( e.getSource() == tabbedPane ) {
      
      ExperimentPanel ep = (ExperimentPanel) tabbedPane.getSelectedComponent();
      
      inputPlots.setChannelTypes( ep.getChannelTypes() );
      
      inputPlots.showDataNeeded( ep.panelsNeeded() );
      DataStore ds = inputPlots.getData();
      boolean canGenerate = ep.hasEnoughData(ds);
      boolean isSet = ep.hasRun();
      generate.setEnabled(canGenerate);
      savePDF.setEnabled(canGenerate && isSet);
    }
  }


}
