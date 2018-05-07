package asl.sensor.gui;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.time.Instant;
import java.time.ZoneOffset;
import java.time.format.DateTimeFormatter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import javax.imageio.ImageIO;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
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
import org.apache.commons.math3.util.Pair;
import org.jfree.chart.ChartColor;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.annotations.XYTitleAnnotation;
import org.jfree.chart.axis.DateAxis;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.IntervalMarker;
import org.jfree.chart.plot.Marker;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.ValueMarker;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.title.TextTitle;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.RectangleAnchor;
import asl.sensor.input.DataBlock;
import asl.sensor.input.DataStore;
import asl.sensor.input.InstrumentResponse;
import asl.sensor.utils.ReportingUtils;
import asl.sensor.utils.TimeSeriesUtils;
import edu.iris.dmc.seedcodec.CodecException;
import edu.sc.seis.seisFile.mseed.SeedFormatException;


/**
 * Panel used to hold the plots for the files taken in as input
 * Handles the UI for loading SEED and RESP files, plotting their data
 * and selecting regions of that data, as well as outputting the plots to PNG
 * The data exists in two forms: one is an otherwise-hidden DataStore object
 * that is used to hold in the loaded data directly, and the other one is also
 * a DataStore object but has its data passed to the panel's input plots and
 * to experiment calculations.
 * Some important functions this object has is to make sure that input data is all trimmed
 * to the same region for use in experiment backends, and display time regions of interest
 *
 * @author akearns - KBRWyle
 */
public class InputPanel
    extends JPanel
    implements ActionListener, ChangeListener {

  private static final long serialVersionUID = -7302813951637543526L;

  /**
   * Default height of image produced by the save-as-image function
   * (each chart is 240 pixels tall)
   */
  private static final int IMAGE_HEIGHT = 240;
  private static final int IMAGE_WIDTH = 640;

  private static final int MAX_UNSCROLLED = 3;

  private static final int PLOTS_PER_PAGE = 3;

  private static final int FILE_COUNT = DataStore.FILE_COUNT;

  /**
   * Minimum space of the two sliders.
   */
  private static final int MARGIN = 10;
  public static final int SLIDER_MAX = 10000;

  /**
   * Gets the value of start or end time from slider value and DataBlock
   *
   * @param dataBlock DataBlock corresponding to one of the plots
   * @param sliderValue Value of starting or ending time slider [0-SLIDER_MAX]
   * @return Long that represents start or end time matching slider's value
   */
  public static long getMarkerLocation(DataBlock dataBlock, int sliderValue) {
    long start = dataBlock.getStartTime();
    long len = (dataBlock.getInterval()) * dataBlock.size();
    return start + (sliderValue * len) / SLIDER_MAX; // start + time offset
  }

  private static int getSliderValue(DataBlock dataBlock, long timeStamp) {
    long start = dataBlock.getStartTime();
    long length = dataBlock.getInterval() * dataBlock.size();
    return (int) ((SLIDER_MAX * (timeStamp - start)) / length);
  }

  private int activePlots = FILE_COUNT; // how much data is being displayed

  private DataStore dataStore; // holds data to be plotted in each chartpanel
  private final ChartPanel[] chartPanels; // show the data for a given input
  private final Color[] defaultColor = {
      ChartColor.LIGHT_RED,
      ChartColor.LIGHT_BLUE,
      ChartColor.LIGHT_GREEN}; // control ordering of colors
  private final JButton save; // save every input of note in a png plot
  private final JButton zoomIn; // select a given window
  private final JButton zoomOut; // revert to full data region
  private final JButton clearAll; // remove all data
  private JFileChooser fileChooser;
  private final JSlider leftSlider;
  private final JSlider rightSlider;
  private final JScrollPane inputScrollPane;

  private final EditableDateDisplayPanel startDate;
  private final EditableDateDisplayPanel endDate;

  private final FileOperationJButton[] seedLoaders;
  private final FileOperationJButton[] seedAppenders;
  private final JTextComponent[] seedFileNames;
  private final JButton[] respLoaders;
  private final JTextComponent[] respFileNames;
  private final JButton[] clearButton;

  private final JLabel[] channelType;

  private final JPanel[] chartSubpanels;

  // used to store current directory locations
  private String seedDirectory = "data";
  private String respDirectory = "responses";

  private int lastRespIndex;

  private String saveDirectory = System.getProperty("user.home");

  /**
   * Creates a new data panel -- instantiates each chart, to be populated with
   * data when a file is loaded in. Also creates a save button for writing all
   * the inputted data plots into a single PNG file.
   */
  public InputPanel() {

    chartPanels = new ChartPanel[FILE_COUNT];
    seedLoaders = new FileOperationJButton[FILE_COUNT];
    seedAppenders = new FileOperationJButton[FILE_COUNT];
    seedFileNames = new JTextComponent[FILE_COUNT];
    respFileNames = new JTextComponent[FILE_COUNT];
    respLoaders = new JButton[FILE_COUNT];
    clearButton = new JButton[FILE_COUNT];
    channelType = new JLabel[FILE_COUNT];
    chartSubpanels = new JPanel[FILE_COUNT];

    this.setLayout(new GridBagLayout());
    GridBagConstraints constraints = new GridBagConstraints();

    dataStore = new DataStore();

    constraints.weightx = 1.0;
    constraints.weighty = 1.0;
    constraints.gridy = 0;
    constraints.gridwidth = 8;
    constraints.anchor = GridBagConstraints.CENTER;
    constraints.fill = GridBagConstraints.BOTH;

    inputScrollPane = new JScrollPane();
    inputScrollPane.setVerticalScrollBarPolicy(
        ScrollPaneConstants.VERTICAL_SCROLLBAR_AS_NEEDED);
    inputScrollPane.setHorizontalScrollBarPolicy(
        ScrollPaneConstants.HORIZONTAL_SCROLLBAR_AS_NEEDED);

    VScrollPanel cont = new VScrollPanel();
    cont.setLayout(new GridBagLayout());
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

      d2.setSize(d2.getWidth(), d.getHeight());

      chartSubpanels[i].setSize(d2);
      cont.add(chartSubpanels[i], contConstraints);
      contConstraints.gridy += 1;
    }

    inputScrollPane.getViewport().setView(cont);
    inputScrollPane.setVisible(true);
    inputScrollPane.setMinimumSize(minDim);

    this.add(inputScrollPane, constraints);
    constraints.gridy += 1;

    // set size so that the result pane isn't distorted on window launch
    Dimension dimension = inputScrollPane.getPreferredSize();
    dimension.setSize(dimension.getWidth() + 5, dimension.getHeight());
    this.setPreferredSize(dimension);

    leftSlider = new JSlider(0, SLIDER_MAX, 0);
    leftSlider.setEnabled(false);
    leftSlider.addChangeListener(this);

    rightSlider = new JSlider(0, SLIDER_MAX, SLIDER_MAX);
    rightSlider.setEnabled(false);
    rightSlider.addChangeListener(this);

    startDate = new EditableDateDisplayPanel();
    startDate.addChangeListener(this);
    endDate = new EditableDateDisplayPanel();
    endDate.addChangeListener(this);

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
    clearAll.setBackground(Color.RED.darker());
    clearAll.addActionListener(this);
    clearAll.setEnabled(false);

    int yOfClear = constraints.gridy;

    constraints.fill = GridBagConstraints.HORIZONTAL;
    constraints.gridy += 2;
    constraints.weighty = 0;
    constraints.weightx = 1;
    constraints.gridwidth = 1;
    constraints.gridx = 0;

    this.add(leftSlider, constraints);

    constraints.gridx += 1;
    this.add(rightSlider, constraints);

    constraints.gridx = 0;
    constraints.gridy += 1;
    this.add(startDate, constraints);
    constraints.gridx += 1;
    this.add(endDate, constraints);

    constraints.fill = GridBagConstraints.NONE;
    constraints.gridy += 1;
    constraints.weighty = 0;
    constraints.weightx = 0;

    constraints.gridx = 0;
    constraints.anchor = GridBagConstraints.EAST;
    this.add(zoomIn, constraints);

    constraints.gridx += 1;
    constraints.anchor = GridBagConstraints.WEST;
    this.add(zoomOut, constraints);

    constraints.gridwidth = 1;
    constraints.anchor = GridBagConstraints.CENTER;
    constraints.gridx = 7;
    constraints.gridy = yOfClear + 1;
    constraints.gridheight = 2;
    constraints.fill = GridBagConstraints.BOTH;

    this.add(save, constraints);

    constraints.gridy += 2;
    constraints.gridheight = GridBagConstraints.REMAINDER;

    this.add(clearAll, constraints);

    fileChooser = new JFileChooser();
    lastRespIndex = -1;

  }

  /**
   * Dispatches commands when interface buttons are clicked.
   * When the save button is clicked, dispatches the command to save plots as
   * an image. When the zoom buttons are clicked, scales the plot to only
   * contain the data within a specific range. Resets plots and removes
   * underlying data when the clear buttons are clicked. Prompts user to
   * load in a file for miniseed and resp data when the corresponding
   * buttons are clicked; because seed files can be large and contain a lot
   * of data to plot, runs seed-loading code backend in a separate thread.
   *
   * When new data is loaded in, this also fires a change event to any listeners
   */
  @Override
  public void actionPerformed(ActionEvent event) {

    for (int i = 0; i < FILE_COUNT; ++i) {
      JButton clear = clearButton[i];
      FileOperationJButton seed = seedLoaders[i];
      JButton resp = respLoaders[i];
      FileOperationJButton append = seedAppenders[i];

      if (event.getSource() == clear) {
        instantiateChart(i);
        dataStore.removeData(i);
        clear.setEnabled(false);
        append.setEnabled(false);
        seedFileNames[i].setText("NO FILE LOADED");
        respFileNames[i].setText("NO FILE LOADED");

        // plot all valid range of loaded-in data or else
        // disable clearAll button if there's no other data loaded in
        clearAll.setEnabled(dataStore.isAnythingSet());

        dataStore.trimToCommonTime();

        showRegionForGeneration();

        fireStateChanged();
      }

      if (event.getSource() == seed) {
        loadData(i, seed);
      }
      if (event.getSource() == append) {
        loadData(i, append);
      }

      if (event.getSource() == resp) {
        // don't need a new thread because resp loading is pretty prompt

        Set<String> respFilenames = InstrumentResponse.parseInstrumentList();

        List<String> names = new ArrayList<>(respFilenames);
        Collections.sort(names);
        names.add("Load custom response...");
        String[] nameArray = new String[names.size()];
        for (int k = 0; k < nameArray.length; ++k) {
          String name = names.get(k);
          name = name.replace("resps/","");
          nameArray[k] = name;
        }

        int index = lastRespIndex;
        if (lastRespIndex < 0) {
          index = names.size() - 1;
        }

        JDialog dialog = new JDialog();
        Object result = JOptionPane.showInputDialog(
            dialog,
            "Select a response to load:",
            "RESP File Selection",
            JOptionPane.PLAIN_MESSAGE,
            null, nameArray,
            nameArray[index]);

        final String resultStr = (String) result;

        // did user cancel operation?
        if (resultStr == null) {
          return;
        }

        String embeddedPath = "resps/" + resultStr;

        // is the loaded string one of the embedded response files?
        if (respFilenames.contains(embeddedPath)) {
          // what was the index of the selected item?
          // used to make sure we default to that choice next round
          lastRespIndex = Collections.binarySearch(names, embeddedPath);
          // final used here in the event of thread weirdness
          try {
            InstrumentResponse instrumentResponse =
                InstrumentResponse.loadEmbeddedResponse(embeddedPath);
            dataStore.setResponse(i, instrumentResponse);

            respFileNames[i].setText(instrumentResponse.getName());
            clear.setEnabled(true);
            clearAll.setEnabled(true);

            fireStateChanged();
          } catch (IOException e) {
            // this really shouldn't be an issue with embedded responses
            responseErrorPopup(resultStr);
            e.printStackTrace();
            return;
          }
        } else {
          lastRespIndex = -1;
          fileChooser.setCurrentDirectory(new File(respDirectory));
          fileChooser.resetChoosableFileFilters();
          fileChooser.setDialogTitle("Load response file...");
          int returnVal = fileChooser.showOpenDialog(resp);
          if (returnVal == JFileChooser.APPROVE_OPTION) {
            File file = fileChooser.getSelectedFile();
            respDirectory = file.getParent();
            try {
              Instant startInst = respEpochChoose(file.getAbsolutePath());
              InstrumentResponse instrumentResponse = new InstrumentResponse(file.getAbsolutePath(),
                  startInst);
              dataStore.setResponse(i, instrumentResponse);
              respFileNames[i].setText(file.getName());
              clear.setEnabled(true);
              clearAll.setEnabled(true);
            } catch (IOException e) {
              responseErrorPopup(file.getName());
              e.printStackTrace();
              return;
            }

            fireStateChanged();
          }
        }

        return;
      }

    }

    if (event.getSource() == clearAll) {
      clearAllData();
      fireStateChanged();
      return;
    }

    if (event.getSource() == save) {
      String ext = ".png";
      fileChooser = new JFileChooser();
      fileChooser.setCurrentDirectory(new File(saveDirectory));
      fileChooser.addChoosableFileFilter(
          new FileNameExtensionFilter("PNG image (.png)", ext));
      fileChooser.setFileFilter(fileChooser.getChoosableFileFilters()[1]);
      int returnVal = fileChooser.showSaveDialog(save);
      if (returnVal == JFileChooser.APPROVE_OPTION) {
        File selFile = fileChooser.getSelectedFile();
        saveDirectory = selFile.getParent();
        if (!selFile.getName().endsWith(ext.toLowerCase())) {
          selFile = new File(selFile.toString() + ext);
        }
        try {
          int height = IMAGE_HEIGHT * activePlots;
          BufferedImage image = getAsImage(IMAGE_WIDTH, height, activePlots);
          ImageIO.write(image, "png", selFile);
        } catch (IOException e) {
          e.printStackTrace();
        }
      }
      return;
    }

    if (event.getSource() == zoomIn) {
      showRegionForGeneration();
      return;
    }

    if (event.getSource() == zoomOut) {
      // restore original loaded dataStore
      dataStore.untrim(activePlots);
      for (int i = 0; i < FILE_COUNT; ++i) {
        if (!dataStore.blockIsSet(i)) {
          continue;
        }
        resetPlotZoom(i);
      }

      leftSlider.setValue(0);
      rightSlider.setValue(SLIDER_MAX);
      setVerticalBars();
      zoomOut.setEnabled(false);
    }
  }

  private void responseErrorPopup(String filename) {
    JDialog errorBox = new JDialog();
    String errorMsg = "Error while loading in response file:\n" + filename
        + "\nCheck that file exists and is formatted correctly.";
    JOptionPane.showMessageDialog(errorBox, errorMsg, "Response Loading Error",
        JOptionPane.ERROR_MESSAGE);
  }

  /**
   * Used to add objects to the list that will be informed when data is loaded
   */
  public void addChangeListener(ChangeListener listener) {
    listenerList.add(ChangeListener.class, listener);
  }

  /**
   * Resets the data and blanks out all charts
   */
  private void clearAllData() {
    dataStore = new DataStore();

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
      seedAppenders[i].setEnabled(false);
      instantiateChart(i);
    }

  }

  /**
   * Informs listening objects that the state of the inputs has changed
   * This is done when new seed or resp data has been loaded in, mainly
   * to tell whether enough data exists to run one of the experiments
   */
  private void fireStateChanged() {
    ChangeListener[] listeners = listenerList.getListeners(ChangeListener.class);
    if (listeners != null && listeners.length > 0) {
      ChangeEvent event = new ChangeEvent(this);
      for (ChangeListener listener : listeners) {
        listener.stateChanged(event);
      }
    }
  }

  /**
   * Return this panel's charts as a single buffered image
   * with specified dimensions
   *
   * @param width Width of returned image
   * @param height Height of returned image
   * @param plotsToShow Plots to be shown in the output image
   * @return Buffered image of the plots, writable to file
   */
  private BufferedImage getAsImage(int width, int height, int plotsToShow) {
    if (plotsToShow <= 0) {
      // should never be called like this but just in case
      // return an empty image
      return new BufferedImage(0, 0, BufferedImage.TYPE_INT_RGB);
    }

    int chartHeight = height / plotsToShow;

    JFreeChart[] chartsToPrint = new JFreeChart[plotsToShow];
    for (int i = 0; i < plotsToShow; ++i) {
      chartsToPrint[i] = chartPanels[i].getChart();
    }

    return ReportingUtils.chartsToImage(width, chartHeight, chartsToPrint);
  }


  /**
   * Produce the visible charts as multiple images, to be used to split
   * across multiple PDF pages when generating reports
   *
   * @param width Width of each plot in the image
   * @param height Height of each plot in the image
   * @param plotsToShow Number of plots that are to be added to the report
   * @return list of buffered images of plots, each image to be made a PDF page
   */
  public BufferedImage[]
  getAsMultipleImages(int width, int height, int plotsToShow) {
    if (plotsToShow <= 0) {
      return new BufferedImage[]{};
    }

    int chartHeight = height / plotsToShow;

    JFreeChart[] chartsToPrint = new JFreeChart[plotsToShow];
    for (int i = 0; i < plotsToShow; ++i) {
      chartsToPrint[i] = chartPanels[i].getChart();
    }

    return ReportingUtils.chartsToImageList(
        PLOTS_PER_PAGE, width, chartHeight, chartsToPrint);
  }

  /**
   * Returns the selected region of underlying DataStore, to be fed
   * into experiments for processing (the results of which will be plotted)
   * When this function is called, the graphs zoom to the currently active
   * range to display selected by the sliders, which is also the range
   * passed into an experiment
   *
   * @return A DataStore object (contains arrays of DataBlocks & Responses)
   */
  public DataStore getData() {
    if (dataStore.numberOfBlocksSet() > 1) {
      dataStore.trimToCommonTime(activePlots);
    }
    return new DataStore(dataStore);
  }

  /**
   * Gets the height of resulting image of plots given default parameters,
   * so that it only needs to fit the plots that have data in them
   *
   * @return height of image to output
   */
  public int getImageHeight(int plotsToShow) {
    return IMAGE_HEIGHT * plotsToShow;
  }

  public String[] getResponseStrings(int[] indices) {
    String[] outStrings = new String[indices.length];
    for (int i = 0; i < indices.length; ++i) {
      int idx = indices[i];
      if (!dataStore.responseIsSet(idx)) {
        System.out.println("ERROR WITH RESP AT INDEX " + idx);
      }
      outStrings[i] = dataStore.getResponse(idx).toString();
    }
    return outStrings;
  }

  /**
   * Instantiates the underlying chart of a chartpanel with default data
   *
   * @param index Index of the chartpanel to instantiate
   */
  private void instantiateChart(int index) {
    JFreeChart chart = ChartFactory.createXYLineChart(
        "SEED input " + (index + 1),
        "Time",
        "Counts",
        new XYSeriesCollection(),
        PlotOrientation.VERTICAL,
        false, false, false);

    if (chartPanels[index] == null) {
      chartPanels[index] = new ChartPanel(chart);
    } else {
      chartPanels[index].setChart(chart);
    }
    chartPanels[index].setMouseZoomable(true);
  }

  /**
   * Load in data for a specified SEED file, to be run in a specific thread.
   * Because loading can be a slow operation, this runs in a background thread
   *
   * @param index Index into datastore/plots this data should be loaded
   * @param seed The JButton to passed into the file loader
   */
  private void loadData(final int index, final FileOperationJButton seed) {

    fileChooser.setCurrentDirectory(new File(seedDirectory));
    fileChooser.resetChoosableFileFilters();
    fileChooser.setDialogTitle("Load SEED file...");
    int returnVal = fileChooser.showOpenDialog(seed);
    if (returnVal == JFileChooser.APPROVE_OPTION) {
      final File file = fileChooser.getSelectedFile();
      seedDirectory = file.getParent();
      String oldName = seedFileNames[index].getText();

      seedFileNames[index].setText("LOADING: " + file.getName());
      final String filePath = file.getAbsolutePath();
      String filterName;
      try {
        Set<String> nameSet = seed.getFilenameSet(dataStore, index, filePath);

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
            // if the user cancelled selecting a subseries
            seedFileNames[index].setText(oldName);
            return;
          }
        } else if (nameSet.size() > 0) {
          // just get the first one; it's the only one in the list
          filterName = new ArrayList<>(nameSet).get(0);
        } else {
          // this shouldn't occur if the seed file is formatted correctly
          // because there will be at least one SNCL described in it
          // unless we are appending from a list without data for the loaded channel, etc.
          seedFileNames[index].setText(oldName);
          seedAppendEmptyPopup(file.getName());
          return;
        }

      } catch (SeedFormatException | IOException | RuntimeException e) {
        e.printStackTrace();
        if (seed instanceof LoadingJButton) {
          seedLoadHandleError(index, file.getName(), e.toString());
        } else {
          seedAppendErrorPopup(file.getName());
        }

        return;
      }

      final String immutableFilter = filterName;

      // create swingworker to load large files in the background
      SwingWorker<Integer, Void> worker = new SwingWorker<Integer, Void>() {

        JFreeChart chart;
        boolean caughtException = false;
        String returnedErrMsg = "";

        @Override
        public Integer doInBackground() {

          try {
            seed.loadInData(dataStore, index, filePath, immutableFilter, activePlots);
          } catch (RuntimeException | SeedFormatException | CodecException |
              FileNotFoundException e) {
            returnedErrMsg = e.toString();
            caughtException = true;
            e.printStackTrace();
            return 1;
          }

          dataStore.untrim(activePlots);

          XYSeries timeSeries = dataStore.getBlock(index).toXYSeries();
          String rateString = " (" + dataStore.getBlock(index).getSampleRate() + " Hz)";
          chart = ChartFactory.createXYLineChart(
              timeSeries.getKey().toString() + rateString,
              "Time",
              "Counts",
              new XYSeriesCollection(timeSeries),
              PlotOrientation.VERTICAL,
              false, false, false);

          XYPlot xyPlot = (XYPlot) chart.getPlot();

          DateAxis dateAxis = new DateAxis();
          dateAxis.setLabel("UTC Time (Year.Day.Hour:Minute)");
          Font bold = dateAxis.getLabelFont();
          bold = bold.deriveFont(Font.BOLD);
          dateAxis.setLabelFont(bold);
          dateAxis.setDateFormatOverride(ExperimentPanel.DATE_TIME_FORMAT.get());
          xyPlot.setDomainAxis(dateAxis);
          int colorIndex = index % defaultColor.length;
          xyPlot.getRenderer().setSeriesPaint(0, defaultColor[colorIndex]);

          return 0;
        }

        @Override
        public void done() {
          if (caughtException) {
            if (seed instanceof AppendingJButton) {
              seedAppendErrorPopup(file.getName());
              return;
            }
            seedLoadHandleError(index, file.getName(), returnedErrMsg);
            return;
          }

          NumberAxis numberAxis = (NumberAxis) chart.getXYPlot().getRangeAxis();
          numberAxis.setAutoRange(true);
          numberAxis.setAutoRangeIncludesZero(false);
          chartPanels[index].setChart(chart);
          chartPanels[index].repaint();
          chartPanels[index].setMouseZoomable(true);

          seedAppenders[index].setEnabled(true);
          clearButton[index].setEnabled(true);

          for (int i = 0; i < FILE_COUNT; ++i) {
            if (!dataStore.blockIsSet(i)) {
              continue;
            }
            resetPlotZoom(i);
          }

          leftSlider.setValue(0);
          rightSlider.setValue(SLIDER_MAX);
          setVerticalBars();

          zoomOut.setEnabled(false);
          zoomIn.setEnabled(true);
          leftSlider.setEnabled(true);
          rightSlider.setEnabled(true);
          save.setEnabled(true);
          clearAll.setEnabled(true);

          seedFileNames[index].setText(file.getName() + ": " + immutableFilter);

          fireStateChanged();
        }
      };

      worker.execute();
    }
  }

  private void seedLoadHandleError(int index, String filename, String error) {
    dataStore.removeBlock(index);
    instantiateChart(index);
    XYPlot xyPlot = (XYPlot) chartPanels[index].getChart().getPlot();
    TextTitle result = new TextTitle();
    result.setText("COULD NOT LOAD IN FILE: " + filename
        + "\n\nThis file is probably not a SEED file or has a formatting error.\n\n"
        + "A full stack trace is in the terminal -- this is the exception:\n"
        + error);
    result.setBackgroundPaint(Color.red);
    result.setPaint(Color.white);
    XYTitleAnnotation titleAnnotation = new XYTitleAnnotation(0.5, 0.5, result,
        RectangleAnchor.CENTER);
    xyPlot.clearAnnotations();
    xyPlot.addAnnotation(titleAnnotation);
    seedFileNames[index].setText("NO FILE LOADED");
    clearButton[index].setEnabled(true);
    fireStateChanged();
  }

  private void seedAppendEmptyPopup(String filename) {
    JDialog errorBox = new JDialog();
    String errorMsg = "Could not load data from seed file: " + filename
        + '\n'
        + "The file appears to be formatted correctly but does not have data for the\n"
        + "SNCL data currently loaded in. Check that the correct file was chosen.\n"
        + "(No data was appended to the currently loaded data.)";
    JOptionPane.showMessageDialog(errorBox, errorMsg, "Response Loading Error",
        JOptionPane.ERROR_MESSAGE);
  }

  private void seedAppendErrorPopup(String filename) {
    JDialog errorBox = new JDialog();
    String errorMsg = "Error while loading in seed file: " + filename
        + "\nCheck that file exists and is formatted correctly.\n"
        + "(No data was appended to the currently loaded data.)";
    JOptionPane.showMessageDialog(errorBox, errorMsg, "Response Loading Error",
        JOptionPane.ERROR_MESSAGE);
  }

  /**
   * Used to construct the panels for loading and displaying SEED data
   * (as well as the corresponding response file)
   *
   * @param index Index of panel to be created, for getting references to the chart
   * panel and the appropriate actionlisteners for the loaders
   * @return composite panel of chart, loaders, and clear button
   */
  private JPanel makeChartSubpanel(int index) {

    JPanel chartSubpanel = new JPanel();
    chartSubpanel.setLayout(new GridBagLayout());
    GridBagConstraints constraints = new GridBagConstraints();

    constraints.weightx = 1.0;
    constraints.weighty = 1.0;
    constraints.gridy = 0;
    constraints.anchor = GridBagConstraints.CENTER;

    instantiateChart(index);

    channelType[index] = new JLabel("");

    chartPanels[index].setMouseZoomable(true);

    seedLoaders[index] = new LoadingJButton("Load SEED file " + (index + 1));
    seedLoaders[index].addActionListener(this);
    seedLoaders[index].setMaximumSize(seedLoaders[index].getMinimumSize());

    seedAppenders[index] = new AppendingJButton("Append SEED");
    seedAppenders[index].addActionListener(this);
    seedAppenders[index].setMaximumSize(seedLoaders[index].getMinimumSize());
    seedAppenders[index].setEnabled(false);

    JTextField text = new JTextField("NO FILE LOADED");
    text.setHorizontalAlignment(SwingConstants.CENTER);
    text.setMaximumSize(text.getPreferredSize());
    seedFileNames[index] = text;
    seedFileNames[index].setEditable(false);

    respLoaders[index] = new JButton("Load RESP file " + (index + 1));
    respLoaders[index].addActionListener(this);
    respLoaders[index].setMaximumSize(respLoaders[index].getMinimumSize());

    text = new JTextField("NO FILE LOADED");
    text.setHorizontalAlignment(SwingConstants.CENTER);
    text.setMaximumSize(text.getPreferredSize());
    respFileNames[index] = text;
    respFileNames[index].setEditable(false);

    clearButton[index] = new JButton("Clear data " + (index + 1));
    clearButton[index].setMaximumSize(clearButton[index].getMinimumSize());

    constraints.gridx = 0;
    constraints.gridy = 0;

    constraints.weightx = 0;
    constraints.weighty = 0;
    constraints.fill = GridBagConstraints.BOTH;

    constraints.fill = GridBagConstraints.NONE;
    chartSubpanel.add(channelType[index], constraints);

    constraints.weightx = 1;
    constraints.weighty = 1;
    constraints.fill = GridBagConstraints.BOTH;
    constraints.gridwidth = 1;
    constraints.gridheight = 6;
    constraints.gridy += 1;
    chartSubpanel.add(chartPanels[index], constraints);

    constraints.fill = GridBagConstraints.BOTH;
    constraints.gridx = 1;
    constraints.gridwidth = 1;
    constraints.gridheight = 1;
    constraints.weightx = 0;
    constraints.weighty = 0.25;
    chartSubpanel.add(seedLoaders[index], constraints);

    constraints.gridy += 1;
    chartSubpanel.add(seedAppenders[index], constraints);

    constraints.fill = GridBagConstraints.BOTH;
    constraints.weighty = 1;
    constraints.gridy += 1;
    JScrollPane scrollPane = new JScrollPane();
    scrollPane.setMaximumSize(scrollPane.getMinimumSize());
    scrollPane.setViewportView(seedFileNames[index]);
    scrollPane.setVerticalScrollBarPolicy(
        ScrollPaneConstants.VERTICAL_SCROLLBAR_NEVER);
    scrollPane.setHorizontalScrollBarPolicy(
        ScrollPaneConstants.HORIZONTAL_SCROLLBAR_ALWAYS);
    chartSubpanel.add(scrollPane, constraints);

    constraints.fill = GridBagConstraints.BOTH;
    constraints.weighty = 0.25;
    constraints.gridy += 1;
    chartSubpanel.add(respLoaders[index], constraints);

    constraints.fill = GridBagConstraints.BOTH;
    constraints.weighty = 1;
    constraints.gridy += 1;
    scrollPane = new JScrollPane();
    scrollPane.setMaximumSize(scrollPane.getMinimumSize());
    scrollPane.setViewportView(respFileNames[index]);
    scrollPane.setVerticalScrollBarPolicy(
        ScrollPaneConstants.VERTICAL_SCROLLBAR_NEVER);
    scrollPane.setHorizontalScrollBarPolicy(
        ScrollPaneConstants.HORIZONTAL_SCROLLBAR_ALWAYS);
    chartSubpanel.add(scrollPane, constraints);

    constraints.fill = GridBagConstraints.HORIZONTAL;
    constraints.weighty = 0;
    constraints.gridy += 1;
    clearButton[index] = new JButton("Clear data " + (index + 1));
    clearButton[index].setOpaque(true);
    clearButton[index].setBackground(Color.RED.darker());
    clearButton[index].addActionListener(this);
    clearButton[index].setEnabled(false);
    chartSubpanel.add(clearButton[index], constraints);

    return chartSubpanel;
  }

  /**
   * Does the work to reset the zoom of a chart when the zoom button is hit
   *
   * @param index Index of appropriate chart/panel
   */
  private void resetPlotZoom(int index) {
    XYPlot xyPlot = chartPanels[index].getChart().getXYPlot();
    XYSeriesCollection timeSeries = new XYSeriesCollection();
    timeSeries.addSeries(dataStore.getBlock(index).toXYSeries());
    xyPlot.setDataset(timeSeries);
    xyPlot.getRenderer().setSeriesPaint(0,
        defaultColor[index % defaultColor.length]);
    xyPlot.getDomainAxis().setAutoRange(true);
    if (xyPlot.getSeriesCount() > 1) {
      throw new RuntimeException("TOO MUCH DATA");
    }
    chartPanels[index].repaint();
  }

  /**
   * Get a selected epoch from a multi-epoch response
   *
   * @param respHandle Name of a given response to read in
   * @return Value of the start instant of the epoch that the user has chosen
   * @throws IOException Error reading resp file
   * @throws FileNotFoundException File does not exist
   */
  private Instant respEpochChoose(String respHandle) throws IOException {

    DateTimeFormatter dateTimeFormatter = InstrumentResponse.RESP_DT_FORMAT
        .withZone(ZoneOffset.UTC);
    List<Pair<Instant, Instant>> epochs = InstrumentResponse.getRespFileEpochs(respHandle);

    if (epochs.size() > 1) {
      // more than one series in the file? prompt user for it
      String[] epochStrings = new String[epochs.size()];
      for (int i = 0; i < epochStrings.length; ++i) {
        String startString = dateTimeFormatter.format(epochs.get(i).getFirst());
        Instant endInstant = epochs.get(i).getSecond();
        String endString;
        if (endInstant != null) {
          endString = dateTimeFormatter.format(epochs.get(i).getSecond());
        } else {
          endString = "(NO EPOCH END SPECIFIED)";
        }
        epochStrings[i] = startString + " | " + endString;
      }
      Arrays.sort(epochStrings);
      JDialog dialog = new JDialog();
      Object result = JOptionPane.showInputDialog(
          dialog,
          "Select the response epoch to load:",
          "Response Epoch Selection",
          JOptionPane.PLAIN_MESSAGE,
          null, epochStrings,
          epochStrings[0]);
      if (result instanceof String) {
        int index = Arrays.binarySearch(epochStrings, result);
        return epochs.get(index).getFirst();
      } else {
        return null;
      }
    } else if (epochs.size() > 0) {
      return epochs.get(0).getFirst();
    } else {
      throw new IOException("RESP file has no epoch data -- check formatting");
    }

  }

  /**
   * Used to get labels for each plot to identify what data they need to
   * contain in order for an experiment to have enough data to run
   *
   * @param channels List of strings to be used as panel title
   */
  public void setChannelTypes(String[] channels) {

    int length = Math.min(channels.length, channelType.length);

    for (int i = 0; i < length; ++i) {
      channelType[i].setText(channels[i]);
      channelType[i].setHorizontalAlignment(SwingConstants.CENTER);
    }
  }

  /**
   * Displays the range set by the sliders using
   * vertical bars at the min and max values
   */
  private void setVerticalBars() {

    if (dataStore.numberOfBlocksSet() < 1) {
      return;
    }

    int leftValue = leftSlider.getValue();
    int rightValue = rightSlider.getValue();
    DataBlock dataBlock = dataStore.getXthLoadedBlock(1);
    long startMarkerLocation = getMarkerLocation(dataBlock, leftValue);
    long endMarkerLocation = getMarkerLocation(dataBlock, rightValue);

    startDate.removeChangeListener(this);
    endDate.removeChangeListener(this);
    startDate.setValues(startMarkerLocation);
    endDate.setValues(endMarkerLocation);
    startDate.addChangeListener(this);
    endDate.addChangeListener(this);

    for (int i = 0; i < FILE_COUNT; ++i) {
      if (!dataStore.blockIsSet(i)) {
        continue;
      }

      XYPlot xyPlot = chartPanels[i].getChart().getXYPlot();
      xyPlot.clearDomainMarkers();

      Marker startMarker = new ValueMarker(startMarkerLocation);
      startMarker.setStroke(new BasicStroke((float) 1.5));
      Marker endMarker = new ValueMarker(endMarkerLocation);
      endMarker.setStroke(new BasicStroke((float) 1.5));

      xyPlot.addDomainMarker(startMarker);
      xyPlot.addDomainMarker(endMarker);

      List<Pair<Long, Long>> gaps = dataStore.getBlock(i).getGapBoundaries();

      XYDataset data = xyPlot.getDataset();
      XYSeriesCollection timeSeries = (XYSeriesCollection) data;
      double min = timeSeries.getDomainLowerBound(false);
      double max = timeSeries.getDomainUpperBound(false);

      for (Pair<Long, Long> gapLocation : gaps) {
        Double gapStart = gapLocation.getFirst().doubleValue();
        Double gapEnd = gapLocation.getSecond().doubleValue();
        if (gapEnd > min || gapStart < max) {
          double start = Math.max(gapStart, min);
          double end = Math.min(gapEnd, max);
          Marker gapMarker = new IntervalMarker(start, end);
          gapMarker.setPaint(Color.ORANGE);
          xyPlot.addDomainMarker(gapMarker);
        }
      }
      chartPanels[i].repaint();
    }

  }

  /**
   * Show the number of panels needed to load in data for a specific experiment
   *
   * @param panelsNeeded Number of panels to show
   */
  public void showDataNeeded(int panelsNeeded) {

    VScrollPanel cont = new VScrollPanel();
    cont.setLayout(new GridBagLayout());
    GridBagConstraints constraints = new GridBagConstraints();
    constraints.weightx = 1.0;
    constraints.weighty = 1.0;
    constraints.gridy = 0;
    constraints.anchor = GridBagConstraints.CENTER;
    constraints.fill = GridBagConstraints.BOTH;

    activePlots = panelsNeeded;

    // get current time range of zoom data for resetting, if any data is loaded
    long start, end;
    if (dataStore.areAnyBlocksSet()) {
      DataBlock dataBlock = dataStore.getXthLoadedBlock(1);
      start = dataBlock.getStartTime();
      end = dataBlock.getEndTime();

      dataStore.trimToCommonTime(activePlots);
      // try to trim to current active time range if possible, otherwise fit
      // as much data as possible
      dataBlock = dataStore.getXthLoadedBlock(1);
      // was the data zoomed in more than it is now?
      if (start > dataBlock.getStartTime() || end < dataBlock.getEndTime()) {
        try {
          // zooms won't be modified if an exception is thrown
          dataStore.trim(start, end, activePlots);
          zoomOut.setEnabled(true);
        } catch (IndexOutOfBoundsException e) {
          // new time range not valid for all current data, show max range
          zoomOut.setEnabled(false);
        }
      } else {
        // common time range was already the max
        zoomOut.setEnabled(false);
      }

    } else {
      // no blocks loaded in, no zooms to handle
      zoomOut.setEnabled(false);
    }

    // reset plot zoom if data is still there; clear out stale data from conflicting time ranges
    for (int i = 0; i < activePlots; ++i) {
      if (dataStore.blockIsSet(i)) {
        resetPlotZoom(i);
      } else {
        instantiateChart(i);
        seedFileNames[i].setText("NO FILE LOADED");
        if (!dataStore.responseIsSet(i)) {
          respFileNames[i].setText("NO FILE LOADED");
          clearButton[i].setEnabled(false);
        }
      }

      cont.add(chartSubpanels[i], constraints);
      constraints.gridy += 1;
    }

    leftSlider.setValue(0);
    rightSlider.setValue(SLIDER_MAX);
    setVerticalBars();

    zoomIn.setEnabled(dataStore.numberOfBlocksSet() > 0);

    // using this test means the panel doesn't try to scroll when it's
    // only got a few inputs to deal with, when stuff is still pretty readable
    cont.setScrollableTracksViewportHeight(activePlots <= MAX_UNSCROLLED);

    inputScrollPane.getViewport().setView(cont);
    inputScrollPane.setPreferredSize(cont.getPreferredSize());
  }

  /**
   * Zooms in on the current range of data, which will be passed into
   * backend functions for experiment calculations
   */
  public void showRegionForGeneration() {

    if (dataStore.numberOfBlocksSet() < 1) {
      return;
    }

    // get (any) loaded data block to map slider to domain boundary
    // all data should have the same range
    DataBlock dataBlock = dataStore.getXthLoadedBlock(1);

    if (leftSlider.getValue() != 0 || rightSlider.getValue() != SLIDER_MAX) {
      long start = getMarkerLocation(dataBlock, leftSlider.getValue());
      long end = getMarkerLocation(dataBlock, rightSlider.getValue());
      dataStore.trim(start, end, activePlots);
      leftSlider.setValue(0);
      rightSlider.setValue(SLIDER_MAX);
      zoomOut.setEnabled(true);
    }

    for (int i = 0; i < activePlots; ++i) {
      if (!dataStore.blockIsSet(i)) {
        continue;
      }
      resetPlotZoom(i);

    }

    setVerticalBars();

  }

  /**
   * Handles changes in value by the sliders below the charts
   **/
  @Override
  public void stateChanged(ChangeEvent event) {

    int leftSliderValue = leftSlider.getValue();
    int rightSliderValue = rightSlider.getValue();

    if (event.getSource() == startDate) {
      // if no data to do windowing on, don't bother
      if (dataStore.numberOfBlocksSet() < 1) {
        return;
      }

      long time = startDate.getTime();
      DataBlock dataBlock = dataStore.getXthLoadedBlock(1);

      long startTime = dataBlock.getStartTime();
      // startValue is current value of left-side slider in ms

      // assume current locations of sliders is valid
      int marginValue = rightSliderValue - MARGIN;
      long marginTime =
          getMarkerLocation(dataBlock, marginValue);

      // fix boundary cases
      if (time < startTime) {
        time = startTime;
      } else if (time > marginTime) {
        time = marginTime;
      }

      startDate.setValues(time);
      int newLeftSliderValue = getSliderValue(dataBlock, time);
      leftSlider.removeChangeListener(this);
      leftSlider.setValue(newLeftSliderValue); // already validated
      leftSlider.addChangeListener(this);
      setVerticalBars();
      return;
    }

    if (event.getSource() == endDate) {
      // if no data to do windowing on, don't bother
      if (dataStore.numberOfBlocksSet() < 1) {
        return;
      }

      long time = endDate.getTime();
      DataBlock dataBlock = dataStore.getXthLoadedBlock(1);

      long endTime = dataBlock.getEndTime();

      int marginValue = leftSliderValue + MARGIN;
      long marginTime = getMarkerLocation(dataBlock, marginValue);

      // fix boundary cases
      if (time > endTime) {
        time = endTime;
      } else if (time < marginTime) {
        time = marginTime;
      }

      endDate.setValues(time);
      int newRightSliderValue = getSliderValue(dataBlock, time);
      rightSlider.removeChangeListener(this);
      rightSlider.setValue(newRightSliderValue); // already validated
      rightSlider.addChangeListener(this);
      setVerticalBars();
      return;
    }

    if (event.getSource() == leftSlider) {
      validateSliderPlacement(true, leftSliderValue);
    } else if (event.getSource() == rightSlider) {
      validateSliderPlacement(false, rightSliderValue);
    }

    setVerticalBars(); // date display object's text gets updated here

  }

  /**
   * Verify that slider locations will not violate restrictions in location
   *
   * @param moveLeft True if left slider needs to move (false if right slider)
   * @param newLocation Value to set slider to if within restrictions
   */
  private void validateSliderPlacement(boolean moveLeft, int newLocation) {

    int leftSliderValue, rightSliderValue;

    if (moveLeft) {
      leftSliderValue = newLocation;
      rightSliderValue = rightSlider.getValue();
    } else {
      leftSliderValue = leftSlider.getValue();
      rightSliderValue = newLocation;
    }

    if (leftSliderValue > rightSliderValue ||
        leftSliderValue + MARGIN > rightSliderValue) {

      // (left slider must stay left of right slider by at least margin)

      if (moveLeft) {
        // move left slider as close to right as possible
        leftSliderValue = rightSliderValue - MARGIN;
        if (leftSliderValue < 0) {
          leftSliderValue = 0;
          rightSliderValue = MARGIN;
        }
      } else {
        // move right slider as close to left as possible
        rightSliderValue = leftSliderValue + MARGIN;
        if (rightSliderValue > SLIDER_MAX) {
          rightSliderValue = SLIDER_MAX;
          leftSliderValue = SLIDER_MAX - MARGIN;
        }
      }

    }

    rightSlider.setValue(rightSliderValue);
    leftSlider.setValue(leftSliderValue);

  }

  private abstract class FileOperationJButton extends JButton {

    private static final long serialVersionUID = -8181485763533504906L;

    FileOperationJButton(String text) {
      super(text);
    }

    protected abstract Set<String> getFilenameSet(DataStore ds, int index, String filePath)
        throws SeedFormatException, IOException;

    protected abstract void loadInData(DataStore dataStore, int index,
        String filePath, String fileFilter, int activePlots)
        throws SeedFormatException, CodecException,
        FileNotFoundException, RuntimeException;
  }

  private class LoadingJButton extends FileOperationJButton {

    private static final long serialVersionUID = -5122928259743764418L;

    LoadingJButton(String text) {
      super(text);
    }

    @Override
    public Set<String> getFilenameSet(DataStore dataStore, int index, String filePath)
        throws SeedFormatException, IOException {
      return TimeSeriesUtils.getMplexNameSet(filePath);
    }

    @Override
    public void loadInData(DataStore dataStore, int index, String filePath,
        String fileFilter, int activePlots) throws SeedFormatException,
        CodecException,
        FileNotFoundException,
        RuntimeException {
      dataStore.setBlock(index, filePath, fileFilter, activePlots);
    }
  }

  private class AppendingJButton extends FileOperationJButton {

    private static final long serialVersionUID = -2300184562192345233L;

    AppendingJButton(String text) {
      super(text);
    }

    @Override
    public Set<String> getFilenameSet(DataStore dataStore, int index, String filePath)
        throws SeedFormatException, IOException {
      String thisName = dataStore.getBlock(index).getName();
      if (!TimeSeriesUtils.getMplexNameSet(filePath).contains(thisName)) {
        return new HashSet<>();
      }
      Set<String> returnSet = new HashSet<>();
      returnSet.add(thisName);
      return returnSet;
    }

    @Override
    public void loadInData(DataStore dataStore, int index, String filePath,
        String fileFilter, int activePlots)
        throws SeedFormatException, CodecException, FileNotFoundException, RuntimeException {
      dataStore.appendBlock(index, filePath, fileFilter, activePlots);
    }
  }

}

