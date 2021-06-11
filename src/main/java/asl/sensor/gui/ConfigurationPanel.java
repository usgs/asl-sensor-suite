package asl.sensor.gui;

import asl.sensor.input.Configuration;
import java.awt.GridLayout;
import javax.swing.JCheckBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.text.AbstractDocument;
import javax.swing.text.AttributeSet;
import javax.swing.text.BadLocationException;
import javax.swing.text.DocumentFilter;

/**
 * Panel that displays key-value pairs for the configuration.
 * Presented as a popup from the main window, which controls saving the data after it is edited.
 */
public class ConfigurationPanel extends JPanel {

  // text fields for file output locations
  private final JTextField dataFolder, respFolder, noiseModelFolder, outputFolder;

  // text fields for FDSN data locations
  private final JTextField fdsnDomain, fdsnProtocol, fdsnPort, fdsnService;
  // text fields for FDSN meta locations
  private final JTextField fdsnMetaDomain, fdsnMetaProtocol, fdsnMetaPort, fdsnMetaService;

  // fields for controlling plot display parameters
  private final JTextField lineWidthOffset;
  private final JCheckBox colorblindColors;

  /**
   * Construct a panel allowing editing of the current configuration parameters.
   */
  public ConfigurationPanel() {
    // populate fields based on current configuration file
    Configuration instance = Configuration.getInstance();
    PositiveIntegerFilter filter = new PositiveIntegerFilter();


    dataFolder = new JTextField();
    dataFolder.setText(instance.getDefaultDataFolder());
    respFolder = new JTextField();
    respFolder.setText(instance.getDefaultRespFolder());
    noiseModelFolder = new JTextField();
    noiseModelFolder.setText(instance.getDefaultNoiseModelFolder());
    outputFolder = new JTextField();
    outputFolder.setText(instance.getDefaultOutputFolder());

    fdsnProtocol = new JTextField();
    fdsnProtocol.setText(instance.getFDSNDataProtocol());
    fdsnDomain = new JTextField();
    fdsnDomain.setText(instance.getFDSNDataDomain());
    fdsnService = new JTextField();
    fdsnService.setText(instance.getFDSNDataPath());
    // enforce port to be a positive integer
    fdsnPort = new JTextField();
    AbstractDocument intDocument = (AbstractDocument) fdsnPort.getDocument();
    intDocument.setDocumentFilter(filter);
    fdsnPort.setText(Integer.toString(instance.getFDSNDataPort()));

    fdsnMetaProtocol = new JTextField();
    fdsnMetaProtocol.setText(instance.getFDSNMetaProtocol());
    fdsnMetaDomain = new JTextField();
    fdsnMetaDomain.setText(instance.getFDSNMetaDomain());
    fdsnMetaService = new JTextField();
    fdsnMetaService.setText(instance.getFDSNMetaPath());
    // enforce port to be a positive integer
    fdsnMetaPort = new JTextField();
    intDocument = (AbstractDocument) fdsnMetaPort.getDocument();
    intDocument.setDocumentFilter(filter);
    fdsnMetaPort.setText(Integer.toString(instance.getFDSNMetaPort()));

    lineWidthOffset = new JTextField();
    intDocument = (AbstractDocument) lineWidthOffset.getDocument();
    intDocument.setDocumentFilter(filter);
    lineWidthOffset.setText(Integer.toString(instance.getLineWidthOffset()));

    colorblindColors = new JCheckBox();
    colorblindColors.setEnabled(true);
    colorblindColors.setSelected(instance.useColorblindColors());

    this.setLayout(new GridLayout(14, 2));

    this.add(new JLabel("Default SEED location:"));
    this.add(dataFolder);
    this.add(new JLabel("Default RESP location:"));
    this.add(respFolder);
    this.add(new JLabel("Default noise model location:"));
    this.add(noiseModelFolder);
    this.add(new JLabel("Default report location:"));
    this.add(outputFolder);

    this.add(new JLabel("FDSN query protocol (i.e., http, https):"));
    this.add(fdsnProtocol);
    this.add(new JLabel("FDSN URL domain:"));
    this.add(fdsnDomain);
    this.add(new JLabel("FDSN service subdomain:"));
    this.add(fdsnService);
    this.add(new JLabel("FDSN service port:"));
    this.add(fdsnPort);

    this.add(new JLabel("FDSN metadata query protocol (i.e., http, https):"));
    this.add(fdsnMetaProtocol);
    this.add(new JLabel("FDSN metadata URL domain:"));
    this.add(fdsnMetaDomain);
    this.add(new JLabel("FDSN metadata service subdomain:"));
    this.add(fdsnMetaService);
    this.add(new JLabel("FDSN metadata service port:"));
    this.add(fdsnMetaPort);

    this.add(new JLabel("Line thickness increase (int):"));
    this.add(lineWidthOffset);
    this.add(new JLabel("Plot colorblind-friendly colors:"));
    this.add(colorblindColors);
  }

  /**
   * Save the values of the current configuration to file.
   */
  public void writeValues() {
    Configuration instance = Configuration.getInstance();

    instance.setDefaultDataFolder(dataFolder.getText());
    instance.setDefaultRespFolder(respFolder.getText());
    instance.setDefaultNoiseModelFolder(noiseModelFolder.getText());
    instance.setDefaultOutputFolder(outputFolder.getText());

    instance.setFDSNDataProtocol(fdsnProtocol.getText());
    instance.setFDSNDataDomain(fdsnDomain.getText());
    instance.setFDSNDataPath(fdsnService.getText());
    instance.setFDSNDataPort(Integer.parseInt(fdsnPort.getText()));

    instance.setFDSNMetaProtocol(fdsnMetaProtocol.getText());
    instance.setFDSNMetaDomain(fdsnMetaDomain.getText());
    instance.setFDSNMetaPath(fdsnMetaService.getText());
    instance.setFDSNMetaPort(Integer.parseInt(fdsnMetaPort.getText()));

    instance.setLineWidthOffset(Integer.parseInt(lineWidthOffset.getText()));
    instance.setUseColorblindColors(colorblindColors.isSelected());

    instance.saveCurrentConfig();
  }

  /**
   * DocumentFilter implementation to force positive integers in relevant text fields
   */
  private static class PositiveIntegerFilter extends DocumentFilter {

    @Override
    public void insertString(DocumentFilter.FilterBypass fb, int offset, String text,
        AttributeSet attr) throws BadLocationException {
      try {
        int test = Integer.parseInt(text);
        test = Math.abs(test);
        fb.insertString(offset, String.valueOf(test), attr);
      } catch (NumberFormatException e) {
        // suppress output to terminal
        // (if user presses '.' that causes the parser to fail and print a stack trace)
      }
    }

    @Override
    public void replace(DocumentFilter.FilterBypass fb, int offset, int length, String text,
        AttributeSet attrs) throws BadLocationException {
      try {
        int test = Integer.parseInt(text);
        test = Math.abs(test);
        fb.replace(offset, length, String.valueOf(test), attrs);
      } catch (NumberFormatException e) {
        // suppress output to terminal as above
      }
    }
  }

}
