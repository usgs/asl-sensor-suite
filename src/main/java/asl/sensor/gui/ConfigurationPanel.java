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

public class ConfigurationPanel extends JPanel {

  // text fields for file output locations
  private JTextField dataFolder, respFolder, outputFolder;

  // text fields for FDSN data locations
  private JTextField fdsnDomain, fdsnProtocol, fdsnPort, fdsnService;

  // fields for controlling plot display parameters
  private JTextField lineWidthOffset;
  private JCheckBox colorblindColors;

  public ConfigurationPanel() {
    // populate fields based on current configuration file
    Configuration instance = Configuration.getInstance();

    dataFolder = new JTextField();
    dataFolder.setText(instance.getDefaultDataFolder());
    respFolder = new JTextField();
    respFolder.setText(instance.getDefaultRespFolder());
    outputFolder = new JTextField();
    outputFolder.setText(instance.getDefaultOutputFolder());

    fdsnProtocol = new JTextField();
    fdsnProtocol.setText(instance.getFDSNProtocol());
    fdsnDomain = new JTextField();
    fdsnDomain.setText(instance.getFDSNDomain());
    fdsnService = new JTextField();
    fdsnService.setText(instance.getFDSNPath());
    // enforce port to be a positive integer
    fdsnPort = new JTextField();
    AbstractDocument intDocument = (AbstractDocument) fdsnPort.getDocument();
    intDocument.setDocumentFilter(new PositiveIntegerFilter());
    fdsnPort.setText(Integer.toString(instance.getFDSNPort()));

    lineWidthOffset = new JTextField();
    intDocument = (AbstractDocument) lineWidthOffset.getDocument();
    intDocument.setDocumentFilter(new PositiveIntegerFilter());
    lineWidthOffset.setText(Integer.toString(instance.getLineWidthOffset()));

    colorblindColors = new JCheckBox();
    colorblindColors.setEnabled(true);
    colorblindColors.setSelected(instance.useColorblindColors());

    this.setLayout(new GridLayout(9, 2));

    this.add(new JLabel("Default SEED location:"));
    this.add(dataFolder);
    this.add(new JLabel("Default RESP location:"));
    this.add(respFolder);
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

    this.add(new JLabel("Line thickness increase (int):"));
    this.add(lineWidthOffset);
    this.add(new JLabel("Plot colorblind-friendly colors:"));
    this.add(colorblindColors);
  }

  public void writeValues() {
    Configuration instance = Configuration.getInstance();

    instance.setDefaultDataFolder(dataFolder.getText());
    instance.setDefaultRespFolder(respFolder.getText());
    instance.setDefaultOutputFolder(outputFolder.getText());

    instance.setFDSNProtocol(fdsnProtocol.getText());
    instance.setFDSNDomain(fdsnDomain.getText());
    instance.setFDSNPath(fdsnService.getText());
    instance.setFDSNPort(Integer.parseInt(fdsnPort.getText()));

    instance.setLineWidthOffset(Integer.parseInt(lineWidthOffset.getText()));
    instance.setUseColorblindColors(colorblindColors.isSelected());

    instance.saveCurrentConfig();
  }


  private class PositiveIntegerFilter extends DocumentFilter {

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
