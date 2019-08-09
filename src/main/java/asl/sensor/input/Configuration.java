package asl.sensor.input;

import org.apache.commons.configuration.ConfigurationException;
import org.apache.commons.configuration.XMLConfiguration;

/**
 * Configuration file including parameters that have been requested for features by program users.
 * These include the thickness of lines being plotted for experiment results,
 * the path of the FDSN data acquisition service,
 * the default location from which to load data and responses,
 * the default folder to which reports are outputted (i.e,. as PDF files),
 * and whether or not to use colorblind-friendly colors in result plots.
 */
public class Configuration {

  private static Configuration instance;

  private static final String CONFIG_PATH = "config.xml";

  private String defaultDataFolder = "data";
  private String defaultRespFolder = "responses";
  private String defaultOutputFolder = System.getProperty("user.home");

  private int lineWidthOffset = 2;
  private boolean useColorblindColors = true;

  private String fdsnProtocol = "http";
  private String fdsnDomain = "service.iris.edu";
  private String fdsnService = "fdsnws";
  private int fdsnPort = 80;

  private Configuration(String configLocation) {
    try {
      XMLConfiguration config = new XMLConfiguration(configLocation);
      config.load();

      String defaultDataFolderParam = config.getString("LocalPaths.DataPath");
      if (defaultDataFolderParam != null) {
        defaultDataFolder = defaultDataFolderParam;
      }
      String defaultRespFolderParam = config.getString("LocalPaths.RespPath");
      if (defaultDataFolderParam != null) {
        defaultRespFolder = defaultRespFolderParam;
      }
      String defaultOutputFolderParam =
          config.getString("LocalPaths.ReportPath");
      if (defaultOutputFolderParam != null) {
        defaultOutputFolder = defaultOutputFolderParam;
      }

      String fdsnProtocolParam = config.getString("FDSNPaths.Protocol");
      if (fdsnProtocolParam != null) {
        fdsnProtocol = fdsnProtocolParam;
      }
      String fdsnDomainParam = config.getString("FDSNPaths.Domain");
      if (fdsnDomainParam != null) {
        fdsnDomain = fdsnDomainParam;
      }
      String fdsnServiceParam = config.getString("FDSNPaths.Service");
      if (fdsnServiceParam != null) {
        fdsnService = fdsnServiceParam;
      }
      int fdsnPort = config.getInt("FDSNPaths.Port", 80);


      useColorblindColors =
          config.getBoolean("VisualOptions.ColorblindFriendly", true);
      lineWidthOffset =
          config.getInt("VisualOptions.LineThicknessIncrease", 2);

    } catch (ConfigurationException e) {
      System.out.println("Error in processing configuration file! Using default settings.");
      System.out.println("Check that config.xml exists in the program directory.");
    }
  }

  /**
   * Gets the current instance of the configuration, or creates one if none exists
   * @return the current configuration instance
   */
  synchronized public static Configuration getInstance() {
    if (instance == null) {
        instance = new Configuration(CONFIG_PATH);
    }
    return instance;
  }

  /**
   * Gets the current instance of the configuration, or creates one from a specified file if none
   * exists.
   * @param configLocation New configuration file location to read from
   * @return the current configuration instance
   */
  synchronized public static Configuration getInstance(String configLocation) {
    if (instance == null) {
      instance = new Configuration(configLocation);
    }
    return instance;
  }

  /**
   * Gets the default data folder. If none was specified in the XML file, the default is the 'data'
   * subdirectory under the current working directory.
   *
   * The property is defined from Configuration.LocalPaths.DataPath
   * @return The folder to start looking for seed files from.
   */
  public String getDefaultDataFolder() {
    return defaultDataFolder;
  }

  /**
   * Gets the default resp folder. If none was specified in the XML file, the default is the 'resps'
   * subdirectory under the current working directory.
   *
   * The property is defined from Configuration.LocalPaths.RespPath
   * @return The folder to start looking for resp files from.
   */
  public String getDefaultRespFolder() {
    return defaultRespFolder;
  }

  /**
   * Gets the default data output folder. If none was specified in the XML file, the default is the
   * user's home directory.
   *
   * The property is defined from Configuration.LocalPaths.ReportPath
   * @return The folder where experiment reports are placed by default.
   */
  public String getDefaultOutputFolder() {
    return defaultOutputFolder;
  }

  /**
   * Gets the current choice of whether or not experiment plots should use colorblind-friendly
   * colors to display data. If true, colorblind-friendly colors are chosen (teal and orange)
   * to replace a default red and green. True is the default value if not specified.
   *
   * The property is defined from Configuration.VisualOptions.ColorblindFriendly as a boolean
   * @return Whether or not to use colorblind-friendly colors.
   */
  public boolean useColorblindColors() {
    return useColorblindColors;
  }

  /**
   * Gets the current line width offset (i.e., value to add to line thickness in experiment plots).
   * This value is set to default to 2 if not specified in the configuration file.
   *
   * The property is defined from Configuration.VisualOptions.LineThicknessIncrease as an integer
   * @return The thickness increase to apply to plots produced as experiment results.
   */
  public int getLineWidthOffset() {
    return lineWidthOffset;
  }

  /**
   * Gets the current protocol to use for FDSN data acquisition. This is a string value.
   * If not set in the configuration file, it defaults to "http".
   *
   * The property is defined from Configuration.FDSNPaths.Protocol
   * @return The FDSN connection protocol (URL prefix).
   */
  public String getFDSNProtocol() {
    return fdsnProtocol;
  }

  /**
   * Gets the current domain to use for FDSN data acquisition. This is a string value.
   * If not set in the configuration file it defaults to "service.iris.edu"
   *
   * The property is defined from Configuration.FDSNPaths.Domain
   * @return The FDSN webservice domain (URL).
   */
  public String getFDSNDomain() {
    return fdsnDomain;
  }

  /**
   * Gets the current path to use for FDSN data acquisition. This is a string value.
   * If not set in the configuration file it defaults to "fdsnws".
   *
   * The property is defined from Configuration.FDSNPaths.Service
   * @return The FDSN webservice path.
   */
  public String getFDSNPath() {
    return fdsnService;
  }

  /**
   * Gets the port to use for FDSN data acquisition. This is an integer value.
   * If not set in the configuration file it defaults to 80.
   *
   * The property is defined from Configuration.FDSNPaths.Port
   * @return The FDSN webservice connection port.
   */
  public int getFDSNPort() {
    return fdsnPort;
  }

}
