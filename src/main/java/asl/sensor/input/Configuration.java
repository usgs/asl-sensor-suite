package asl.sensor.input;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import org.apache.commons.configuration.ConfigurationException;
import org.apache.commons.configuration.XMLConfiguration;
import org.apache.log4j.Logger;

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

  private static final String DEFAULT_CONFIG_PATH = "sensor-suite-config.xml";
  private static final Logger logger = Logger.getLogger(Configuration.class);

  private String loadedConfigPath = DEFAULT_CONFIG_PATH;

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
    logger.info("Attempting reading in config file from " + configLocation);
    try {
      XMLConfiguration config = new XMLConfiguration(configLocation);
      config.load();
      config.setAutoSave(true);

      String defaultDataFolderParam = config.getString("LocalPaths.DataPath");
      if (defaultDataFolderParam != null) {
        defaultDataFolder = defaultDataFolderParam;
      }
      String defaultRespFolderParam = config.getString("LocalPaths.RespPath");
      if (defaultRespFolderParam != null) {
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
      fdsnPort = config.getInt("FDSNPaths.Port", 80);


      useColorblindColors =
          config.getBoolean("VisualOptions.ColorblindFriendly", true);
      lineWidthOffset =
          config.getInt("VisualOptions.LineThicknessIncrease", 2);

      try {
        loadedConfigPath = config.getFile().getCanonicalPath();
        logger.info("Succesfully loaded in configuration: " + loadedConfigPath);
      } catch (IOException e) {
        e.printStackTrace();
      }
    } catch (ConfigurationException e) {
      logger.error("Error encountered while reading XML file, load failed, using defaults", e);
    }
  }

  private static boolean copyEmbedXML(String pathToPlaceFile){
    File fileOut = new File(pathToPlaceFile);
    try (InputStream stream =
        Configuration.class.getClassLoader().getResourceAsStream(DEFAULT_CONFIG_PATH)) {
      try {
        logger.info("Copying over embedded jar file to absolute path " + fileOut.getAbsolutePath());
        Files.copy(stream, fileOut.toPath());
        return true;
      } catch (IOException e) {
        logger.warn("Could not copy over the file...", e);
      }
    } catch (IOException e) {
      logger.error("Major error: config XML file not part of resources!!");
    }
    return false;
  }

  /**
   * Gets the current instance of the configuration, or creates one if none exists
   * @return the current configuration instance
   */
  synchronized public static Configuration getInstance() {
    return getInstance(System.getProperty("user.dir") + File.separator + DEFAULT_CONFIG_PATH);
  }

  /**
   * Gets the current instance of the configuration, or creates one from a specified file if none
   * exists.
   * @param configLocation New configuration file location to read from
   * @return the current configuration instance
   */
  synchronized public static Configuration getInstance(String configLocation) {
    if (instance == null) {
      File config = new File(configLocation);
      if (!config.exists()) {
        boolean success = copyEmbedXML(configLocation);
        if (!success) {
          logger.warn("Could not find or write to specified config location: " + configLocation);
          logger.warn("Will attempt to initialize config file at user home directory.");
          configLocation = System.getProperty("user.home") + File.separator + DEFAULT_CONFIG_PATH;
          config = new File(configLocation);
          if (!config.exists()) {
            success = copyEmbedXML(configLocation);
            if (!success) {
              logger.warn("Could not find or write to user home directory either!");
            }
          }
        }
      }
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

  public void setDefaultDataFolder(String replacement) {
    defaultDataFolder = replacement;
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

  public void setDefaultRespFolder(String replacement) {
    defaultRespFolder = replacement;
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

  public void setDefaultOutputFolder(String replacement) {
    defaultOutputFolder = replacement;
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

  public void setUseColorblindColors(boolean trueIfUsed) {
    useColorblindColors = trueIfUsed;
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

  public void setLineWidthOffset(int replacementOffset) {
    lineWidthOffset = replacementOffset;
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

  public void setFDSNProtocol(String replacement) {
    fdsnProtocol = replacement;
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

  public void setFDSNDomain(String replacement) {
    fdsnDomain = replacement;
  }

  /**
   * Gets the current path (service) to use for FDSN data acquisition. This is a string value.
   * If not set in the configuration file it defaults to "fdsnws".
   *
   * The property is defined from Configuration.FDSNPaths.Service
   * @return The FDSN webservice path.
   */
  public String getFDSNPath() {
    return fdsnService;
  }

  public void setFDSNPath(String replacement) {
    fdsnService = replacement;
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

  public void setFDSNPort(int replacementPort) {
    fdsnPort = replacementPort;
  }

  public void saveCurrentConfig() {
    XMLConfiguration config;
    try {
      config = new XMLConfiguration(loadedConfigPath);
      config.load();

      config.setProperty("LocalPaths.DataPath", defaultDataFolder);
      config.setProperty("LocalPaths.RespPath", defaultRespFolder);
      config.setProperty("LocalPaths.ReportPath", defaultOutputFolder);
      config.setProperty("FDSNPaths.Protocol", fdsnProtocol);
      config.setProperty("FDSNPaths.Domain", fdsnDomain);
      config.setProperty("FDSNPaths.Service", fdsnService);
      config.setProperty("FDSNPaths.Port", fdsnPort);
      config.setProperty("VisualOptions.ColorblindFriendly", useColorblindColors);
      config.setProperty("VisualOptions.LineThicknessIncrease", lineWidthOffset);

      config.save();
    } catch (ConfigurationException e) {
      logger.error(e);
    }
  }

}
