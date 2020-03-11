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
 * the path of the FDSN data and metadata acquisition services,
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

  private String fdsnDataProtocol = "http";
  private String fdsnDataDomain = "service.iris.edu";
  private String fdsnDataService = "fdsnws";
  private int fdsnDataPort = -1; // determined based on given protocol
    // will default to 80 if other terms are not set from these default values

  private String fdsnMetaProtocol = "http";
  private String fdsnMetaDomain = "service.iris.edu";
  private String fdsnMetaService = "fdsnws";
  private int fdsnMetaPort = -1;

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

      String fdsnProtocolParam = config.getString("FDSNData.Protocol");
      if (fdsnProtocolParam != null) {
        fdsnDataProtocol = fdsnProtocolParam;
      }
      String fdsnDomainParam = config.getString("FDSNData.Domain");
      if (fdsnDomainParam != null) {
        fdsnDataDomain = fdsnDomainParam;
      }
      String fdsnServiceParam = config.getString("FDSNData.Service");
      if (fdsnServiceParam != null) {
        fdsnDataService = fdsnServiceParam;
      }
      int fdsnPortParam = config.getInt("FDSNData.Port", -1);
      if (fdsnPortParam < 0) {
        // use value of 443 if https, use value of 80 if http
        fdsnPortParam = fdsnDataProtocol.equals("https") ? 443 : 80;
      }
      fdsnDataPort = fdsnPortParam;

      fdsnProtocolParam = config.getString("FDSNMeta.Protocol");
      if (fdsnProtocolParam != null) {
        fdsnMetaProtocol = fdsnProtocolParam;
      }
      fdsnDomainParam = config.getString("FDSNMeta.Domain");
      if (fdsnDomainParam != null) {
        fdsnMetaDomain = fdsnDomainParam;
      }
      fdsnServiceParam = config.getString("FDSNMeta.Service");
      if (fdsnServiceParam != null) {
        fdsnMetaService = fdsnServiceParam;
      }
      fdsnPortParam = config.getInt("FDSNMeta.Port", -1);
      if (fdsnPortParam < 0) {
        // use value of 443 if https, use value of 80 if http
        fdsnPortParam = fdsnMetaProtocol.equals("https") ? 443 : 80;
      }
      fdsnMetaPort = fdsnPortParam;

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

  /**
   * Set the starting directory for SEED file loading.
   * Used when selecting a new folder to load data from but only written
   * out to file when the configuration is edited from the main window.
   * @param replacement Path (relative or absolute) to the folder where data is kept.
   */
  public void setDefaultDataFolder(String replacement) {
    File config = new File(replacement);
    if (config.exists() && config.canRead()) {
      defaultDataFolder = replacement;
    }
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
   * Set the starting directory for resp loading.
   * Used when selecting a new folder to load data from but only written
   * out to file when the configuration is edited from the main window.
   * @param replacement Path (relative or absolute) to the folder where resps are kept.
   */
  public void setDefaultRespFolder(String replacement) {
    File config = new File(replacement);
    if (config.exists() && config.canRead()) {
      defaultRespFolder = replacement;
    }
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
   * Set the default report folder.
   * @param replacement The folder to write reports to.
   */
  public void setDefaultOutputFolder(String replacement) {
    File config = new File(replacement);
    if (config.exists() && config.canWrite()) {
      defaultOutputFolder = replacement;
    }
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
   * Used to toggle whether or not plots should use colorblind-friendly colors.
   * When true, red and green are replaced with teal and orange.
   * @param trueIfUsed True if colorblind-friendly colors should be used
   */
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

  /**
   * Set the thickness to apply to output plots.
   * This will not affect plots that have already been produced by the program.
   * @param replacementOffset New thickness increase to apply to plots.
   */
  public void setLineWidthOffset(int replacementOffset) {
    lineWidthOffset = replacementOffset;
  }

  /**
   * Gets the current protocol to use for FDSN data acquisition. This is a string value.
   * If not set in the configuration file, it defaults to "http".
   *
   * The property is defined from Configuration.FDSNData.Protocol
   * @return The FDSN connection protocol (URL prefix).
   */
  public String getFDSNDataProtocol() {
    return fdsnDataProtocol;
  }

  /**
   * Set the connection protocol for FDSN data services.
   * @param replacement Protocol to use (i.e., http, https)
   */
  public void setFDSNDataProtocol(String replacement) {
    fdsnDataProtocol = replacement;
  }

  /**
   * Gets the current domain to use for FDSN data acquisition. This is a string value.
   * If not set in the configuration file it defaults to "service.iris.edu"
   *
   * The property is defined from Configuration.FDSNData.Domain
   * @return The FDSN webservice domain (URL).
   */
  public String getFDSNDataDomain() {
    return fdsnDataDomain;
  }

  /**
   * Set the FDSN data services domain url.
   * @param replacement New URL to use webservices from
   */
  public void setFDSNDataDomain(String replacement) {
    fdsnDataDomain = replacement;
  }

  /**
   * Gets the current path (service) to use for FDSN data acquisition. This is a string value.
   * If not set in the configuration file it defaults to "fdsnws".
   *
   * The property is defined from Configuration.FDSNData.Service
   * @return The FDSN webservice path.
   */
  public String getFDSNDataPath() {
    return fdsnDataService;
  }

  /**
   * Set the path (service) to use for FDSN data acquisition
   * @param replacement The FDSN webservice path to use
   */
  public void setFDSNDataPath(String replacement) {
    fdsnDataService = replacement;
  }

  /**
   * Gets the port to use for FDSN data acquisition. This is an integer value.
   * If not set in the configuration file it defaults to 80.
   *
   * The property is defined from Configuration.FDSNData.Port
   * @return The FDSN webservice connection port.
   */
  public int getFDSNDataPort() {
    return fdsnDataPort;
  }

  /**
   * Set the port to use for FDSN data acquisition.
   * @param replacementPort New port to use for FDSN Meta.
   */
  public void setFDSNDataPort(int replacementPort) {
    fdsnDataPort = replacementPort;
  }

  /**
   * Gets the current protocol to use for FDSN metadata acquisition. This is a string value.
   * If not set in the configuration file, it defaults to "http".
   *
   * The property is defined from Configuration.FDSNMeta.Protocol
   * @return The FDSN connection protocol (URL prefix).
   */
  public String getFDSNMetaProtocol() {
    return fdsnMetaProtocol;
  }

  /**
   * Set the connection protocol for FDSN metadata services.
   * @param replacement Protocol to use (i.e., http, https)
   */
  public void setFDSNMetaProtocol(String replacement) {
    fdsnMetaProtocol = replacement;
  }

  /**
   * Gets the current domain to use for FDSN metadata acquisition. This is a string value.
   * If not set in the configuration file it defaults to "service.iris.edu"
   *
   * The property is defined from Configuration.FDSNMeta.Domain
   * @return The FDSN webservice domain (URL).
   */
  public String getFDSNMetaDomain() {
    return fdsnMetaDomain;
  }

  /**
   * Set the FDSN metadata services domain url.
   * @param replacement New URL to use webservices from
   */
  public void setFDSNMetaDomain(String replacement) {
    fdsnMetaDomain = replacement;
  }

  /**
   * Gets the current path (service) to use for FDSN metadata acquisition. This is a string value.
   * If not set in the configuration file it defaults to "fdsnws".
   *
   * The property is defined from Configuration.FDSNMeta.Service
   * @return The FDSN webservice path.
   */
  public String getFDSNMetaPath() {
    return fdsnMetaService;
  }

  /**
   * Set the path (service) to use for FDSN metadata acquisition
   * @param replacement The FDSN webservice path to use
   */
  public void setFDSNMetaPath(String replacement) {
    fdsnMetaService = replacement;
  }

  /**
   * Gets the port to use for FDSN metadata acquisition. This is an integer value.
   * If not set in the configuration file it defaults to 80.
   *
   * The property is defined from Configuration.FDSNMeta.Port
   * @return The FDSN webservice connection port.
   */
  public int getFDSNMetaPort() {
    return fdsnMetaPort;
  }

  /**
   * Set the port to use for FDSN metadata acquisition.
   * @param replacementPort New port to use for FDSN Meta.
   */
  public void setFDSNMetaPort(int replacementPort) {
    fdsnMetaPort = replacementPort;
  }

  /**
   * Writes out the current configuration to file.
   * This is called when the configuration is saved via the GUI utils.
   */
  public void saveCurrentConfig() {
    XMLConfiguration config;
    try {
      config = new XMLConfiguration(loadedConfigPath);
      config.load();

      config.setProperty("LocalPaths.DataPath", defaultDataFolder);
      config.setProperty("LocalPaths.RespPath", defaultRespFolder);
      config.setProperty("LocalPaths.ReportPath", defaultOutputFolder);
      config.setProperty("FDSNData.Protocol", fdsnDataProtocol);
      config.setProperty("FDSNData.Domain", fdsnDataDomain);
      config.setProperty("FDSNData.Service", fdsnDataService);
      config.setProperty("FDSNData.Port", fdsnDataPort);
      config.setProperty("FDSNMeta.Protocol", fdsnMetaProtocol);
      config.setProperty("FDSNMeta.Domain", fdsnMetaDomain);
      config.setProperty("FDSNMeta.Service", fdsnMetaService);
      config.setProperty("FDSNMeta.Port", fdsnMetaPort);
      config.setProperty("VisualOptions.ColorblindFriendly", useColorblindColors);
      config.setProperty("VisualOptions.LineThicknessIncrease", lineWidthOffset);

      config.save();
    } catch (ConfigurationException e) {
      logger.error(e);
    }
  }

}
