package asl.sensor.input;

import org.apache.commons.configuration.ConfigurationException;
import org.apache.commons.configuration.XMLConfiguration;

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

  private Configuration() {
    try {
      XMLConfiguration config = new XMLConfiguration(CONFIG_PATH);
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

      useColorblindColors =
          config.getBoolean("VisualOptions.ColorblindFriendly", true);
      lineWidthOffset =
          config.getInt("VisualOptions.LineThicknessIncrease", 2);

    } catch (ConfigurationException e) {
      System.out.println("Error in processing configuration file! Using default settings.");
      System.out.println("Check that config.xml exists in the program directory.");
    }
  }

  synchronized public static Configuration getInstance() {
    if (instance == null) {
        instance = new Configuration();
    }
    return instance;
  }

  public String getDefaultDataFolder() {
    return defaultDataFolder;
  }

  public String getDefaultRespFolder() {
    return defaultRespFolder;
  }

  public String getDefaultOutputFolder() {
    return defaultOutputFolder;
  }

  public boolean useColorblindColors() {
    return useColorblindColors;
  }

  public int getLineWidthOffset() {
    return lineWidthOffset;
  }

  public String getFDSNProtocol() {
    return fdsnProtocol;
  }

  public String getFDSNDomain() {
    return fdsnDomain;
  }

  public String getFDSNPath() {
    return fdsnService;
  }

}
