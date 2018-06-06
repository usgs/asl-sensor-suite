package asl.sensor.output;

import java.util.HashMap;
import java.util.Map;

/**
 * Abstract class that forms an easy interface by which external programs can interact with
 * all the data provided from a calibration solver experiment. CalResult contains two maps, one of
 * which is a map from string descriptors to a set of binary objects representing images as PNGs
 * (these are stored as byte arrays to be more easily imported into, say, a Django database backend)
 * and the other of which is a map from string descriptors to the variables fit by the solver,
 * given as a list of doubles (which has more than one entry in the case of, say, poles and zeros
 * returned by a randomized cal experiment).
 * Implementing classes don't need to add additional functions but must populate the maps with
 * actual data to be returned, which varies in content depending on the type of calibration done.
 * This class is not useful for the GUI interface, as the results there are contained within the
 * panel and saved to PDF reports as desired. In practice this is (currently) only used by the
 * ASL calibration tracking database.
 */
public abstract class CalResult {

  Map<String, double[]> numerMap;
  Map<String, byte[]> imageMap;

  CalResult() {
    numerMap = new HashMap<>();
    imageMap = new HashMap<>();
  }

  /**
   * Return the map of images
   * @return
   */
  public Map<String, byte[]> getImageMap() {
    return imageMap;
  }

  public Map<String, double[]> getNumerMap() {
    return numerMap;
  }
}
