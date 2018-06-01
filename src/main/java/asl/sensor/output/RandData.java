package asl.sensor.output;

import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.TimeZone;

/**
 * Output format for randomized calibration solver as used in ASL's automated test tracking database
 * This class retains more granular access of its fields besides the map for compatibility with
 * older versions of the code in that database.
 */
public class RandData extends CalResult {

  private String[] gapNameIdentifiers;
  private Date[][] gapStarts;
  private Date[][] gapEnds;

  public RandData(double[] fitPoles, double[] fitZeros, double[] initialPoles,
      double[] initialZeros,
      byte[][] images,
      String[] gapNames, Date[][] gapStartTimes, Date[][] gapEndTimes) {
    super();
    numerMap.put("Best_fit_poles", fitPoles);
    numerMap.put("Best_fit_zeros", fitZeros);
    numerMap.put("Initial_poles", initialPoles);
    numerMap.put("Initial_zeros", initialZeros);
    imageMap.put("Response_amplitudes", images[0]);
    imageMap.put("Response_phases", images[1]);
    imageMap.put("Amplitude_error", images[2]);
    imageMap.put("Phase_error", images[3]);
    gapNameIdentifiers = gapNames;
    gapStarts = gapStartTimes;
    gapEnds = gapEndTimes;
  }

  public byte[] getAmpErrorImage() {
    return imageMap.get("Amplitude_error");
  }

  public byte[] getAmpImage() {
    return imageMap.get("Response_amplitudes");
  }

  public double[] getFitPoles() {
    return numerMap.get("Best_fit_poles");
  }

  public double[] getFitZeros() {
    return numerMap.get("Best_fit_zeros");
  }

  public Date[][] getGapEndDates() {
    return gapEnds;
  }

  public String[] getGapIdentifiers() {
    return gapNameIdentifiers;
  }

  public String getGapInfoAsString() {
    SimpleDateFormat sdf = new SimpleDateFormat("DD.HH:m:s");
    sdf.setTimeZone(TimeZone.getTimeZone("UTC"));
    return getGapInfoAsString(sdf);
  }

  public String getGapInfoAsString(DateFormat df) {
    StringBuilder sb = new StringBuilder();
    for (int j = 0; j < gapNameIdentifiers.length; ++j) {
      sb.append(gapNameIdentifiers[j]);
      sb.append(":\n");
      for (int i = 0; i < gapStarts[j].length; ++i) {
        sb.append("\t");
        Date start = gapStarts[j][i];
        Date end = gapEnds[j][i];
        sb.append(df.format(start));
        sb.append("\t");
        sb.append(df.format(end));
        sb.append("\n");
      }
      sb.append("\n");
    }
    return sb.toString();
  }

  public Date[][] getGapStartDates() {
    return gapStarts;
  }

  public double[] getInitPoles() {
    return numerMap.get("Initial_poles");
  }

  public double[] getInitZeros() {
    return numerMap.get("Initial_zeros");
  }

  public byte[] getPhaseErrorImage() {
    return imageMap.get("Phase_error");
  }

  public byte[] getPhaseImage() {
    return imageMap.get("Response_phases");
  }
}
