package asl.sensor.experiment;

import asl.sensor.input.DataStore;
import asl.sensor.input.InstrumentResponse;
import asl.sensor.utils.NumericUtils;
import java.time.ZoneOffset;
import java.time.format.DateTimeFormatter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;
import org.apache.commons.math3.complex.Complex;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

/**
 * Produces plots of response curves' magnitudes (Bode plot) and angle of
 * rotation in complex space. 1-3 responses can be plotted at a time.
 * No timeseries data (that is, miniSEED) is used in this calculation.
 * Response curves can be plotted in either frequency or interval space
 * (units of Hz or seconds respectively).
 *
 * @author akearns
 */
public class ResponseExperiment extends Experiment {

  public static final String MAGNITUDE = "Response amplitude";
  public static final String ARGUMENT = "Response phase";

  private boolean freqSpace; // choose between units of Hz or seconds (time between samples)

  private Set<InstrumentResponse> responses;

  public ResponseExperiment() {
    super();
    freqSpace = false;
  }

  @Override
  protected void backend(DataStore dataStore) {

    responses = new HashSet<>();

    double lowFreq = .0001;
    double highFreq = 200;

    int pointCount = 100000;
    double linearChange = (highFreq - lowFreq) / pointCount;
    // find logarithmic parameters for linear components
    double b = Math.log10(lowFreq / highFreq) / (lowFreq - highFreq);
    double a = lowFreq / Math.pow(10, b * lowFreq);

    // hard-code length here because the limits of the calculated range are fixed
    double[] freqArray = new double[pointCount];

    double currentFreq = lowFreq;
    for (int i = 0; i < freqArray.length; ++i) {
      freqArray[i] = currentFreq;
      currentFreq = a * Math.pow(10, b * (i * linearChange));
    }

    // used to prevent issues with duplicate response plotting / XYSeries names
    Set<String> respNames = new HashSet<>();

    XYSeriesCollection args = new XYSeriesCollection();
    XYSeriesCollection mags = new XYSeriesCollection();

    for (int responseIndex = 0; responseIndex < 3; ++responseIndex) {
      if (!dataStore.responseIsSet(responseIndex)) {
        continue;
      }

      InstrumentResponse instrumentResponse = dataStore.getResponse(responseIndex);
      String name = instrumentResponse.getName() + " [" +
          DateTimeFormatter.ofPattern("uuuu.DDD").withZone(ZoneOffset.UTC)
              .format(instrumentResponse.getEpochStart()) + ']';
      if (respNames.contains(name)) {
        continue;
      } else {
        respNames.add(name);
        responses.add(instrumentResponse);
      }

      Complex[] result = instrumentResponse.applyResponseToInput(freqArray);

      double phiPrev = 0; // use with unwrapping
      XYSeries magnitude = new XYSeries(name + " " + MAGNITUDE);
      XYSeries argument = new XYSeries(name + " " + ARGUMENT);
      for (int i = 0; i < freqArray.length; ++i) {
        Complex tmp = result[i];
        double phi = NumericUtils.atanc(tmp);
        phi = NumericUtils.unwrap(phi, phiPrev);
        phiPrev = phi;
        phi = Math.toDegrees(phi);
        double magAccel = 10 * Math.log10(tmp.abs());
        double xVal = 1 / freqArray[i];
        if (freqSpace) {
          xVal = freqArray[i];
        }
        magnitude.add(xVal, magAccel);
        argument.add(xVal, phi);
      }

      mags.addSeries(magnitude);
      args.addSeries(argument);

    }

    xySeriesData.add(mags);
    xySeriesData.add(args);

    dataNames = new ArrayList<>(respNames);

  }

  @Override
  public int blocksNeeded() {
    return 0;
  }

  @Override
  public long getEnd() {
    return 0L;
  }

  public InstrumentResponse[] getResponses() {
    return responses.toArray(new InstrumentResponse[0]);
  }

  @Override
  public long getStart() {
    // there is no actual timeseries data used in this experiment
    // only response data, so start and end times are not defined
    return 0L;
  }

  @Override
  public boolean hasEnoughData(DataStore dataStore) {
    for (int i = 0; i < 3; ++i) {
      if (dataStore.responseIsSet(i)) {
        return true;
      }
    }

    return false;
  }

  /**
   * Used to set the x-axis over which the response curve is plotted,
   * either frequency (Hz) units or sample-interval (s) units
   *
   * @param freqSpace True if the plot should use units of Hz
   */
  public void setFreqSpace(boolean freqSpace) {
    this.freqSpace = freqSpace;
  }

}
