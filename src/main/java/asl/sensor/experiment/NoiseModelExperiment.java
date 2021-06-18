package asl.sensor.experiment;

import asl.sensor.input.DataStore;
import java.util.ArrayList;

public class NoiseModelExperiment extends SpectralAnalysisExperiment {

  public NoiseModelExperiment() {
    xySeriesData = new ArrayList<>();
  }

  @Override
  protected void backend(DataStore dataStore) {
    // doesn't do anything, nor does it need to be called
  }

  @Override
  public int blocksNeeded() {
    return 0;
  }

  @Override
  public boolean hasEnoughData(DataStore dataStore) {
    return false;
  }
}
