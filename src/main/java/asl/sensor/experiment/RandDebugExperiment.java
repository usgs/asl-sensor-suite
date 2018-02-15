package asl.sensor.experiment;

import asl.sensor.input.DataStore;

/**
 * Debugging class; outputs inputted and estimated response curves.
 * see the RandomizedExperiment class.
 * @author akearns
 *
 */
public class RandDebugExperiment extends RandomizedExperiment {

  public final boolean SKIP_SOLVING = true;
  
  @Override
  protected void backend(DataStore ds) {
    super.backend(ds);
  }
  
  @Override
  public boolean getSolverState() {
    return SKIP_SOLVING;
  }

}
