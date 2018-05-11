package asl.sensor.experiment;

import asl.sensor.input.DataStore;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

public class MockExperiment extends Experiment {
  boolean backendCalled = false;
  boolean hasEnoughDataCalled = false;
  boolean blocksNeededCalled = false;
  /**
   * Counts the total number of times fireStateChange was called.
   * It is used for verifing that methods are still calling when we expect them too.
   */
  int numberOfChangesFired = 0;

  MockExperiment() {
    super();
    this.addChangeListener(new changeCountingListener());
  }

  @Override
  protected void backend(final DataStore dataStore) {
    backendCalled = true;
  }

  @Override
  public int blocksNeeded() {
    blocksNeededCalled = true;
    return 3;
  }

  @Override
  public boolean hasEnoughData(DataStore dataStore) {
    hasEnoughDataCalled = true;
    return true;
  }

  private class changeCountingListener implements ChangeListener {
    public void stateChanged(ChangeEvent event) {
      numberOfChangesFired++;
    }
  }

}
