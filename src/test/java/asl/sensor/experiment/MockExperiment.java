package asl.sensor.experiment;

import asl.sensor.input.DataStore;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

class MockExperiment extends Experiment {

  boolean backendCalled = false;
  boolean hasEnoughDataCalled = false;
  boolean blocksNeededCalled = false;

  boolean setHasEnoughData = true;
  private int numberBlocksNeeded = 0;
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
    return numberBlocksNeeded;
  }

  @Override
  public boolean hasEnoughData(DataStore dataStore) {
    hasEnoughDataCalled = true;
    return setHasEnoughData;
  }

  private class changeCountingListener implements ChangeListener {

    public void stateChanged(ChangeEvent event) {
      numberOfChangesFired++;
    }
  }

}
