package asl.sensor.experiment;

import asl.sensor.input.DataStore;

/**
 * Mock class that is used to represent experiments that handle timeseries data for testing
 * output formatting routines
 */
public class MockExperimentThatNeedsBlocks extends MockExperiment{

  @Override
  public int blocksNeeded() {
    return 1;
  }

  @Override
  public boolean hasEnoughData(DataStore dataStore) {
    return dataStore.areAnyBlocksSet();
  }

}
