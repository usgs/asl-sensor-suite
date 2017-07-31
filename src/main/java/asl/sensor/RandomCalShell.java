package asl.sensor;

import java.io.FileNotFoundException;

import asl.sensor.experiment.RandomizedExperiment;
import asl.sensor.input.DataBlock;
import asl.sensor.input.DataStore;
import asl.sensor.input.InstrumentResponse;
import asl.sensor.utils.TimeSeriesUtils;

public class RandomCalShell {

  /**
   * Shell interface for sensor suite, runnable from command line
   * @param args
   */
  public static void main(String[] args) {
    // TODO parse arguments and set up data structures to run operations
    
    // command line arg structure, by index:
    // 0. denotes whether or not calibration is high or low period
    // 1. location of file with cal input data
    // 2. location of file with sensor output
    // 3. location of RESP of sensor
    // 4. denote whether or not calibration is high or low period
    // 5. start time to trim data down to (assume ddd.hh:mm:ss.ms format)
    // 6. end time to trim data down to (same format)
    
    RandomizedExperiment re = new RandomizedExperiment();
    if ( args[0].toUpperCase().startsWith("H") ) {
      re.setLowFreq(false);
    } else {
      re.setLowFreq(true);
    }
    
    try {
      DataStore ds = new DataStore();
      String calFilt = TimeSeriesUtils.getMplexNameList(args[1]).get(0);
      DataBlock calBlock = TimeSeriesUtils.getTimeSeries(args[1], calFilt);
      String outFilt = TimeSeriesUtils.getMplexNameList(args[2]).get(0);
      DataBlock outBlock = TimeSeriesUtils.getTimeSeries(args[1], calFilt);
      InstrumentResponse ir;
      
    } catch (FileNotFoundException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
    
    // what sort of autocorrection
  }

}
