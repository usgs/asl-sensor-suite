package asl.sensor;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.List;

import org.apache.commons.math3.complex.Complex;
import org.jfree.data.xy.XYSeriesCollection;

import asl.sensor.experiment.RandomizedExperiment;
import asl.sensor.input.DataBlock;
import asl.sensor.input.DataStore;
import asl.sensor.input.InstrumentResponse;
import asl.sensor.utils.TimeSeriesUtils;
import py4j.GatewayServer;

public class RandomCalShell {

  public RandomCalShell() {
    
  }

  /**
   * Acquire data and run calibration over it.
   * Returns the experiment (all data kept locally to maintain thread safety)
   * @param calFileName Filename of calibration signal
   * @param outFileName Filename of sensor output 
   * @param respName Filename of response to load in
   * @param respEmbd True if response is an embedded response in program
   * @param startTime Long representing ms-since-epoch of data start time
   * @param endTime Long representing ms-since-epoch of data end time
   * @param lowFreq True if a low-freq cal should be run
   * @return The experiment the data was run on (for getting results)
   */
  public RandomizedExperiment populateDataAndRun(String calFileName, 
      String outFileName, String respName, boolean respEmbd, long startTime, 
      long endTime, boolean lowFreq) {

    RandomizedExperiment re = new RandomizedExperiment();
    
    try {
      DataStore ds = new DataStore();
      String calFilt = TimeSeriesUtils.getMplexNameList(calFileName).get(0);
      DataBlock calBlock = TimeSeriesUtils.getTimeSeries(calFileName, calFilt);
      String outFilt = TimeSeriesUtils.getMplexNameList(outFileName).get(0);
      DataBlock outBlock = TimeSeriesUtils.getTimeSeries(outFileName, outFilt);
      InstrumentResponse ir;
      if (respEmbd) {
        ir = InstrumentResponse.loadEmbeddedResponse(respName);
      } else{
        ir = new InstrumentResponse(respName);
      }

      ds.setBlock(0, calBlock);
      ds.setBlock(1, outBlock);
      ds.setResponse(1, ir);
      ds.trim(startTime, endTime);

      re.setLowFreq(lowFreq);
      re.runExperimentOnData(ds);
    } catch (FileNotFoundException e) {
      e.printStackTrace();
    } catch (IOException e) {
      e.printStackTrace();
    }

    return re;
    
  }
  
  public static void main(String[] args) {
    GatewayServer gatewayServer = new GatewayServer( new RandomCalShell() );
    gatewayServer.start();
    System.out.println("Gateway Server Started");
  }

}
