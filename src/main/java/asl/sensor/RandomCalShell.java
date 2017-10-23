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

  private RandomizedExperiment re;
  private DataStore ds;

  public RandomCalShell() {
    re = new RandomizedExperiment();
    ds = new DataStore();
  }

  /**
   * Acquire data and run calibration over it
   * @param calFileName Filename of calibration signal
   * @param outFileName Filename of sensor output 
   * @param respName Filename of response to load in
   * @param respEmbd True if response is an embedded response in program
   * @param startTime Long representing ms-since-epoch of data start time
   * @param endTime Long representing ms-since-epoch of data end time
   * @param lowFreq True if a low-freq cal should be run
   */
  public void populateDataAndRun(String calFileName, String outFileName,
      String respName, boolean respEmbd, long startTime, long endTime,
      boolean lowFreq) {

    try {
      ds = new DataStore();
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

  }
  
  /**
   * Return the randomized experiment data is being run on.
   * This should not be called until populateDataAndRun(..) has been.
   * @return Randomized experiment, to enable reading the results.
   */
  public RandomizedExperiment getExperiment() {
    return re;
  }
  
  public static void main(String[] args) {
    GatewayServer gatewayServer = new GatewayServer( new RandomCalShell() );
    gatewayServer.start();
    System.out.println("Gateway Server Started");
  }

}
