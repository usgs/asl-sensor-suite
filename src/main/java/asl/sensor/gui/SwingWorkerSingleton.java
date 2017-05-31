package asl.sensor.gui;

import java.util.concurrent.ExecutionException;

import javax.swing.SwingWorker;

import asl.sensor.input.DataStore;

/**
 * Used as singleton instance of swingworker to prevent clogging the program
 * with many dead threads and make sure experiments cancel gracefully
 * @author akearns
 *
 */
public class SwingWorkerSingleton {

  private static SwingWorker<Boolean, Void> worker;
  private static ExperimentPanel epHandle;
  
  private SwingWorkerSingleton() {
    // empty constructor; worker is built when experiment is passed into it
  }
  
  /**
   * Used to run experiment panel's backend in the background as it should be
   * As long as this singleton is used to run an experiment, any currently
   * running experiment will be cancelled, meaning only one is run at a time.
   * This should improve performance and prevent experiments from competing
   * with each other for processing time, especially multiple runs of the same
   * experiment.
   * The returned result is a boolean used to determine if the data in the 
   * experiment's display panel was properly set.
   * @param active ExperimentPanel to run calculations from
   * @param ds DataStore whose data will be used in the calculations
   */
  public static void setInstance(ExperimentPanel active, DataStore ds) 
      {
    
    if (worker != null) {
      // clear out any old data in the chart
      // since we only have one worker thread for experiment calculations
      // if we run a new experiment while another one was calculating,
      // the result won't actually complete, so we should make it clear that
      // other panel was cancelled, and thus clear the chart / unset data
      if ( !worker.isDone() ) {
        worker.cancel(true); // cancel worker, set it to the new task
      }
    }
     
    epHandle = active;
    epHandle.clearChartAndSetProgressData();
    
    worker = new SwingWorker<Boolean, Void>() {
      @Override
      protected Boolean doInBackground() {
        epHandle.updateData(ds); 
        // calculate backend and get chart, insets to show
        return epHandle.set;
      }
      
      @Override
      protected void done() {
        try {
          boolean set = get();
          if (set) {
            // display the results of experiment in the panel
            epHandle.setDone(); 
          } 
        } catch (Exception ex) {
          epHandle.displayErrorMessage( ex.getMessage() );
          ex.printStackTrace();
        }
      }
    };

    worker.execute();
  }
  
  public static SwingWorker<Boolean, Void> getInstance() {
    return worker;
  }
  
  /**
   * Get the result of running the swingworker, that is, completion status
   * @return True if the worker (i.e., experiment) completed successfully
   * @throws InterruptedException
   * @throws ExecutionException
   */
  public static boolean getCompleted() 
      throws InterruptedException, ExecutionException {
    return worker.get();
  }
  
}
