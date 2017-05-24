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
   * running experiment will be cancelled, meaning only one is run at a time
   * This should improve performance and prevent experiments from competing
   * with each other for 
   * @param active ExperimentPanel to run calculations from
   * @param ds DataStore whose data will be used in the calculations
   * @return true if the 
   * @throws ExecutionException 
   * @throws InterruptedException 
   */
  public static boolean setInstance(ExperimentPanel active, DataStore ds) 
      throws InterruptedException, ExecutionException {
    
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
        boolean set = false;
        try {
          set = get();
        } catch (InterruptedException e1) {
          epHandle.displayErrorMessage( e1.getMessage() );
        } catch (ExecutionException e2) {
          epHandle.displayErrorMessage( e2.getMessage() );
        }
        if (set) {
          // display the results of experiment in the panel
          epHandle.drawCharts(); 
        }
      }
    };

    worker.execute();
    
    return worker.get();
  }
  
  
  
}
