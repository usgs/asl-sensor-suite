package asl.sensor.gui;

import asl.sensor.input.DataStore;
import java.util.concurrent.CancellationException;
import java.util.concurrent.ExecutionException;
import javax.swing.SwingWorker;
import org.apache.commons.math3.exception.ConvergenceException;
import org.apache.commons.math3.exception.NoDataException;

/**
 * Used as singleton instance of swingworker to prevent clogging the program
 * with many dead threads and make sure experiments cancel gracefully
 *
 * @author akearns
 */
public class SwingWorkerSingleton {

  private static SwingWorker<Boolean, Void> worker;
  private static ExperimentPanel epHandle;

  private SwingWorkerSingleton() {
    // empty constructor; worker is built when experiment is passed into it
  }

  public static SwingWorker<Boolean, Void> getInstance() {
    return worker;
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
   *
   * @param active ExperimentPanel to run calculations from
   * @param ds DataStore whose data will be used in the calculations
   */
  public static void setInstance(ExperimentPanel active, DataStore ds) {

    if (worker != null) {
      // clear out any old data in the chart
      // since we only have one worker thread for experiment calculations
      // if we run a new experiment while another one was calculating,
      // the result won't actually complete, so we should make it clear that
      // other panel was cancelled, and thus clear the chart / unset data
      if (!worker.isDone()) {
        try {
          worker.cancel(true); // cancel worker, set it to the new task
        } catch (CancellationException ignore) {
        }
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
        } catch (ExecutionException ex) {
          Throwable cause = ex.getCause();
          StringBuilder text = new StringBuilder();
          if (cause instanceof NoDataException) {
            text.append("Solver was not given proper input to solve against.\n");
            text.append("If you are running high-frequency cal, check the data sample rate.\n");
            text.append("(Is the data sampled at a rate less than the frequencies of the");
            text.append(" poles defining the rolloff?)\n");
            text.append("A more detailed explanation has been output to terminal,\n");
            text.append("but here is the error message returned by the backend:\n");
            text.append(cause.toString());
          } else if (cause instanceof ConvergenceException) {
            text.append("Solver was unable to converge on a specific value.\n");
            text.append("A common cause of this is having too little timeseries data to");
            text.append(" process.\n");
            text.append("(Are you running a low-frequency cal on high-frequency data?)\n");
            text.append("A more detailed explanation has been output to terminal,\n");
            text.append("but here is the error message returned by the backend:\n");
            text.append(cause.toString());
          } else if (cause instanceof ArrayIndexOutOfBoundsException) {
            text.append("Solver attempted to access nonexistent data\n");
            text.append("A common cause of this is having too little timeseries data to");
            text.append(" process.\n");
            text.append("(Are you running a low-frequency cal on a small amount of data?)\n");
            text.append("A more detailed explanation has been output to terminal,\n");
            text.append("but here is the error message returned by the backend:\n");
            text.append(cause.toString());
          } else {
            if (ex.getMessage() == null) {
              text.append("CANCELLED");
            } else {
              text.append(cause.getMessage());
            }
          }
          epHandle.displayErrorMessage(text.toString());
          cause.printStackTrace();
        } catch (InterruptedException ex) {
          String text;
          if (ex.getMessage() == null) {
            text = "CANCELLED";
          } else {
            text = ex.getMessage();
          }
          epHandle.displayErrorMessage(text);
          ex.getCause().printStackTrace();
        }
      }
    };

    worker.execute();
  }

}
