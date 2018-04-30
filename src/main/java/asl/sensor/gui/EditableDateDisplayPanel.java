package asl.sensor.gui;

import asl.sensor.utils.TimeSeriesUtils;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.time.Instant;
import java.time.OffsetDateTime;
import java.time.ZoneOffset;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JSpinner;
import javax.swing.SpinnerNumberModel;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.event.EventListenerList;

/**
 * This panel presents an alternate means of range-trimming to the standard
 * slider. Given timeseries data, this produces x-axis values as time longs that
 * can be used to set data boundaries. The data is edited by manipulating
 * calendar values from inside this object's spinners, each one corresponding
 * to a date component between the year and the millisecond inclusive.
 *
 * @author akearns - KBRWyle
 */
public class EditableDateDisplayPanel extends JPanel implements ChangeListener {

  private static final long serialVersionUID = -1649983797482938586L;

  private final JSpinner year;
  private final JSpinner day;
  private final JSpinner hour;
  private final JSpinner minute;
  private final JSpinner second;
  private final JSpinner millisecond;
  private final EventListenerList listeners;
  private OffsetDateTime currentTime; // holds the data for the current time

  /**
   * Constructor initializes components, layouts, and listeners.
   */
  EditableDateDisplayPanel() {

    listeners = new EventListenerList();

    currentTime = OffsetDateTime.now(ZoneOffset.UTC);

    SpinnerNumberModel model = new SpinnerNumberModel();
    model.setStepSize(1);
    model.setMinimum(0);
    model.setMaximum(9999); // TODO: in the year 9000, change this?
    year = new JSpinner(model);
    // remove commas from date display
    JSpinner.NumberEditor editor = new JSpinner.NumberEditor(year, "#");
    year.setEditor(editor);
    year.addChangeListener(this);
    JLabel yLabel = new JLabel("(Y)");

    model = new SpinnerNumberModel();
    model.setStepSize(1);
    model.setMinimum(0);
    model.setMaximum(366);
    day = new JSpinner(model);
    day.addChangeListener(this);
    JLabel dLabel = new JLabel("(D)");

    model = new SpinnerNumberModel();
    model.setStepSize(1);
    model.setMinimum(0);
    model.setMaximum(23);
    hour = new JSpinner(model);
    hour.addChangeListener(this);
    JLabel hLabel = new JLabel("(H-24)");

    model = new SpinnerNumberModel();
    model.setStepSize(1);
    model.setMinimum(0);
    model.setMaximum(59);
    minute = new JSpinner(model);
    minute.addChangeListener(this);
    JLabel mLabel = new JLabel("(M)");

    model = new SpinnerNumberModel();
    model.setStepSize(1);
    model.setMinimum(0);
    model.setMaximum(59);
    second = new JSpinner(model);
    second.addChangeListener(this);
    JLabel sLabel = new JLabel("(S)");

    model = new SpinnerNumberModel();
    model.setStepSize(1);
    model.setMinimum(0);
    model.setMaximum(999);
    millisecond = new JSpinner(model);
    millisecond.addChangeListener(this);
    JLabel msLabel = new JLabel("(ms)");

    // build panel
    GridBagConstraints constraints = new GridBagConstraints();
    this.setLayout(new GridBagLayout());

    constraints.weightx = 1.0;
    constraints.gridx = 0;
    constraints.gridy = 0;
    this.add(year, constraints);
    constraints.weightx = 0.;
    constraints.gridy = 1;
    this.add(yLabel, constraints);
    constraints.weightx = 1.0;
    constraints.gridy = 0;
    constraints.gridx += 1;
    this.add(day, constraints);
    constraints.weightx = 0.;
    constraints.gridy = 1;
    this.add(dLabel, constraints);
    constraints.weightx = 1.0;
    constraints.gridy = 0;
    constraints.gridx += 1;
    this.add(hour, constraints);
    constraints.weightx = 0.;
    constraints.gridy = 1;
    this.add(hLabel, constraints);
    constraints.weightx = 1.0;
    constraints.gridy = 0;
    constraints.gridx += 1;
    this.add(minute, constraints);
    constraints.weightx = 0.;
    constraints.gridy = 1;
    this.add(mLabel, constraints);
    constraints.weightx = 1.0;
    constraints.gridy = 0;
    constraints.gridx += 1;
    this.add(second, constraints);
    constraints.weightx = 0.;
    constraints.gridy = 1;
    this.add(sLabel, constraints);
    constraints.weightx = 1.0;
    constraints.gridy = 0;
    constraints.gridx += 1;
    this.add(millisecond, constraints);
    constraints.weightx = 0.;
    constraints.gridy = 1;
    this.add(msLabel, constraints);

  }

  /**
   * Register a Swing component to be notified when this object changes
   * (add it to the callback list)
   *
   * @param listener swing component to be notified of this panel's state
   */
  public void addChangeListener(ChangeListener listener) {
    listeners.add(ChangeListener.class, listener);
  }

  /**
   * Notifies listening components that the state of this object has changed
   */
  private void fireStateChanged() {
    ChangeListener[] listeners = this.listeners.getListeners(ChangeListener.class);
    if (listeners != null && listeners.length > 0) {
      ChangeEvent evt = new ChangeEvent(this);
      for (ChangeListener listener : listeners) {
        listener.stateChanged(evt);
      }
    }
  }

  /**
   * Get the time specified by this panel's components in milliseconds
   *
   * @return The time based on the calendar components
   */
  public long getTime() {
    return currentTime.toInstant().toEpochMilli();
  }

  /**
   * Remove a Swing component from the list of active listeners to this panel
   *
   * @param listener Swing component no longer listening to this object
   */
  public void removeChangeListener(ChangeListener listener) {
    listeners.remove(ChangeListener.class, listener);
  }

  /**
   * Convert time in milliseconds to calendar components to populate panel
   *
   * @param timeStamp Time to set this panel to
   */
  public void setValues(long timeStamp) {

    if (timeStamp == getTime()) {
      return; // don't do anything if no change is necessary
    }

    currentTime = OffsetDateTime.ofInstant(Instant.ofEpochMilli(timeStamp), ZoneOffset.UTC);
    year.setValue(currentTime.getYear());
    day.setValue(currentTime.getDayOfYear());
    hour.setValue(currentTime.getHour());
    minute.setValue(currentTime.getMinute());
    second.setValue(currentTime.getSecond());
    millisecond.setValue(currentTime.getNano() / TimeSeriesUtils.TO_MILLI_FACTOR);

  }

  @Override
  public void stateChanged(ChangeEvent event) {

    if (event.getSource() == year) {
      currentTime = currentTime.withYear((int) year.getValue());
    } else if (event.getSource() == day) {
      currentTime = currentTime.withDayOfYear((int) day.getValue());
    } else if (event.getSource() == hour) {
      currentTime = currentTime.withHour((int) hour.getValue());
    } else if (event.getSource() == minute) {
      currentTime = currentTime.withMinute((int) minute.getValue());
    } else if (event.getSource() == second) {
      currentTime = currentTime.withSecond((int) second.getValue());
    } else if (event.getSource() == millisecond) {
      currentTime = currentTime.withNano((int) millisecond.getValue() * TimeSeriesUtils.TO_MILLI_FACTOR);
    }

    fireStateChanged(); // percolate change in component up to any containers
  }

}
