package asl.sensor.gui;

import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.util.Calendar;
import java.util.TimeZone;

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
 * @author akearns
 *
 */
public class EditableDateDisplayPanel extends JPanel implements ChangeListener {

  /**
   * 
   */
  private static final long serialVersionUID = -1649983797482938586L;
  
  private JSpinner year, day, hour, minute, second, millisecond;
  private JLabel yLabel, dLabel, hLabel, mLabel, sLabel, msLabel;
  private EventListenerList listeners;
  private Calendar c; // holds the data for the current time
  
  /**
   * Constructor initializes components, layouts, and listeners.
   */
  public EditableDateDisplayPanel() {
    
    listeners = new EventListenerList();
    
    c = Calendar.getInstance();
    c.setTimeZone( TimeZone.getTimeZone("UTC") );
    
    SpinnerNumberModel model = new SpinnerNumberModel();
    model.setStepSize(1);
    model.setMinimum(0);
    model.setMaximum(9999); // TODO: in the year 9000, change this
    year = new JSpinner(model);
    // remove commas from date display
    JSpinner.NumberEditor editor = new JSpinner.NumberEditor(year, "#");
    year.setEditor(editor);
    year.addChangeListener(this);
    yLabel = new JLabel("(Y)");
    
    model = new SpinnerNumberModel();
    model.setStepSize(1);
    model.setMinimum(0);
    model.setMaximum(366);
    day = new JSpinner(model);
    day.addChangeListener(this);
    dLabel = new JLabel("(D)");
    
    model = new SpinnerNumberModel();
    model.setStepSize(1);
    model.setMinimum(0);
    model.setMaximum(23);
    hour = new JSpinner(model);
    hour.addChangeListener(this);
    hLabel = new JLabel("(H-24)");
    
    model = new SpinnerNumberModel();
    model.setStepSize(1);
    model.setMinimum(0);
    model.setMaximum(59);
    minute = new JSpinner(model);
    minute.addChangeListener(this);
    mLabel = new JLabel("(M)");
    
    model = new SpinnerNumberModel();
    model.setStepSize(1);
    model.setMinimum(0);
    model.setMaximum(59);
    second = new JSpinner(model);
    second.addChangeListener(this);
    sLabel = new JLabel("(S)");
    
    model = new SpinnerNumberModel();
    model.setStepSize(1);
    model.setMinimum(0);
    model.setMaximum(999);
    millisecond = new JSpinner(model);
    millisecond.addChangeListener(this);
    msLabel = new JLabel("(ms)");
    
    // build panel
    GridBagConstraints gbc = new GridBagConstraints();
    this.setLayout( new GridBagLayout() );
    
    gbc.weightx = 1.0;
    gbc.gridx = 0;
    gbc.gridy = 0;
    this.add(year, gbc);
    gbc.weightx = 0.;
    gbc.gridy = 1;
    this.add(yLabel, gbc);
    gbc.weightx = 1.0;
    gbc.gridy = 0;
    gbc.gridx += 1;
    this.add(day, gbc);
    gbc.weightx = 0.;
    gbc.gridy = 1;
    this.add(dLabel, gbc);
    gbc.weightx = 1.0;
    gbc.gridy = 0;
    gbc.gridx += 1;
    this.add(hour, gbc);
    gbc.weightx = 0.;
    gbc.gridy = 1;
    this.add(hLabel, gbc);
    gbc.weightx = 1.0;
    gbc.gridy = 0;
    gbc.gridx += 1;
    this.add(minute, gbc);
    gbc.weightx = 0.;
    gbc.gridy = 1;
    this.add(mLabel, gbc);
    gbc.weightx = 1.0;
    gbc.gridy = 0;
    gbc.gridx += 1;
    this.add(second, gbc);
    gbc.weightx = 0.;
    gbc.gridy = 1;
    this.add(sLabel, gbc);
    gbc.weightx = 1.0;
    gbc.gridy = 0;
    gbc.gridx += 1;
    this.add(millisecond, gbc);
    gbc.weightx = 0.;
    gbc.gridy = 1;
    this.add(msLabel, gbc);
    
  }
  
  /**
   * Initialize a panel with a specific initial time value
   * @param timeStamp Java epoch time (ms)
   */
  public EditableDateDisplayPanel(long timeStamp) {
    this();
    setValues(timeStamp);
  }

  /**
   * Register a Swing component to be notified when this object changes
   * (add it to the callback list)
   * @param listener swing component to be notified of this panel's state
   */
  public void addChangeListener(ChangeListener listener) {
    listeners.add(ChangeListener.class, listener);
  }
  
  /**
   * Notifies listening components that the state of this object has changed
   */
  private void fireStateChanged() {
    ChangeListener[] lsners = listeners.getListeners(ChangeListener.class);
    if (lsners != null && lsners.length > 0) {
      ChangeEvent evt = new ChangeEvent(this);
      for (ChangeListener lsnr : lsners) {
        lsnr.stateChanged(evt);
      }
    }
  }
  
  /**
   * Get the time specified by this panel's components in milliseconds
   * @return The time based on the calendar componets
   */
  public long getTime() {
    return c.getTimeInMillis();
  }
  
  /**
   * Remove a Swing component from the list of active listeners to this panel
   * @param listener Swing component no longer listening to this object
   */
  public void removeChangeListener(ChangeListener listener) {
    listeners.remove(ChangeListener.class, listener);
  }
  
  /**
   * Convert time in milliseconds to calendar components to populate panel
   * @param timeStamp Time to set this panel to
   */
  public void setValues(long timeStamp) {
    
    if ( timeStamp == getTime() ) {
      return; // don't do anything if no change is necessary
    }
    
    c.setTimeInMillis(timeStamp);
    year.setValue( c.get(Calendar.YEAR) );
    day.setValue( c.get(Calendar.DAY_OF_YEAR) );
    hour.setValue( c.get(Calendar.HOUR_OF_DAY) );
    minute.setValue( c.get(Calendar.MINUTE) );
    second.setValue( c.get(Calendar.SECOND) );
    millisecond.setValue( c.get(Calendar.MILLISECOND) );
    
  }
  
  @Override
  public void stateChanged(ChangeEvent e) {
    
    if ( e.getSource() == year ) {
      c.set( Calendar.YEAR, (int) year.getValue() );
    } else if ( e.getSource() == day ) {
      c.set( Calendar.DAY_OF_YEAR, (int) day.getValue() );
    } else if ( e.getSource() == hour ) {
      c.set( Calendar.HOUR_OF_DAY, (int) hour.getValue() );
    } else if ( e.getSource() == minute ) {
      c.set( Calendar.MINUTE, (int) minute.getValue() );
    } else if ( e.getSource() == second ) {
      c.set( Calendar.SECOND, (int) second.getValue() );
    } else if ( e.getSource() == millisecond ) {
      c.set( Calendar.MILLISECOND, (int) millisecond.getValue() );
    }
    
    fireStateChanged(); // percolate change in component up to any containers
  }
  
}
