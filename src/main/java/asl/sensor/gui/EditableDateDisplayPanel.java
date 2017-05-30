package asl.sensor.gui;

import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.util.Calendar;
import java.util.TimeZone;

import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JSpinner;
import javax.swing.SpinnerNumberModel;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.event.EventListenerList;

public class EditableDateDisplayPanel extends JPanel implements ChangeListener {

  /**
   * 
   */
  private static final long serialVersionUID = -1649983797482938586L;
  
  private JSpinner year, day, hour, minute, second, millisecond;
  private JLabel yLabel, dLabel, hLabel, mLabel, sLabel, msLabel;
  private EventListenerList listeners;
  private Calendar c;
  
  public EditableDateDisplayPanel() {
    
    listeners = new EventListenerList();
    
    c = Calendar.getInstance();
    c.setTimeZone( TimeZone.getTimeZone("UTC") );
    
    SpinnerNumberModel model = new SpinnerNumberModel();
    model.setStepSize(1);
    model.setMinimum(0);
    model.setMaximum(9999); // TODO: change this by year 9000 or so
    year = new JSpinner(model);
    year.addChangeListener(this);
    // remove commas from date display
    JSpinner.NumberEditor editor = new JSpinner.NumberEditor(year, "#");
    year.setEditor(editor);
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
  
  public EditableDateDisplayPanel(long timeStamp) {
    this();
    setValues(timeStamp);
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
  
  public void addChangeListener(ChangeListener listener) {
    listeners.add(ChangeListener.class, listener);
  }
  
  public void removeChangeListener(ChangeListener listener) {
    listeners.remove(ChangeListener.class, listener);
  }
  
  public void fireStateChanged() {
    ChangeListener[] lsners = listeners.getListeners(ChangeListener.class);
    if (lsners != null && lsners.length > 0) {
      ChangeEvent evt = new ChangeEvent(this);
      for (ChangeListener lsnr : lsners) {
        lsnr.stateChanged(evt);
      }
    }
  }
  
  public long getTime() {
    return c.getTimeInMillis();
  }
  
  public void setValues(long timeStamp) {
    c.setTimeInMillis(timeStamp);
    year.setValue( c.get(Calendar.YEAR) );
    day.setValue( c.get(Calendar.DAY_OF_YEAR) );
    hour.setValue( c.get(Calendar.HOUR_OF_DAY) );
    minute.setValue( c.get(Calendar.MINUTE) );
    second.setValue( c.get(Calendar.SECOND) );
    millisecond.setValue( c.get(Calendar.MILLISECOND) );
    fireStateChanged();
  }
  
}
