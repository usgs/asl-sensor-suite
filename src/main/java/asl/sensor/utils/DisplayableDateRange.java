package asl.sensor.utils;

import java.text.DateFormat;
import java.util.Calendar;
import java.util.TimeZone;

import org.apache.commons.math3.util.Pair;

import asl.sensor.gui.InputPanel;

public class DisplayableDateRange {

  private Pair<Long, Long> dateRangeInternals;
  private String note;
  
  public DisplayableDateRange(Pair<Long, Long> times, String note) {
    long start = Math.min(times.getFirst(), times.getSecond());
    long end = Math.max(times.getFirst(), times.getSecond());
    this.note = note;
    dateRangeInternals = new Pair<Long, Long>(start, end);
  }
  
  public Pair<Long, Long> getDateRange() {
    return dateRangeInternals;
  }
  
  @Override
  public String toString() {
    DateFormat df = InputPanel.SDF;
    df.setTimeZone( TimeZone.getTimeZone("UTC") );
    Calendar cCal = Calendar.getInstance( df.getTimeZone() );
    
    cCal.setTimeInMillis( dateRangeInternals.getFirst() / 1000 );
    StringBuilder sb = new StringBuilder( df.format( cCal.getTime() ) );
    cCal.setTimeInMillis( dateRangeInternals.getSecond() / 1000 );
    sb.append(" | ");
    sb.append( df.format( cCal.getTime() ) );
    sb.append("  ");
    sb.append(note);
    return sb.toString();
  }
}
