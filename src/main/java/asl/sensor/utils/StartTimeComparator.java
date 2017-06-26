package asl.sensor.utils;

import java.util.Comparator;

import org.apache.commons.math3.util.Pair;

public class StartTimeComparator implements Comparator<Pair<Long, Long>> {

  @Override
  public int compare(Pair<Long, Long> p1, Pair<Long, Long> p2) {
    return p1.getFirst().compareTo( p2.getFirst() );
  }

}
