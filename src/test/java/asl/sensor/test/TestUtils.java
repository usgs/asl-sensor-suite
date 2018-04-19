package asl.sensor.test;

import java.time.Instant;
import java.time.LocalDateTime;
import java.time.OffsetDateTime;
import java.time.ZoneOffset;
import java.time.format.DateTimeFormatter;
import asl.sensor.input.DataStore;

public class TestUtils {

  // may need to change this to deal with eventual migration to main usgs github
  static String SUBPAGE = "tests/";
  static String TEST_DATA_LOCATION = "src/test/resources/seismic-test-data/";
  static DateTimeFormatter DATE_TIME_FORMAT = DateTimeFormatter.ofPattern("uuuu-DDD'T'HH:mm:ss.S");

  public static OffsetDateTime getEndCalendar(DataStore ds) {
    Instant time = Instant.ofEpochMilli( ds.getBlock(0).getEndTime() );
    return OffsetDateTime.ofInstant(time, ZoneOffset.UTC);
  }

  public static OffsetDateTime getStartCalendar(DataStore ds) {
    Instant time = Instant.ofEpochMilli( ds.getBlock(0).getStartTime() );
    return OffsetDateTime.ofInstant(time, ZoneOffset.UTC);
  }

  public static long timeStringToEpochMilli(String time) {
    return LocalDateTime.parse(time, DATE_TIME_FORMAT).toInstant(ZoneOffset.UTC).toEpochMilli();
  }

}
