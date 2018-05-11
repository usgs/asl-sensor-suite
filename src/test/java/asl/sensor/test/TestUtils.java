package asl.sensor.test;

import java.time.Instant;
import java.time.LocalDateTime;
import java.time.OffsetDateTime;
import java.time.ZoneOffset;
import java.time.format.DateTimeFormatter;
import asl.sensor.input.DataStore;

public class TestUtils {

  // may need to change this to deal with eventual migration to main usgs github
  public static String SUBPAGE = "tests/";
  public static String TEST_DATA_LOCATION = "src/test/resources/seismic-test-data/";
  public static String RESP_LOCATION = TEST_DATA_LOCATION + "RESPs/";

  public static String getSeedFolder(String network, String station, String year, String dayOfYear) {
    return String.format(TEST_DATA_LOCATION + "seed_data/%s_%s/%s/%s/", network, station, year, dayOfYear);
  }

  public static OffsetDateTime getEndCalendar(DataStore ds) {
    Instant time = Instant.ofEpochMilli( ds.getBlock(0).getEndTime() );
    return OffsetDateTime.ofInstant(time, ZoneOffset.UTC);
  }

  public static OffsetDateTime getStartCalendar(DataStore ds) {
    Instant time = Instant.ofEpochMilli( ds.getBlock(0).getStartTime() );
    return OffsetDateTime.ofInstant(time, ZoneOffset.UTC);
  }

  public static long timeStringToEpochMilli(String time) {
    DateTimeFormatter DATE_TIME_FORMAT = DateTimeFormatter.ofPattern("uuuu-DDD'T'HH:mm:ss.S");
    return LocalDateTime.parse(time, DATE_TIME_FORMAT).toInstant(ZoneOffset.UTC).toEpochMilli();
  }

}