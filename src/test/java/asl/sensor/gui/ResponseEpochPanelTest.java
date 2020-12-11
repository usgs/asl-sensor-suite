package asl.sensor.gui;

import static asl.sensor.gui.ResponseEpochPanel.compareNamesPossiblyMissingComponents;
import static asl.utils.response.ChannelMetadata.RESP_DT_FORMAT;
import static org.junit.Assert.assertEquals;

import asl.sensor.test.TestUtils;
import asl.utils.response.ResponseParser;
import asl.utils.response.ResponseParser.EpochIdentifier;
import java.io.IOException;
import java.time.Instant;
import java.util.List;
import org.junit.Before;
import org.junit.Test;

public class ResponseEpochPanelTest {

  private static final String folder = TestUtils.TEST_DATA_LOCATION + "metadata/rdseed/";
  private static List<EpochIdentifier> epochs;

  @Before
  public void setUp() throws IOException {
    String filename = folder + "/IU-ANMO-ascii.txt";
    epochs = ResponseParser.listEpochsForSelection(filename);
  }

  @Test
  public void testGetsCorrectIndex_withinEpochInList() throws IOException {
    ResponseEpochPanel panel = new ResponseEpochPanel(epochs);
    long start = Instant.from(RESP_DT_FORMAT.parse("2015,206,00:00:00")).toEpochMilli();
    long end = Instant.from(RESP_DT_FORMAT.parse("2015,207,00:00:00")).toEpochMilli();
    panel.setIndexOfClosestEpoch(start, end, "IU.ANMO.00.LH1");
    assertEquals("IU.ANMO.00.LH1: 2014,351,18:40:00 | 2018,190,20:45:00",
        panel.getSelectedEpoch().toString());
  }

  @Test
  public void testGetsCorrectIndex_beforeEpochInList() throws IOException {
    ResponseEpochPanel panel = new ResponseEpochPanel(epochs);
    long start = Instant.from(RESP_DT_FORMAT.parse("1997,206,00:00:00")).toEpochMilli();
    long end = Instant.from(RESP_DT_FORMAT.parse("1997,207,00:00:00")).toEpochMilli();
    panel.setIndexOfClosestEpoch(start, end, "IU.ANMO.00.LH1");
    assertEquals("IU.ANMO.00.LH1: 1998,299,20:00:00 | 2000,293,16:00:00",
        panel.getSelectedEpoch().toString());
  }

  @Test
  public void testGetsCorrectIndex_afterLastEpochInList() throws IOException {
    ResponseEpochPanel panel = new ResponseEpochPanel(epochs);
    long start = Instant.from(RESP_DT_FORMAT.parse("2020,206,00:00:00")).toEpochMilli();
    long end = Instant.from(RESP_DT_FORMAT.parse("2020,207,00:00:00")).toEpochMilli();
    panel.setIndexOfClosestEpoch(start, end, "IU.ANMO.00.LH1");
    assertEquals("IU.ANMO.00.LH1: 2018,190,20:45:00 | ",
        panel.getSelectedEpoch().toString());
  }

  @Test
  public void testGetsCorrectIndex_notInListWithinBoundary() throws IOException {
    ResponseEpochPanel panel = new ResponseEpochPanel(epochs);
    long start = Instant.from(RESP_DT_FORMAT.parse("2020,206,00:00:00")).toEpochMilli();
    long end = Instant.from(RESP_DT_FORMAT.parse("2020,207,00:00:00")).toEpochMilli();
    panel.setIndexOfClosestEpoch(start, end, "IU.ANMO.00.LLL");
    assertEquals("IU.ANMO. .BC0: 2000,293,16:00:00 | 2002,323,21:07:00",
        panel.getSelectedEpoch().toString());
  }

  @Test
  public void testGetsCorrectIndex_notInListBeforeBoundary() throws IOException {
    ResponseEpochPanel panel = new ResponseEpochPanel(epochs);
    long start = Instant.from(RESP_DT_FORMAT.parse("2020,206,00:00:00")).toEpochMilli();
    long end = Instant.from(RESP_DT_FORMAT.parse("2020,207,00:00:00")).toEpochMilli();
    panel.setIndexOfClosestEpoch(start, end, "AA.AAAA.AA.AAA");
    assertEquals("IU.ANMO. .BC0: 2000,293,16:00:00 | 2002,323,21:07:00",
        panel.getSelectedEpoch().toString());
  }

  @Test
  public void testGetsCorrectIndex_notInListPastBoundary() throws IOException {
    ResponseEpochPanel panel = new ResponseEpochPanel(epochs);
    long start = Instant.from(RESP_DT_FORMAT.parse("2020,206,00:00:00")).toEpochMilli();
    long end = Instant.from(RESP_DT_FORMAT.parse("2020,207,00:00:00")).toEpochMilli();
    panel.setIndexOfClosestEpoch(start, end, "IU.ANTO.00.LH1");
    assertEquals("IU.ANMO. .BC0: 2000,293,16:00:00 | 2002,323,21:07:00",
        panel.getSelectedEpoch().toString());
  }

  @Test
  public void compareStringsMissingComponents_stringsShouldMatch() {
    String fullName = "IU.ANMO.00.LH1";
    String truncated = "00.LH1";
    assertEquals(0, compareNamesPossiblyMissingComponents(fullName, truncated));
  }

  @Test
  public void compareStringsMissingComponents_fullStringComesFirst() {
    String fullName = "IU.ANMO.00.LH1";
    String truncated = "00.LH2";
    assertEquals(-1, compareNamesPossiblyMissingComponents(fullName, truncated));
  }

  @Test
  public void compareStringsMissingComponents_fullStringComesSecond() {
    String fullName = "IU.ANMO.00.LH2";
    String truncated = "00.LH1";
    assertEquals(1, compareNamesPossiblyMissingComponents(fullName, truncated));
  }
}
