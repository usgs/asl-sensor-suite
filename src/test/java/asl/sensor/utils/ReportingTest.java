package asl.sensor.utils;

import org.apache.pdfbox.pdmodel.PDDocument;
import org.junit.Test;

public class ReportingTest {

  @Test
  public void testStringListConversion() {
    String[] toReport = new String[]{};
    PDDocument pdf = new PDDocument();
    ReportingUtils.textListToPDFPages(pdf, toReport);
  }

}
