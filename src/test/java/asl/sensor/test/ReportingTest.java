package asl.sensor.test;

import org.apache.pdfbox.pdmodel.PDDocument;
import org.junit.Test;

import asl.sensor.utils.ReportingUtils;

public class ReportingTest {

  @Test
  public void testStringListConversion() {
    String[] toReport = new String[]{};
    PDDocument pdf = new PDDocument();
    ReportingUtils.textListToPDFPages(pdf, toReport);
  }
  
}
