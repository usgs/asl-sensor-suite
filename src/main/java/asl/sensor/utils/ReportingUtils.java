package asl.sensor.utils;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import org.apache.pdfbox.pdmodel.PDDocument;
import org.apache.pdfbox.pdmodel.PDPage;
import org.apache.pdfbox.pdmodel.PDPageContentStream;
import org.apache.pdfbox.pdmodel.common.PDRectangle;
import org.apache.pdfbox.pdmodel.font.PDFont;
import org.apache.pdfbox.pdmodel.font.PDType1Font;
import org.apache.pdfbox.pdmodel.graphics.image.LosslessFactory;
import org.apache.pdfbox.pdmodel.graphics.image.PDImageXObject;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;

/**
 * This class defines functions relevant to creating output files,
 * such as images of plots and PDF reports. These methods are all static.
 * Input, output, and test classes all use these functions, so by placing them
 * in a utility function with static methods the code inside those classes can
 * be simplified and redundant calls or procedures reduced.
 *
 * @author akearns - KBRWyle
 */
public class ReportingUtils {

  /**
   * Add a buffered image to a PDDocument page
   *
   * @param bi BufferedImage to be added to PDF
   * @param pdf PDF to have BufferedImage appended to
   */
  private static void
  bufferedImageToPDFPage(BufferedImage bi, PDDocument pdf) {

    PDRectangle rec =
        new PDRectangle(bi.getWidth(),
            bi.getHeight());
    PDPage page = new PDPage(rec);

    try {
      PDImageXObject pdImageXObject =
          LosslessFactory.createFromImage(pdf, bi);
      pdf.addPage(page);
      PDPageContentStream contentStream =
          new PDPageContentStream(pdf, page,
              PDPageContentStream.AppendMode.OVERWRITE,
              true, false);

      contentStream.drawImage(pdImageXObject, 0, 0,
          bi.getWidth(), bi.getHeight());
      contentStream.close();

    } catch (IOException e) {
      e.printStackTrace();
    }
  }

  /**
   * Converts a series of charts into a buffered image. Each chart has the
   * dimensions given as the width and height parameters, and so the resulting
   * image has width given by that parameter and height equal to height
   * multiplied by the number of charts passed in
   * (that is, the charts are concatenated vertically)
   *
   * @param width width of each chart plot
   * @param height height of each chart plot
   * @param charts series of charts to be plotted in
   * @return buffered image consisting of the concatenation of the given charts
   */
  public static BufferedImage
  chartsToImage(int width, int height, JFreeChart... charts) {

    BufferedImage[] bis = new BufferedImage[charts.length];

    for (int i = 0; i < charts.length; ++i) {
      ChartPanel cp = new ChartPanel(charts[i]);
      cp.setSize(new Dimension(width, height));
      BufferedImage temp = new BufferedImage(
          cp.getWidth(),
          cp.getHeight(),
          BufferedImage.TYPE_INT_ARGB);
      Graphics2D g = temp.createGraphics();
      cp.printAll(g);
      g.dispose();
      bis[i] = temp;
    }

    return mergeBufferedImages(bis);
  }

  /**
   * Create a list of buffered images from a series of charts, with a specified
   * number of charts included on each image. This is used to write the charts
   * to a series of pages in a PDF report
   *
   * @param perImg Number of charts' plots to write to a single page
   * @param width Width to set each chart's output image
   * @param height Height to set each chart's output image
   * @param charts List of charts to be compiled into images
   * @return A list of buffered images with no more than perImg plots in
   * each image
   */
  public static BufferedImage[]
      chartsToImageList(int perImg, int width, int height, JFreeChart... charts) {


    int totalNumber = charts.length;

    if (totalNumber < perImg) {
      // if we can fit them all on a single page, then we'll do so
      return new BufferedImage[] {chartsToImage(width, height, charts)};
    }

    // want to keep all charts the same size;
    // if we can't fit them all on a single page, how many pages will have
    // complete charts?
    int numFilledPages = totalNumber / perImg;
    BufferedImage[] imageList = new BufferedImage[numFilledPages];
    // if we can't fill the last page with plots, how many will it have
    int lastPageChartCount = totalNumber % perImg;
    // how many chart-size blank spaces to keep plots on last page same size
    int spacerCount = perImg - lastPageChartCount;

    // handle all the pages with complete data here
    for (int i = 0; i < numFilledPages-1; ++i) {
      int start = perImg * i;
      JFreeChart[] onOnePage = Arrays.copyOfRange(charts, start, start + perImg);
      imageList[i] = chartsToImage(width, height, onOnePage);
    }

    int lastIndex = numFilledPages * perImg;
    JFreeChart[] lastPage = Arrays.copyOfRange(charts, lastIndex, charts.length);
    BufferedImage lastPageImage = chartsToImage(width, height, lastPage);
    // special case for a non-evenly dividing plot series (append blank space to last page)
    if (spacerCount > 0) {
      BufferedImage space = createWhitespace(width, height * spacerCount);
      lastPageImage = mergeBufferedImages(lastPageImage, space);
    }
    imageList[numFilledPages - 1] = lastPageImage;

    return imageList;
  }

  /**
   * Takes in a series of charts and produces a PDF page of those charts.
   * For more details on this method, see the chartsToImage function, which
   * this method uses to produce the image to be added to the PDF
   *
   * @param width Width of each chart to be added to the PDF
   * @param height Height of each chart to be added to the PDF
   * @param pdf PDF document to have the data appended to
   * @param charts series of charts to place in the PDF
   */
  public static void
  chartsToPDFPage(int width, int height, PDDocument pdf, JFreeChart... charts) {

    BufferedImage bi = chartsToImage(width, height, charts);
    bufferedImageToPDFPage(bi, pdf);
  }

  private static BufferedImage createWhitespace(int width, int height) {
    BufferedImage out = new BufferedImage(width, height,
        BufferedImage.TYPE_INT_RGB);
    Graphics2D g = out.createGraphics();

    g.setPaint(Color.WHITE);
    g.fillRect(0, 0, out.getWidth(), out.getHeight());
    g.dispose();
    return out;
  }

  /**
   * Write a list of images to a pdf document, each image its own page
   *
   * @param pdf PDF document to write to
   * @param images List of buffered images to write. Each image is written to its
   * own PDF page.
   */
  public static void
  imageListToPDFPages(PDDocument pdf, BufferedImage... images) {
    for (BufferedImage bi : images) {
      bufferedImageToPDFPage(bi, pdf);
    }
  }

  /**
   * Utility function to combine a series of buffered images into a single
   * buffered image. Images are concatenated vertically and centered
   * horizontally into an image as wide as the widest passed-in image
   *
   * @param images Buffered images to send in
   * @return Single concatenated buffered image
   */
  private static BufferedImage mergeBufferedImages(BufferedImage... images) {

    int maxWidth = 0;
    int totalHeight = 0;
    for (BufferedImage bi : images) {
      if (maxWidth < bi.getWidth()) {
        maxWidth = bi.getWidth();
      }
      totalHeight += bi.getHeight();
    }

    BufferedImage out =
        new BufferedImage(maxWidth, totalHeight, BufferedImage.TYPE_INT_RGB);
    Graphics2D g = out.createGraphics();

    int heightIndex = 0;
    for (BufferedImage bi : images) {
      int centeringOffset = 0; // need to center the component?
      if (bi.getWidth() < maxWidth) {
        centeringOffset = (maxWidth - bi.getWidth()) / 2;
      }
      g.drawImage(bi, null, centeringOffset, heightIndex);
      heightIndex += bi.getHeight();
    }

    return out;

  }

  /**
   * Add pages to a PDF document consisting of textual data with a series of
   * strings, where each string is written to a separate page
   *
   * @param pdf Document to append pages of text to
   * @param toWrite Series of strings to write to PDF
   */
  public static void
  textListToPDFPages(PDDocument pdf, String... toWrite) {

    for (String onePage : toWrite) {
      textToPDFPage(onePage, pdf);
    }

  }


  /**
   * Add a page to a PDF document consisting of textual data
   *
   * @param toWrite String to add to a new PDF page
   * @param pdf Document to append the page to
   */
  public static void textToPDFPage(String toWrite, PDDocument pdf) {

    if (toWrite.length() == 0) {
      return;
    }

    PDPage page = new PDPage();
    pdf.addPage(page);

    PDFont pdfFont = PDType1Font.COURIER;
    float fontSize = 12;
    float leading = 1.5f * fontSize;

    PDRectangle mediaBox = page.getMediaBox();
    float margin = 72;
    float width = mediaBox.getWidth() - 2 * margin;
    float startX = mediaBox.getLowerLeftX() + margin;
    float startY = mediaBox.getUpperRightY() - margin;

    List<String> lines = new ArrayList<>();

    for (String text : toWrite.split("\n")) {

      int lastSpace = -1;
      while (text.length() > 0) {

        int spaceIndex = text.indexOf(' ', lastSpace + 1);
        if (spaceIndex < 0) {
          spaceIndex = text.length();

        }
        String subString = text.substring(0, spaceIndex);
        float size;
        try {
          size = fontSize * pdfFont.getStringWidth(subString) / 1000;
          if (size > width) {
            if (lastSpace < 0) {
              lastSpace = spaceIndex;
            }
            subString = text.substring(0, lastSpace);
            lines.add(subString);
            text = text.substring(lastSpace).trim();
            // System.out.printf("'%s' is line\n", subString);
            lastSpace = -1;
          } else if (spaceIndex == text.length()) {
            lines.add(text);
            // System.out.printf("'%s' is line\n", text);
            text = "";
          } else {
            lastSpace = spaceIndex;
          }
        } catch (IOException e) {
          e.printStackTrace();
        }
      }
    }

    try {
      PDPageContentStream contentStream;
      contentStream = new PDPageContentStream(pdf, page);
      contentStream.beginText();
      contentStream.setFont(pdfFont, fontSize);
      contentStream.newLineAtOffset(startX, startY);
      for (String line : lines) {
        contentStream.showText(line);
        contentStream.newLineAtOffset(0, -leading);
      }
      contentStream.endText();
      contentStream.close();
    } catch (IOException e) {
      e.printStackTrace();
    }
  }
}
