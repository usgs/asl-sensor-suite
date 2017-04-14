package asl.sensor.utils;

import java.awt.Dimension;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.io.IOException;
import java.util.ArrayList;
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
 * @author akearns
 *
 */
public class ReportingUtils {

  /**
   * Utility function to combine a series of buffered images into a single
   * buffered image. Images are concatenated vertically and centered 
   * horizontally into an image as wide as the widest passed-in image
   * @param bis Buffered images to send in
   * @return Single concatenated buffered image
   */
  public static BufferedImage mergeBufferedImages(BufferedImage... bis) {
   
    int maxWidth = 0;
    int totalHeight = 0;
    for (BufferedImage bi : bis) {
      if ( maxWidth < bi.getWidth() ) {
        maxWidth = bi.getWidth();
      }
      totalHeight += bi.getHeight();
    }
    
    BufferedImage out = 
        new BufferedImage(maxWidth, totalHeight, BufferedImage.TYPE_INT_RGB);
    Graphics2D g = out.createGraphics();
    
    int heightIndex = 0;
    for (BufferedImage bi : bis) {
      int centeringOffset = 0; // need to center the component?
      if (bi.getWidth() < maxWidth) {
        centeringOffset = ( maxWidth - bi.getWidth() ) / 2;
      }
      g.drawImage(bi, null, centeringOffset, heightIndex);
      heightIndex += bi.getHeight();
    }
    
    return out;
    
  }
  
  /**
   * Converts a series of charts into a buffered image. Each chart has the
   * dimensions given as the width and height parameters, and so the resulting
   * image has width given by that parameter and height equal to height
   * multiplied by the number of charts passed in 
   * (that is, the charts are concatenated vertically)
   * @param width width of each chart plot
   * @param height height of each chart plot
   * @param jfcs series of charts to be plotted in
   * @return buffered image consisting of the concatenation of the given charts
   */
  public static BufferedImage 
  chartsToImage(int width, int height, JFreeChart... jfcs) {
    
    BufferedImage[] bis = new BufferedImage[jfcs.length];
    
    for (int i = 0; i < jfcs.length; ++i) {
      ChartPanel cp = new ChartPanel(jfcs[i]);
      cp.setSize( new Dimension(width, height) );
      BufferedImage temp = new BufferedImage(
          (int) cp.getWidth(),
          (int) cp.getHeight(),
          BufferedImage.TYPE_INT_ARGB);
      Graphics2D g = temp.createGraphics();
      g = temp.createGraphics();
      cp.printAll(g);
      g.dispose();
      bis[i] = temp;
    }
    
    return mergeBufferedImages(bis);
  }
  
  /**
   * Add a buffered image to a PDDocument page
   * @param bi BufferedImage to be added to PDF
   * @param pdf PDF to have BufferedImage appended to
   * @return PDF with the image included on the page
   */
  public static PDDocument 
  bufferedImageToPDFPage(BufferedImage bi, PDDocument pdf) {
    
    PDRectangle rec = 
        new PDRectangle( (float) bi.getWidth(), 
                         (float) bi.getHeight() );
    PDPage page = new PDPage(rec);
    
    try {
      PDImageXObject  pdImageXObject = 
          LosslessFactory.createFromImage(pdf, bi);
      pdf.addPage(page);
      PDPageContentStream contentStream = 
          new PDPageContentStream(pdf, page, 
                                  PDPageContentStream.AppendMode.OVERWRITE, 
                                  true, false);

      contentStream.drawImage( pdImageXObject, 0, 0, 
          bi.getWidth(), bi.getHeight() );
      contentStream.close();
      
    } catch (IOException e) {
      e.printStackTrace();
    }

    return pdf;
    
  }
  
  /**
   * Add a page to a PDF document consisting of textual data
   * @param toWrite String to add to a new PDF page
   * @param pdf Document to append the page to
   * @return PDF with text page appended
   */
  public static PDDocument
  textToPDFPage(String toWrite, PDDocument pdf) {
    if ( toWrite.length() == 0 ) {
      return pdf;
    }
    
      PDPage page = new PDPage();
      pdf.addPage(page);

      PDFont pdfFont = PDType1Font.COURIER;
      float fontSize = 14;
      float leading = 1.5f * fontSize;

      PDRectangle mediabox = page.getMediaBox();
      float margin = 72;
      float width = mediabox.getWidth() - 2*margin;
      float startX = mediabox.getLowerLeftX() + margin;
      float startY = mediabox.getUpperRightY() - margin;

      List<String> lines = new ArrayList<String>();

      for (String text : toWrite.split("\n") ) {

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
            } else if ( spaceIndex == text.length() ) {
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
        contentStream.beginText();contentStream.setFont(pdfFont, fontSize);
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


      return pdf;
  }
  
  /**
   * Takes in a series of charts and produces a PDF page of those charts.
   * For more details on this method, see the chartsToImage function, which
   * this method uses to produce the image to be added to the PDF
   * @param width Width of each chart to be added to the PDF
   * @param height Height of each chart to be added to the PDF
   * @param pdf PDF document to have the data appended to
   * @param jfcs series of charts to place in the PDF
   * @return PDF document with new page added consisting of the included charts
   */
  public static PDDocument 
  chartsToPDFPage(int width, int height, PDDocument pdf, JFreeChart... jfcs) {
    
    BufferedImage bi = chartsToImage(width, height, jfcs);
    return bufferedImageToPDFPage(bi, pdf);
    
  }

}
