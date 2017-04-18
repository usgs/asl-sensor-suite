package asl.sensor.utils;

import java.awt.Color;
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
   * Add a buffered image to a PDDocument page
   * @param bi BufferedImage to be added to PDF
   * @param pdf PDF to have BufferedImage appended to
   */
  public static void 
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

    return;
    
  }
  
  /**
   * Merge a series of buffered images and write them to a single PDF page
   * @param pdf PDF document to append the page onto
   * @param bis Series of buffered images to merge onto a single page
   */
  public static void 
  bufferedImagesToPDFPage(PDDocument pdf, BufferedImage... bis) {
    BufferedImage toPDF = mergeBufferedImages(bis);
    bufferedImageToPDFPage(toPDF, pdf);
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
   * Takes in a series of charts and produces a PDF page of those charts.
   * For more details on this method, see the chartsToImage function, which
   * this method uses to produce the image to be added to the PDF
   * @param width Width of each chart to be added to the PDF
   * @param height Height of each chart to be added to the PDF
   * @param pdf PDF document to have the data appended to
   * @param jfcs series of charts to place in the PDF
   */
  public static void 
  chartsToPDFPage(int width, int height, PDDocument pdf, JFreeChart... jfcs) {
    
    BufferedImage bi = chartsToImage(width, height, jfcs);
    bufferedImageToPDFPage(bi, pdf);
    return;
    
  }
  
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
   * Add pages to a PDF document consisting of textual data with a series of
   * strings, where each string is written to a separate page
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
   * Writes multiple pages of charts to a PDF file, with the number of charts
   * to display per page set according to a parameter 
   * @param perPage Number of charts to put in a page at a time
   * @param width Width of each chart to write to file
   * @param height Height of each chart to write to file
   * @param pdf Document to append pages of chart plots to
   * @param charts Charts whose plots will be written to PDF
   */
  public static void
  groupChartsToPDFPages(int perPage, int width, int height,
      PDDocument pdf, JFreeChart... charts) {
    
    imageListToPDFPages( pdf, 
        chartsToImageList(perPage, width, height, charts) );
    
  }
  
  /**
   * Write a list of images to a pdf document, each image its own page
   * @param pdf PDF document to write to
   * @param bis List of buffered images to write. Each image is written to its
   * own PDF page.
   */
  public static void
  imageListToPDFPages(PDDocument pdf, BufferedImage... bis) {
    for (BufferedImage bi : bis) {
      bufferedImageToPDFPage(bi, pdf);
    }
  }
  
  /**
   * Create a list of buffered images from a series of charts, with a specified
   * number of charts included on each image. This is used to write the charts
   * to a series of pages in a PDF report
   * @param perImg Number of charts' plots to write to a single page
   * @param width Width to set each chart's output image
   * @param height Height to set each chart's output image
   * @param charts List of charts to be compiled into images
   * @return A list of buffered images with no more than perImg plots in 
   * each image
   */
  public static BufferedImage[]
  chartsToImageList(int perImg, int width, int height, JFreeChart... charts) {
    
    List<BufferedImage> imageList = new ArrayList<BufferedImage>();
    int totalNumber = charts.length;
    
    if (totalNumber < perImg) {
      // if we can fit them all on a single page, then we'll do so
      imageList.add( chartsToImage(width, height, charts) );
      return imageList.toArray( new BufferedImage[]{} );
    }
    
    // want to keep all charts the same size; 
    // if we can't fit them all on a single page, how many pages will have
    // complete charts?
    int numFilledPages = totalNumber / perImg;
    // if we can't fill the last page with plots, how many will it have
    int lastPageChartCount = totalNumber % perImg;
    // how many chart-size blank spaces to keep plots on last page same size
    int spacerCount = perImg - lastPageChartCount;
    // Note that the above variable is not used if lastPageChartCount is zero
    
    // handle all the pages with complete data here
    for (int i = 0; i < numFilledPages; ++i) {
      JFreeChart[] onOnePage = new JFreeChart[perImg];
      for (int j = 0; j < perImg; ++j) {
        onOnePage[j] = charts[(perImg * i) + j];
      }
      imageList.add( chartsToImage(width, height, onOnePage) );
    }
    
    // special case for a non-evenly dividing plot series
    if (lastPageChartCount != 0) {
      int lastIndex = numFilledPages * perImg;
      JFreeChart[] lastPage = new JFreeChart[lastPageChartCount];
      for (int j = 0; j < lastPageChartCount; ++j) {
        lastPage[j] = charts[lastIndex + j];
      }
      BufferedImage lastPageImage = chartsToImage(width, height, lastPage);
      BufferedImage space = createWhitespace(width, height * spacerCount);
      imageList.add( mergeBufferedImages(lastPageImage, space) );
    }
    
    return imageList.toArray( new BufferedImage[]{} );
  }
  
  public static BufferedImage createWhitespace(int width, int height) {
    BufferedImage out = new BufferedImage(width, height, 
        BufferedImage.TYPE_INT_RGB);
    Graphics2D g = out.createGraphics();

    g.setPaint(Color.WHITE);
    g.fillRect ( 0, 0, out.getWidth(), out.getHeight() );
    g.dispose();
    return out;
  }
  
  
  /**
   * Add a page to a PDF document consisting of textual data
   * @param toWrite String to add to a new PDF page
   * @param pdf Document to append the page to
   */
  public static void textToPDFPage(String toWrite, PDDocument pdf) {
    
    if ( toWrite.length() == 0 ) {
      return;
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


      return;
  }

}
