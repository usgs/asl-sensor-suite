package asl.sensor.gui;

import java.awt.Dimension;
import java.awt.Rectangle;

import javax.swing.JPanel;
import javax.swing.Scrollable;

/**
 * Simple class to create a composable JPanel that is scrollable only
 * in the vertical direction, as necessary. Used to create an input panel that
 * doesn't get shrunk when too many inputs are needed for an experiment
 * such as for 9-component self-noise.
 * @author akearns
 *
 */
public class VScrollPanel extends JPanel implements Scrollable {

  
  /**
   * 
   */
  private static final long serialVersionUID = -4358405572297138013L;

  boolean scaleVerticallyNotScroll = false;
  
  @Override
  public Dimension getPreferredScrollableViewportSize() {
    return getPreferredSize();
  }

  @Override
  public int getScrollableBlockIncrement(Rectangle visibleRect, int orientation,
      int direction) {
    return 1;
  }

  @Override
  public boolean getScrollableTracksViewportHeight() {
    return scaleVerticallyNotScroll;
  }

  @Override
  public boolean getScrollableTracksViewportWidth() {
    // never set this panel to scroll vertically
    return true;
  }

  @Override
  public int getScrollableUnitIncrement(Rectangle visibleRect, int orientation,
      int direction) {
    return 10;
  }

  /**
   * Sets whether or not to enable scrolling in this panel vs. fitting to size
   * @param bool True if the panel should NOT have scrolling set
   */
  public void setScrollableTracksViewportHeight(boolean bool) {
    scaleVerticallyNotScroll = bool;
  }
  
}
