package asl.sensor.gui;

import static asl.utils.response.ChannelMetadata.RESP_DT_FORMAT;

import asl.utils.response.ResponseParser.EpochIdentifier;
import java.time.Instant;
import java.time.ZoneOffset;
import java.time.format.DateTimeFormatter;
import java.util.Comparator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import javax.swing.BoxLayout;
import javax.swing.DefaultComboBoxModel;
import javax.swing.JComboBox;
import javax.swing.JPanel;
import org.apache.commons.math3.util.Pair;

/**
 * Class to organize the GUI components necessary for choosing epochs for RESP and XML files.
 * For XML files we expect that files may contain multiple channels and so we allow for the
 * channel of interest to be selected, and the epochs of that channel auto-populating the selection.
 * RESP files are expected to only contain one channel (this is ASL convention), and so we have
 * methods to select an epoch without needing to specify the channel in that case.
 * We attempt to print out the epoch time ranges in a way that is more human-readable than the
 * default method for doing so.
 * The list of epochs is automatically sorted according to start time, increasing.
 * @author akearns
 */
public class ReponseEpochPanel extends JPanel {

  // list of epochs associated with the currently selected channel
  // we just keep this separate from the combobox and look up by index
  // because it's easier than trying to create a completely new renderer
  private List<EpochIdentifier> epochs;
  private final JComboBox<String> selectableEpochBox;

  /**
   * Create a chooser for epoch data under the assumption only one channel is specified in metadata.
   * @param gottenEpochs List of epochs
   */
  public ReponseEpochPanel(List<EpochIdentifier> gottenEpochs) {
    epochs = gottenEpochs;
    selectableEpochBox = new JComboBox<>(new EpochComboBoxModel(epochs));
    this.setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));
    this.add(selectableEpochBox);
  }

  public ReponseEpochPanel(Map<String, List<Pair<Instant, Instant>>> epochMap) {
    List<EpochIdentifier> gottenEpochs = new LinkedList<>();
    for (String channel : epochMap.keySet()) {
      List<Pair<Instant, Instant>> startsAndEnds = epochMap.get(channel);
      for (Pair<Instant, Instant> startAndEnd : startsAndEnds) {
        gottenEpochs.add(new EpochIdentifier(startAndEnd.getFirst(), startAndEnd.getSecond(),
            channel, 0)); // file pointer not used for xml parsing
      }
    }
    epochs = gottenEpochs;
    selectableEpochBox = new JComboBox<>(new EpochComboBoxModel(epochs));
    this.setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));
    this.add(selectableEpochBox);
  }

  /**
   * Get the selected epoch (for the selected channel)
   * @return The start and end instants of the desired epoch
   */
  public EpochIdentifier getSelectedEpoch() {
    return epochs.get(selectableEpochBox.getSelectedIndex());
  }

  /**
   * Populates the epoch selector with formatted date-times (listed in sorted order by start)
   */
  private static class EpochComboBoxModel extends DefaultComboBoxModel<String> {
      public EpochComboBoxModel(List<EpochIdentifier> epochs) {
        super();
        // get the string formatter for the epochs
        DateTimeFormatter dateTimeFormatter = RESP_DT_FORMAT.withZone(ZoneOffset.UTC);
        // make sure the given epochs are sorted!
        epochs.sort(EpochStartComparator.instance);
        // now turn the epochs into strings that are more easily human-read
        for (EpochIdentifier epoch : epochs) {
          String string = epoch.toString();
          if (epoch.endInstant == null) {
            string += " (OPEN)";
          }
          this.addElement(epoch.toString());
        }
      }

    /**
     * Sorts the lists of epochs by start time increasing
     */
    private static class EpochStartComparator implements Comparator<EpochIdentifier> {

      static final EpochStartComparator instance = new EpochStartComparator();

      private EpochStartComparator() {
      }

      @Override
      public int compare(EpochIdentifier o1, EpochIdentifier o2) {
        return o1.startInstant.compareTo(o2.startInstant);
      }
    }
  }
}
