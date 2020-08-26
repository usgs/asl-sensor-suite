package asl.sensor.gui;

import asl.utils.input.InstrumentResponse;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.time.Instant;
import java.time.ZoneOffset;
import java.time.format.DateTimeFormatter;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import javax.swing.BoxLayout;
import javax.swing.DefaultComboBoxModel;
import javax.swing.JComboBox;
import javax.swing.JLabel;
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
public class ReponseEpochPanel extends JPanel implements ActionListener {

  private String[] channels;
  private Map<String, List<Pair<Instant, Instant>>> epochMap;
  private JComboBox<String> channelIDBox;
  private JLabel channelIDName;
  // list of epochs associated with the currently selected channel
  // we just keep this separate from the combobox and look up by index
  // because it's easier than trying to create a completely new renderer
  private List<Pair<Instant, Instant>> epochs;
  private JComboBox<String> selectableEpochBox;

  /**
   * Create a chooser for both channel and epoch from metadata. If there is only one channel
   * specified in the metadata, only the epoch will be selectable
   * @param gottenEpochs Map of epochs keyed to channel names
   */
  public ReponseEpochPanel(Map<String, List<Pair<Instant, Instant>>> gottenEpochs) {
    epochMap = gottenEpochs;
    channels = epochMap.keySet().toArray(new String[]{});
    Arrays.sort(channels);
    channelIDBox = new JComboBox<>(channels);
    channelIDBox.setSelectedIndex(0);
    channelIDBox.addActionListener(this);
    epochs = epochMap.get(channels[0]);
    selectableEpochBox = new JComboBox<>(new EpochComboBoxModel(epochs));

    this.setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));
    if (channels.length > 1) {
      this.add(channelIDBox);
    }
    this.add(selectableEpochBox);
  }

  /**
   * Create a chooser for epoch data under the assumption only one channel is specified in metadata.
   * @param gottenEpochs List of epochs
   */
  public ReponseEpochPanel(List<Pair<Instant, Instant>> gottenEpochs) {
    epochMap = null;
    channels = new String[]{""};
    channelIDBox = new JComboBox<>(channels);
    channelIDBox.setSelectedIndex(0);
    epochs = gottenEpochs;
    selectableEpochBox = new JComboBox<>(new EpochComboBoxModel(epochs));
    this.setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));
    this.add(selectableEpochBox);
  }

  /**
   * Return the name of the currently selected channel based on the map keys
   * @return String representation of the channel's identifier
   */
  public String getSelectedChannelName() {
    return (String) channelIDBox.getSelectedItem();
  }

  /**
   * Get the selected epoch (for the selected channel)
   * @return The start and end instants of the desired epoch
   */
  public Pair<Instant, Instant> getSelectedEpoch() {
    return epochs.get(selectableEpochBox.getSelectedIndex());
  }

  @Override
  public void actionPerformed(ActionEvent e) {
    // actionlistener should only be enabled if there are multiple channels
    if (e.getSource() == channelIDBox) {
      String channel = channels[channelIDBox.getSelectedIndex()];
      selectableEpochBox.removeAllItems();
      epochs = epochMap.get(channel);
      selectableEpochBox.setModel(new EpochComboBoxModel(epochMap.get(channel)));
    }
  }

  /**
   * Populates the epoch selector with formatted date-times (listed in sorted order by start)
   */
  private static class EpochComboBoxModel extends DefaultComboBoxModel<String> {
      public EpochComboBoxModel(List<Pair<Instant, Instant>> epochs) {
        super();
        // get the string formatter for the epochs
        DateTimeFormatter dateTimeFormatter = InstrumentResponse.RESP_DT_FORMAT
            .withZone(ZoneOffset.UTC);
        // make sure the given epochs are sorted!
        epochs.sort(EpochStartComparator.instance);
        // now turn the epochs into strings that are more easily human-read
        for (Pair<Instant, Instant> epoch : epochs) {
          Instant startInstant = epoch.getFirst();
          String startString = dateTimeFormatter.format(startInstant);
          Instant endInstant = epoch.getSecond();
          String endString;
          if (endInstant != null) {
            endString = dateTimeFormatter.format(endInstant);
          } else {
            endString = "(NO EPOCH END DEFINED)";
          }
          this.addElement(startString + " | " + endString);
        }
      }

    /**
     * Sorts the lists of epochs by start time increasing
     */
    private static class EpochStartComparator implements Comparator<Pair<Instant, Instant>> {

      static final EpochStartComparator instance = new EpochStartComparator();

      private EpochStartComparator() {
      }

      @Override
      public int compare(Pair<Instant, Instant> o1, Pair<Instant, Instant> o2) {
        return o1.getFirst().compareTo(o2.getFirst());
      }
    }
  }
}
