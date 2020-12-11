package asl.sensor.gui;

import asl.utils.response.ResponseParser.EpochIdentifier;
import java.time.Instant;
import java.util.Arrays;
import java.util.Collections;
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
public class ResponseEpochPanel extends JPanel {

  // list of epochs associated with the currently selected channel
  // we just keep this separate from the combobox and look up by index
  // because it's easier than trying to create a completely new renderer
  private final List<EpochIdentifier> epochs;
  private final JComboBox<String> selectableEpochBox;

  /**
   * Create a chooser for epoch data under the assumption only one channel is specified in metadata.
   * @param gottenEpochs List of epochs
   */
  public ResponseEpochPanel(List<EpochIdentifier> gottenEpochs) {
    epochs = gottenEpochs;
    selectableEpochBox = new JComboBox<>(new EpochComboBoxModel(epochs));
    this.setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));
    this.add(selectableEpochBox);
  }

  public ResponseEpochPanel(Map<String, List<Pair<Instant, Instant>>> epochMap) {
    List<EpochIdentifier> gottenEpochs = new LinkedList<>();
    for (String channel : epochMap.keySet()) {
      List<Pair<Instant, Instant>> startsAndEnds = epochMap.get(channel);
      for (Pair<Instant, Instant> startAndEnd : startsAndEnds) {
        gottenEpochs.add(new EpochIdentifier(startAndEnd.getFirst(), startAndEnd.getSecond(),
            channel, 0)); // file pointer not used for xml parsing
      }
    }
    epochs = gottenEpochs;
    epochs.sort(EpochComparator.instance);
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

  public void setIndexOfClosestEpoch(long start, long end, String channelName) {
    // first of all, let's see if we can find the channel name in the epoch list
    EpochIdentifier epochDummy = new EpochIdentifier(Instant.ofEpochMilli(start),
        Instant.ofEpochMilli(end), channelName, 0);
    int goodGuess = Collections.binarySearch(epochs, epochDummy, EpochComparator.instance);
    if (goodGuess >= 0) {
      // this is almost certainly never going to happen, because response epochs can be years long
      selectableEpochBox.setSelectedIndex(goodGuess);
      return;
    }

    // OK, so the epoch isn't in the list, so let's see how many (if any) epochs match the
    // given channel name
    goodGuess = (goodGuess + 1) * -1;
    // subtract one because otherwise we'll get the epoch after the one it slots into
    --goodGuess;
    if (goodGuess < 0) return; // but of course it has to be
    if (!epochs.get(goodGuess).channelIdentifier.equals(channelName)) {
      // the epoch we expected doesn't actually match the channel name, which might be because it
      // comes before even the first epoch in the response file, so let's try that one again
      ++goodGuess; // also, this can push us past the last entry in the list
      if (goodGuess >= epochs.size() ||
          !epochs.get(goodGuess).channelIdentifier.equals(channelName)) {
        // if this expected index doesn't contain the current channel, then it's not in the list
        // and so we'll just go with the default option (i.e., 0)
        return;
      }
    }
    selectableEpochBox.setSelectedIndex(goodGuess);
  }

  /**
   * Populates the epoch selector with formatted date-times (listed in sorted order by start)
   */
  private static class EpochComboBoxModel extends DefaultComboBoxModel<String> {
      public EpochComboBoxModel(List<EpochIdentifier> epochs) {
        super();
        // make sure the given epochs are sorted!
        epochs.sort(EpochComparator.instance);
        // now turn the epochs into strings that are more easily human-read
        for (EpochIdentifier epoch : epochs) {
          this.addElement(epoch.toString());
        }
      }
  }

  static int compareNamesPossiblyMissingComponents(String firstName, String secondName) {
    // in cases where a resp file doesn't include, say, network and station information,
    // we will still want to try to match what information (i.e., channel name) that we have
    if (firstName != null && secondName != null) {
      String[] wordsFromFirst = firstName.split("\\.");
      String[] wordsFromSecond = secondName.split("\\.");
      if (wordsFromFirst.length != wordsFromSecond.length) {
        int shorterLength = Math.min(wordsFromFirst.length, wordsFromSecond.length);
        int firstOffset = wordsFromFirst.length - shorterLength;
        int secondOffset = wordsFromSecond.length - shorterLength;
        wordsFromFirst = Arrays.copyOfRange(wordsFromFirst, firstOffset, wordsFromFirst.length);
        wordsFromSecond = Arrays.copyOfRange(wordsFromSecond, secondOffset, wordsFromSecond.length);
        StringBuilder buildFirst = new StringBuilder();
        StringBuilder buildSecond = new StringBuilder();
        for (int i = 0; i < wordsFromFirst.length; ++i) {
          buildFirst.append(wordsFromFirst[i]).append('.');
          buildSecond.append(wordsFromSecond[i]).append('.');
        }
        firstName = buildFirst.toString();
        secondName = buildSecond.toString();
      }
    } else if (firstName == null) {
      return 1;
    } else {
      // in this case, second name is null
      return -1;
    }
    return firstName.compareTo(secondName);
  }

  /**
   * Sorts the lists of epochs by SNCL identifier, then start time increasing
   */
  private static class EpochComparator implements Comparator<EpochIdentifier> {

    static final EpochComparator instance = new EpochComparator();

    private EpochComparator() {
    }

    @Override
    public int compare(EpochIdentifier o1, EpochIdentifier o2) {
      if (o1 == null && o2 == null) {
        return 0;
      }
      else if (o1 == null) {
        return 1;
      }
      else if (o2 == null) {
        return -1;
      }

      int check = compareNamesPossiblyMissingComponents(o1.channelIdentifier, o2.channelIdentifier);
      if (check != 0) {
        return check;
      }
      return o1.startInstant.compareTo(o2.startInstant);
    }
  }
}
