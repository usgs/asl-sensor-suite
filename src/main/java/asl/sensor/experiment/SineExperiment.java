package asl.sensor.experiment;

import java.util.ArrayList;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import asl.sensor.input.DataStore;
import asl.sensor.utils.TimeSeriesUtils;

public class SineExperiment extends Experiment {

  private double calSDev, outSDev, peakPeakFreq;

  public SineExperiment() {
    super();
    calSDev = 0.;
    outSDev = 0.;
    peakPeakFreq = 0.;
  }

  public double getCalAmplitude() {
    return calSDev;
  }

  public double getOutAmplitude() {
    return outSDev;
  }

  public double getEstSineFreq() {
    return peakPeakFreq;
  }

  @Override
  protected void backend(DataStore ds) {
    calSDev = 0.;
    outSDev = 0.;
    peakPeakFreq = 0.;

    double[] calTimeSeries = ds.getBlock(0).getData().clone();
    double[] outTimeSeries = ds.getBlock(1).getData().clone();

    dataNames.add(ds.getBlock(0).getName());
    dataNames.add(ds.getBlock(1).getName());


    // get the sine wave frequency by measuring wavelengths (get distance between peaks)
    int currentPeakDistance = 0;
    int totalPeakDistance = 0;
    int peakCount = -1; // start at negative 1 so first peak counts as 0
    for (int i = 1; i < outTimeSeries.length - 1; ++i) {
      if (outTimeSeries[i] > outTimeSeries[i - 1] && outTimeSeries[i] > outTimeSeries[i + 1]) {
        ++peakCount;
        totalPeakDistance += currentPeakDistance;
        currentPeakDistance = 0;
      } else {
        ++currentPeakDistance;
      }
    }
    if (peakCount > 0) {
      peakPeakFreq = totalPeakDistance / (double) peakCount;
    } else {
      peakPeakFreq = 0;
    }

    // we would like to have the signals centered on 0 for lining them up
    calTimeSeries = TimeSeriesUtils.demean(calTimeSeries);
    outTimeSeries = TimeSeriesUtils.demean(outTimeSeries);
    // standard deviation is a good
    calSDev = new DescriptiveStatistics(calTimeSeries).getStandardDeviation();
    outSDev = new DescriptiveStatistics(outTimeSeries).getStandardDeviation();

    // add plots sine waves
    XYSeriesCollection xysc = new XYSeriesCollection();
    XYSeries cal = new XYSeries(ds.getBlock(0).getName() + " [cal]");
    XYSeries out = new XYSeries(ds.getBlock(1).getName() + " [out, scaled]");
    double interval = ds.getBlock(0).getInterval() / (double) TimeSeriesUtils.ONE_HZ_INTERVAL;
    for (int i = 0; i < calTimeSeries.length; ++i) {
      cal.add(interval * i, calTimeSeries[i]);
      out.add(interval * i, outTimeSeries[i] * calSDev / outSDev);
    }
    xysc.addSeries(cal);
    xysc.addSeries(out);
    xySeriesData = new ArrayList<XYSeriesCollection>();
    xySeriesData.add(xysc);
    // produce linearity plots
    xysc = new XYSeriesCollection();
    // booleans: don't autosort, do allow duplicate values
    XYSeries lin = new XYSeries(ds.getBlock(1).getName() + " linearity", false, true);
    for (int i = 0; i < calTimeSeries.length; ++i) {
      lin.add(calTimeSeries[i], outTimeSeries[i]);
    }
    xysc.addSeries(lin);
    xySeriesData.add(xysc);
  }

  @Override
  public int blocksNeeded() {
    return 2;
  }

  @Override
  public boolean hasEnoughData(DataStore ds) {
    for (int i = 0; i < blocksNeeded(); ++i) {
      if (!ds.blockIsSet(i)) {
        return false;
      }
    }
    return true;
  }


}
