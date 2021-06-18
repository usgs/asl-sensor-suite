package asl.sensor.experiment;

import static asl.sensor.test.TestUtils.TEST_DATA_LOCATION;
import static org.junit.Assert.assertEquals;

import asl.sensor.input.NoiseModel;
import asl.sensor.input.NoiseModel.NoiseModelFormatException;
import java.io.IOException;
import org.jfree.data.xy.XYSeries;
import org.junit.Test;

public class NoiseModelExperimentTest {


  @Test
  public void verifyModelPlottableIsCorrect() throws IOException, NoiseModelFormatException {
    String noiseModelFile = TEST_DATA_LOCATION + "station_noise_models/IU.ANMO.00.LH1.csv";
    NoiseModel model = new NoiseModel(noiseModelFile);
    double[][] modelValues = model.getPeriodAndMean();

    NoiseModelExperiment exp = new NoiseModelExperiment();
    exp.loadNoiseModel(noiseModelFile);
    XYSeries plottable = exp.getPlottableNoiseModelData(false);

    for (int i = 0; i < plottable.getItemCount(); ++i) {
      assertEquals(modelValues[0][i], plottable.getDataItem(i).getXValue(), 0.);
      assertEquals(modelValues[1][i], plottable.getDataItem(i).getYValue(), 0.);
    }
  }
}
