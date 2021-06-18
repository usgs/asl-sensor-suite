package asl.sensor.input;

import static asl.sensor.test.TestUtils.TEST_DATA_LOCATION;
import static org.junit.Assert.assertArrayEquals;

import asl.sensor.input.NoiseModel.NoiseModelFormatException;
import java.io.IOException;
import org.junit.Test;

public class NoiseModelTest {

  @Test
  public void checkNoiseModelValues() throws IOException, NoiseModelFormatException {
    double[] xValues = {2, 2.18, 2.38, 2.59, 2.83, 3.08, 3.36, 3.67, 4, 4.36, 4.76, 5.19, 5.66,
        6.17, 6.73, 7.34, 8, 8.72, 9.51, 10.37, 11.31, 12.34, 13.45, 14.67, 16, 17.45, 19.03,
        20.75, 22.63, 24.68, 26.91, 29.34, 32, 34.9, 38.05, 41.5, 45.25, 49.35, 53.82, 58.69, 64,
        69.79, 76.11, 83, 90.51, 98.7, 107.63, 117.38, 128, 139.58, 152.22, 166, 181.02, 197.4,
        215.27, 234.75, 256, 279.17, 304.44, 331.99, 362.04, 394.81, 430.54, 469.51, 512, 558.34,
        608.87, 663.98, 724.08, 789.61, 861.08, 939.01, 1024};
    double[] yValues = {-134.1, -138.6, -142.74, -141.53, -140.4, -138.64, -136.84, -135.19,
        -133.58, -132.68, -132.54, -132.64, -132.62, -132.5, -132.63, -133.29, -135.61, -140.7,
        -146.26, -150.29, -152.11, -152.18, -151.58, -151.73, -152.57, -154.58, -157.4, -160.86,
        -164.14, -166.71, -168.99, -170.3, -171.39, -172.35, -173.05, -173.87, -174.36, -174.76,
        -175.33, -175.52, -175.47, -175.13, -174.85, -174.53, -174.22, -173.88, -173.64, -173.56,
        -173.55, -173.65, -173.65, -173.55, -173.55, -173.31, -173.31, -172.4, -172.4, -50.5, -50.5,
        -170.44, -170.44, -50.5, -50.5, -166.98, -166.98, -50.5, -50.5, -50.5, -50.5, -50.5, -50.5,
        -161.98, -161.98};

    NoiseModel model =
        new NoiseModel(TEST_DATA_LOCATION + "station_noise_models/IU.ANMO.00.LH1.csv");
    double[][] modelValues = model.getPeriodAndMean();

    assertArrayEquals(xValues, modelValues[0], 0.);
    assertArrayEquals(yValues, modelValues[1], 0.);
  }
}
