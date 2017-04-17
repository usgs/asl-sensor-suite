package asl.sensor.test;

import java.util.Arrays;

import org.apache.commons.math3.fitting.leastsquares.LeastSquaresBuilder;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresOptimizer;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresProblem;
import org.apache.commons.math3.fitting.leastsquares.LevenbergMarquardtOptimizer;
import org.apache.commons.math3.fitting.leastsquares.MultivariateJacobianFunction;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;
import org.junit.Test;

public class SolverTest {

  @Test
  public void solverWorks() {

    MultivariateJacobianFunction jacobian = new MultivariateJacobianFunction() {
      
      private double[] doCalculation(double[] point) {
        
        double[] result = new double[point.length];
        
        for (int i = 0; i < result.length; ++i) {
          double first = point[i] - 3;
          double second = point[i] + 4;
          result[i] = first * second;
        }
        
        return result;
        
      }
      
      public Pair<RealVector, RealMatrix> value(final RealVector point) {
        
        double[] pointArr = point.toArray();
        double[] resultArr = doCalculation(pointArr);
        
        RealVector result = MatrixUtils.createRealVector(resultArr);
        
        // System.out.println(init + "," + res);
        
        double[][] jacobianArr = new double[resultArr.length][pointArr.length];

        for (int i = 0; i < pointArr.length; ++i) {
          double change = 1E-5;
          double initVal = pointArr[i];
          pointArr[i] += change;
          double[] forwardDiff = doCalculation(pointArr);
          
          for (int j = 0; j < forwardDiff.length; ++j) {
            jacobianArr[j][i] = (forwardDiff[j] - resultArr[j]) / change;
          }
          
          // set value back to what it initially was
          pointArr[i] = initVal;
          
        }
        
        RealMatrix jacobian = MatrixUtils.createRealMatrix(jacobianArr);
        
        return new Pair<RealVector, RealMatrix>(result, jacobian);        
      }
    };
    
    RealVector initialGuess1 = MatrixUtils.createRealVector(new double[]{-10});
    RealVector initialGuess2 = MatrixUtils.createRealVector(new double[]{20});
    
    RealVector obsResVector = MatrixUtils.createRealVector(new double[]{0});
    
    LeastSquaresProblem lsp = new LeastSquaresBuilder().
        start(initialGuess1).
        target(obsResVector).
        model(jacobian).
        lazyEvaluation(false).
        maxEvaluations(Integer.MAX_VALUE).
        maxIterations(Integer.MAX_VALUE).
        build();
    
    LeastSquaresOptimizer optimizer = new LevenbergMarquardtOptimizer();
   
    LeastSquaresOptimizer.Optimum optimum = optimizer.optimize(lsp);
    
    System.out.println( Arrays.toString( optimum.getPoint().toArray() ) );
    
    lsp = new LeastSquaresBuilder().
        start(initialGuess2).
        target(obsResVector).
        model(jacobian).
        lazyEvaluation(false).
        maxEvaluations(Integer.MAX_VALUE).
        maxIterations(Integer.MAX_VALUE).
        build();
    
    optimum = optimizer.optimize(lsp);
    
    System.out.println( Arrays.toString( optimum.getPoint().toArray() ) );
    
  }
  
}
