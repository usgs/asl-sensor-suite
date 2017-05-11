package asl.sensor.test;

import static org.junit.Assert.*;

import org.apache.commons.math3.fitting.leastsquares.LeastSquaresBuilder;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresOptimizer;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresProblem;
import 
org.apache.commons.math3.fitting.leastsquares.LevenbergMarquardtOptimizer;
import 
org.apache.commons.math3.fitting.leastsquares.MultivariateJacobianFunction;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;
import org.junit.Test;

public class SolverTest {

  @Test
  public void solverWorks() {

    MultivariateJacobianFunction jacobian = new MultivariateJacobianFunction() {
      
      // evaluate the result at a specific point
      private double[] doCalculation(double[] point) {
        
        double[] result = new double[point.length];
        
        for (int i = 0; i < result.length; ++i) {
          double first = point[i] - 3;
          double second = point[i] + 4;
          result[i] = first * second;
        }
        
        return result;
        
      }
      
      // computes value by direct evaluation, jacobian by forward difference
      public Pair<RealVector, RealMatrix> value(final RealVector point) {
        
        double[] pointArr = point.toArray();
        double[] resultArr = doCalculation(pointArr);
        
        RealVector result = MatrixUtils.createRealVector(resultArr);
        
        // System.out.println(init + "," + res);
        
        double[][] jacobianArr = new double[resultArr.length][pointArr.length];

        for (int i = 0; i < pointArr.length; ++i) {
          // create a new array, change the relevant variable, return
          double change = 1E-5;
          double[] fwdDiffArr = point.toArray();

          fwdDiffArr[i] += change;
          // fwdDiffArr[1] += change; test should fail if uncommented
          
          for (int k = 0; k < pointArr.length; ++k) {
            if (i != k) {
              assertEquals(fwdDiffArr[k], pointArr[k], change);
            }
          }
          
          // System.out.println( Arrays.toString(fwdDiffArr) );
          double[] forwardDiffRes = doCalculation(fwdDiffArr);
          
          for (int j = 0; j < forwardDiffRes.length; ++j) {
            jacobianArr[j][i] = (forwardDiffRes[j] - resultArr[j]) / change;
          }

        }
        
        RealMatrix jacobian = MatrixUtils.createRealMatrix(jacobianArr);
        
        return new Pair<RealVector, RealMatrix>(result, jacobian);        
      }
    };
    
    RealVector initialGuess1 = 
        MatrixUtils.createRealVector(new double[]{-10, -10});
    RealVector initialGuess2 = 
        MatrixUtils.createRealVector(new double[]{20, 20});
    
    RealVector obsResVector = 
        MatrixUtils.createRealVector(new double[]{0, 0});
    
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
    
    double[] firstOptimum = optimum.getPoint().toArray();
    
    // System.out.println( Arrays.toString( optimum.getPoint().toArray() ) );
    
    lsp = new LeastSquaresBuilder().
        start(initialGuess2).
        target(obsResVector).
        model(jacobian).
        lazyEvaluation(false).
        maxEvaluations(Integer.MAX_VALUE).
        maxIterations(Integer.MAX_VALUE).
        build();
    
    optimum = optimizer.optimize(lsp);
    
    double secondOptimum[] = optimum.getPoint().toArray();
    
    // System.out.println( Arrays.toString( optimum.getPoint().toArray() ) );
    
    assertEquals(firstOptimum[0], -4, 0.1);
    assertEquals(secondOptimum[0], 3, 0.1);
    
    for (int i = 0; i < firstOptimum.length; ++i) {
      assertEquals(firstOptimum[0], firstOptimum[i], 0.1);
      assertEquals(secondOptimum[0], secondOptimum[i], 0.1);
    }
  }
  
}
