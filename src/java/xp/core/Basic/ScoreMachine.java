
package xp.core.Basic;
/**
 *  Created on 2013-4-1  
 */

public interface ScoreMachine {
  public double getScore(double[] x);
  public double getScore(double[] x, double[] background);
  public double getScore(double[] x, double[] globalBackground, double[] localBackGround);
  public double getThreshold();
  public void setThreshold(double t);
}
