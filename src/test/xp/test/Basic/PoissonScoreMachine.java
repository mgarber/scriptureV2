
package xp.test.Basic;

import org.apache.log4j.Logger;
import org.apache.commons.math3.distribution.PoissonDistribution;
/**
 *  Created on 2013-8-5  
 */

public class PoissonScoreMachine implements ScoreMachine {
	private double MAXSCORE=300;
	private double threshold=8.0;
	private static final Logger logger = Logger
			.getLogger(PoissonScoreMachine.class);
	public double getScore(double x, double lambda)
	{
		
		return getScore((int)Math.round(x),lambda);
	
	}
	public double getScore(int x, double lambda)
	{
		PoissonDistribution dist = new PoissonDistribution(lambda);
		double p = 1.0 - dist.cumulativeProbability(x);
		double score=-Math.log10(p);
    	if( score < MAXSCORE)
    		return score;
    	else
    		return MAXSCORE;
	}
	
	
	@Override
	public double getScore(double[] x) {
		// TODO Auto-generated method stub
	    	
		return getScore(x[0],x[1]);// assume that x[0] is x and x[1] is lambda
	}

	@Override
	public double getScore(double[] x, double[] background) {
		// TODO Auto-generated method stub
		return getScore(x[0],background[0]);
	}

	@Override
	public double getScore(double[] x, double[] globalBackground,
			double[] localBackGround) {
		// TODO Auto-generated method stub
		return Math.min(getScore(x,globalBackground), getScore(x,localBackGround));
	}
	public double getScore(double x, double global, double local)
	{
		return Math.min(getScore(x,global), getScore(x,local));
	}

	@Override
	public double getThreshold() {
		// TODO Auto-generated method stub
		return threshold;
	}

	@Override
	public void setThreshold(double t) {
		// TODO Auto-generated method stub
       threshold=t;
	}

}
