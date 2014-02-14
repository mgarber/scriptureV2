package xp.core.Basic;

import broad.core.math.ScanStatistics;

public class PoissonScanStatisticScoreMachine implements ScoreMachine {
	private double effectiveGenomeSize;
	private double window;
	
	private double MAXSCORE=300;

	public PoissonScanStatisticScoreMachine(double window, double effectiveGenomeSize) {
		this.window = window;
		this.effectiveGenomeSize = effectiveGenomeSize;
	}
	
	public double getScore(double x, double lambda)
	{
		
		return getScore((int)Math.round(x),lambda);
	
	}
	
	public double getScore(int x, double lambda)
	{
		double p = ScanStatistics.calculatePVal(x, lambda, window, effectiveGenomeSize);
		double score=-Math.log10(p);
		
		return score < MAXSCORE ? score : MAXSCORE;
	}
	
	

	public double getScore(double[] x) {
		return getScore(x[0], x[1]);// assume that x[0] is x and x[1] is lambda
	}

	@Override
	public double getScore(double[] x, double[] background) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public double getScore(double[] x, double[] globalBackground,
			double[] localBackGround) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public double getThreshold() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public void setThreshold(double t) {
		// TODO Auto-generated method stub

	}

}
