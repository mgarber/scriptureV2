
package xp.test.Basic;

import org.apache.log4j.Logger;

import broad.core.math.Distribution;

/**
 *  Created on 2013-4-2  
 */

public class SkellamScoreMachine implements ScoreMachine {
	private double threshold = 8.0;
	private static double MAXSCORE=30.0;
	private static final Logger logger = Logger
			.getLogger(SkellamScoreMachine.class);

	@Override
	public double getScore(double[] x) {
		// TODO Auto-generated method stub
		// NOT IMPLEMENT
		return 0;
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
	    int sampleNumber=x.length;
	    if(sampleNumber==2)
	    {
	    	return getScore2(x,globalBackground,localBackGround);
	    	
	    }
	    else if (sampleNumber==4)
	    {
	    	return getScore4(x,globalBackground,localBackGround);
	    }
	    else
	    {
	    	logger.warn("Wrong Sample Number, not implement yet");
	    }
		return 0;
	}

	@Override
	public double getThreshold() {
		 return threshold;
	}

	@Override
	public void setThreshold(double t) {
	     this.threshold=t;	
	}
	
	
	private double getScore2(double[] x, double[] globalBackground,double[] localBackground)
	{
		/*
		 * ORDER:
		 * Case , Control
		 *   x[0] is Case SignalCount
		 *   x[1] is Control backgroundCount
		 *   
		 *   globalBackground[0] is Case global background lambda
		 *   globalBackground[1] is Control global background lambda
		 *   
		 *   localBackground[0] is Case local background lambda
		 *   localBackground[1] is Control local background lambda
		 *   
		 */
		Double globalScore;
		Double localScore;
		globalScore=getSkellamScore((int)Math.round(x[0]),globalBackground[0],(int)Math.round(x[1]),globalBackground[1]);
		localScore=getSkellamScore((int)Math.round(x[0]),localBackground[0],(int)Math.round(x[1]),localBackground[1]);
		double score;
		if(localScore==null)
			score=globalScore;
		else
			score=Math.min(globalScore, localScore);
		return score;
	}
	private double getScore4(double[] x,double[] globalBackground,double[] localBackground)
	{
		/*
		 *  order:
		 *  CaseA, ControlA, CaseB, ControlB
		 *  
		 *  x[0] caseA count
		 *  x[1] controlA count
		 *  x[2] caseB count
		 *  x[3] conttrolB count
		 *  
		 *  globalBackground[0]  caseA global lambda
		 *  globalBackground[1]  controlA global lambda
		 *  globalBackground[2]  caseB global lambda
		 *  globalBackground[3]  controlB global lambda
		 *  
		 *  localBackground[0]   caseA local lambda
		 *  localBackground[1]   controlA local lambda
		 *  localBackground[2]   caseB local lambda
		 *  localBackground[3]   controlB local lambda
		 *  
		 *  
		 *  Define Score is Max Score of  
		 *          caseA vs controlA
		 *          caseB vs controlB
		 *          caseA + caseB  vs  controlA + controlB
		 */
        
		double[] tX= new double[2];
		double[] tGlobalBackground = new double[2];
		double[] tLocalBackground = new double[2];
		
		double scoreA;
		double scoreB;
		double scoreAB;
		
		// CaseA vs Control A
		tX[0]=x[0];
		tX[1]=x[1];
		tGlobalBackground[0]=globalBackground[0];
		tGlobalBackground[1]=globalBackground[1];
        tLocalBackground[0]=localBackground[0];
        tLocalBackground[1]=localBackground[1];
		scoreA=getScore2(tX,tGlobalBackground,tLocalBackground);
		// CaseB vs Control B
		tX[0]=x[2];
		tX[1]=x[3];
		tGlobalBackground[0]=globalBackground[2];
		tGlobalBackground[1]=globalBackground[3];
        tLocalBackground[0]=localBackground[2];
        tLocalBackground[1]=localBackground[3];
		scoreB=getScore2(tX,tGlobalBackground,tLocalBackground);
			
		// CaseA + CaseB vs Control A + Control B
		tX[0]+=x[0];
		tX[1]+=x[1];
		tGlobalBackground[0]+=globalBackground[0];
		tGlobalBackground[1]+=globalBackground[1];
        tLocalBackground[0]+=localBackground[0];
        tLocalBackground[1]+=localBackground[1];
		scoreAB=getScore2(tX,tGlobalBackground,tLocalBackground);
		
		double score;
		score=Math.max(scoreA, scoreB);
		score=Math.max(score, scoreAB);
		return score;
		
	}
	
	
	private static Double getSkellamScore(int signalCount, double signalLambda, int controlCount, double controlLambda)  {
		int k=signalCount-controlCount;
		if(controlLambda==0 || signalLambda==0)
		{
			return null;
		}
		
		double rightPvalue = Distribution.skellamRightTail(k, signalLambda, controlLambda);
		double score=-Math.log10(rightPvalue);
	    	if( score < MAXSCORE)
	    		return score;
	    	else
	    		return MAXSCORE;
	    
	  }
	

}
