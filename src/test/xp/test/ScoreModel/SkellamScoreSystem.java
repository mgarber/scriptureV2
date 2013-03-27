package xp.test.ScoreModel;

import java.io.IOException;
import java.util.List;

import org.apache.log4j.Level;
import org.apache.log4j.Logger;

import xp.test.ChIPSeqPeakCaller;
import xp.test.DBI.AlignmentDBI;

import broad.core.math.Distribution;

import nextgen.core.annotation.Annotation;

/**
 *  Created on 2013-3-5  
 * @author zhuxp
 *
 */
public class SkellamScoreSystem implements ScoreSystem {

	/* (non-Javadoc)
	 * @see xp.test.ScoreSystem#getScore(nextgen.core.annotation.Annotation)
	 */
	static Logger logger = Logger.getLogger(ChIPSeqPeakCaller.class.getName());
    static double MAXSCORE=300.0; //if pvalue is 0.0	
	
	private AlignmentDBI caseA;
	private AlignmentDBI caseB;
	private AlignmentDBI controlA;
	private AlignmentDBI controlB;
	
	private boolean hasControlA = false;
	private boolean hasControlB = false;
	
	private boolean hasGlobalCaseA = false;
	private double globalCaseA;
	
	private boolean hasGlobalCaseB = false;
	private double globalCaseB;
	
	private boolean hasGlobalControlA = false;
	private double  globalControlA;
	
	private boolean hasGlobalControlB = false;
	private double globalControlB;
	
	private int nearby = 50000; // 10K  +- 5K
	private boolean isInit = false;
	
	private double globalLambdaAIndex;
	private double globalLambdaBIndex;
	
	
    public SkellamScoreSystem(AlignmentDBI caseA, AlignmentDBI caseB, AlignmentDBI controlA, AlignmentDBI controlB)
    {
    	 this.caseA=caseA;
    	 this.caseB=caseB;
    	 this.controlA=controlA;
    	 this.controlB=controlB;
    	 this.hasControlA=true;
    	 this.hasControlB=true;
    	 initGlobalCount();
    	 //logger.setLevel(Level.DEBUG);
    	
    }
    
    public SkellamScoreSystem(AlignmentDBI caseA, AlignmentDBI caseB, AlignmentDBI controlA, AlignmentDBI controlB, int nearby)

    {
    	this(caseA,caseB,controlA,controlB);
    	setNearby(nearby);
    }
    public SkellamScoreSystem(AlignmentDBI caseA, AlignmentDBI caseB)
    {
    	 this.caseA=caseA;
    	 this.caseB=caseB;
    	 initGlobalCount();
    	
    }
	/**
	 * BEANS
	 */
	
	public AlignmentDBI getCaseA() {
		return caseA;
	}

	public void setCaseA(AlignmentDBI caseA) {
		this.caseA = caseA;
	}

	public AlignmentDBI getCaseB() {
		return caseB;
	}

	public void setCaseB(AlignmentDBI caseB) {
		this.caseB = caseB;
	}

	public AlignmentDBI getControlA() {
		return controlA;
	}

	public void setControlA(AlignmentDBI controlA) {
		this.controlA = controlA;
	}

	public AlignmentDBI getControlB() {
		return controlB;
	}

	public void setControlB(AlignmentDBI controlB) {
		this.controlB = controlB;
	}

	public int getNearby() {
		return nearby;
	}

	public void setNearby(int nearby) {
		this.nearby = nearby;
	}

	
	
	
	public void initGlobalCount()
	{
		if (!hasGlobalCaseA) {globalCaseA=caseA.getGlobalCount(); hasGlobalCaseA=true;}
		if (!hasGlobalCaseB) {globalCaseB=caseB.getGlobalCount(); hasGlobalCaseB=true;}
		if (hasControlA & !hasGlobalControlA) {globalControlA=controlA.getGlobalCount(); hasGlobalControlA=true;}
		if (hasControlB & !hasGlobalControlB) {globalControlB=controlB.getGlobalCount(); hasGlobalControlB=true;}
		long length=0;
		List<String> chrs=caseA.getChrs();
		for(String chr: chrs)
		{length+=caseA.getChrLength(chr);};
		globalLambdaAIndex=globalCaseA/length;
	    globalLambdaBIndex=globalCaseB/length;
		isInit=true;
	}
	
	@Override
	public Double getScore(Annotation a) {
		return getScore(a.getChr(),a.getStart(),a.getEnd());
	}

	/* (non-Javadoc)
	 * @see xp.test.ScoreSystem#getScore(java.lang.String, int, int)
	 */
	@Override
	public Double getScore(String chr, int start, int end) {
		double score;
		logger.debug(chr); 
		logger.debug(start); 
		logger.debug(end); 
		score = _getScore(chr,start,end);
		 
		 return score;
	}
	
	
	private Double _getScore(String chr,int start, int end)
	{
	   
	   int chrLen=(int)caseA.getChrLength(chr);
	   
	   if (end>chrLen) end=chrLen;
	   if (start<0)  start=0;
	   int caseLen=end-start;
	   
	   if (caseLen < 0)
	   {
		   logger.warn("query position is wrong ( start > end )");
		   return null;
	   }
	   
	   
	   
	   double countCaseA=caseA.getLocalCount(chr, start, end);
	   double countCaseB=caseB.getLocalCount(chr, start, end);
	   
	   
	   int controlStart=start-nearby;
	   int controlEnd=end+nearby;
	   
	   if(controlEnd > chrLen) controlEnd=chrLen;
	   if(controlStart < 0) controlStart=0;
	   
	   
	   int controlLen=controlEnd-controlStart; 
	  double countControlA;
	  double countControlB;
	  double lambdaA;
	  double lambdaB;
	  
	   if (hasControlA & hasControlB)
	   {
	   countControlA=controlA.getLocalCount(chr,controlStart,controlEnd);
	   countControlB=controlB.getLocalCount(chr,controlStart,controlEnd);
	   lambdaA = (globalCaseA/globalControlA) * countControlA/controlLen * caseLen;
	   lambdaB = (globalCaseB/globalControlB) * countControlB/controlLen * caseLen;
	   }
	   else
	   {
	   countControlA=caseA.getLocalCount(chr,controlStart,controlEnd);
	   countControlB=caseB.getLocalCount(chr,controlStart,controlEnd);
	   lambdaA = countControlA/controlLen * caseLen;
	   lambdaB = countControlB/controlLen * caseLen;
	   }
	   
	   logger.debug("lambdaA"+lambdaA);
	   logger.debug("lambdaB"+lambdaB);
	   if(lambdaA==0)
	   {
		   // TO DOO.
		  lambdaA=globalLambdaAIndex*caseLen;
		
	   }
	   if(lambdaB==0)
	   {
		   lambdaB=globalLambdaBIndex*caseLen;
	   }
	   if(countCaseA==0 & countCaseB==0)
	   {
		   return 0.0; //no need to calculate score if there's no reads.
	   }
	   double score;
		score=getSkellamScore(lambdaA,lambdaB,(int)Math.round(countCaseA),(int)Math.round(countCaseB));
        	  
		
		
		return score;	
	}
	
	
	protected static double getSkellamScore(double backgroundLambda, double signalLambda, int backgroundCount, int signalCount)  {
		int k=signalCount-backgroundCount;
	    	
		double leftPvalue=Distribution.skellamLeftTail(k, signalLambda, backgroundLambda);
		
		double rightPvalue = Distribution.skellamRightTail(k, signalLambda, backgroundLambda);
		
	    if (rightPvalue < leftPvalue)
	    {
	    	if(rightPvalue!=0)
	    		return -Math.log10(rightPvalue);
	    	else
	    		return MAXSCORE;
	    }
	    else
	    {
	    	if(leftPvalue!=0)
	    		return Math.log10(leftPvalue);
	    	else
	    		return -MAXSCORE;
	    }
	    
	  }
	}
	
	


