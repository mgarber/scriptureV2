package xp.test.command;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;

import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.Usage;
import net.sf.samtools.util.SortingCollection;
import nextgen.core.alignment.Alignment;
import nextgen.core.coordinatesystem.GenomicSpace;
import broad.core.annotation.ShortBED;
import broad.core.math.Distribution;

import org.apache.log4j.Level;
import org.apache.log4j.Logger;

import xp.test.Basic.BedGraphMultiScore;
import xp.test.Basic.FocusAndEnv;
import xp.test.Converter.BamIteratorFactory;
import xp.test.Converter.BedGraphMultiScoreReader;
import xp.test.Converter.LocalEnvReader;
import xp.test.DBI.ShortBEDTabixDBI;
import xp.test.Utils.JieCode;
import xp.test.Utils.JieCodeSortingCollection;

/**
 *  Created on 2013-3-11
 *  
 *    update to bam readers.
 *    
 */


public class CallPeak3 extends CommandLineProgram{
	static Logger logger = Logger.getLogger(CallPeak3.class.getName());	
	
	private static final String PROGRAM_VERSION = "0.01";
	private static int READS_THRESHOLD=5;
	private static int MAXGAP=200;
    private static int NEARBY=50000;
    private static int MAXSCORE=300;
    private static int MINPEAKSCORE=5;
	@Usage

    @Option(doc="This is caseA bam", shortName="A") public File CASEA;
    @Option(doc="This is caseB bam", shortName="B") public File CASEB;
    @Option(doc="This is control A bam", shortName="X") public File CONTROLA;
    @Option(doc="This is control B bam", shortName="Y") public File CONTROLB;
    @Option(doc="Chromosome Size",shortName="G") public String CHROMSIZES;
    @Option(doc="output file",shortName="o") public String FOUT;
    @Option(doc="paired end",shortName="p") public boolean ISPAIRED;
	
    public static void main(String[] argv){
    	System.exit(new CallPeak3().instanceMain(argv));
        
    }

	@Override
	protected int doWork() {
		// TODO Auto-generated method stub
		 logger.setLevel(Level.DEBUG);
    	 logger.info("initialize genomic space");
    	 GenomicSpace GENOMESPACE = new GenomicSpace(CHROMSIZES);
    	 logger.info("initialize genomic space done");
    	
    	 logger.info(GENOMESPACE.getReferenceNames());
    	 try {
    	 String isPaired="s";
    	 if(ISPAIRED) isPaired="p";
    	 
    	 Iterator<Alignment> caseAIter = BamIteratorFactory.makeIterator(CASEA,isPaired); 
    	 Iterator<Alignment> caseBIter = BamIteratorFactory.makeIterator(CASEB,isPaired); 
    	 Iterator<Alignment> controlAIter = BamIteratorFactory.makeIterator(CONTROLA,isPaired); 
    	 Iterator<Alignment> controlBIter = BamIteratorFactory.makeIterator(CONTROLB,isPaired); 
    	
    	
    	 JieCodeSortingCollection jie = new JieCodeSortingCollection(GENOMESPACE);
    	 jie.add(caseAIter,0);
    	 logger.info("Coverage CaseA "+jie.getCoverage(0));
    	 logger.info("Reads Number CaseA "+jie.getReadsNumber(0));
    	 double globalLambdaCaseA=(double)jie.getCoverage(0)/GENOMESPACE.getLength();
    	 long caseACoverage=jie.getCoverage(0);
    	 
    	 
    	 
    	 jie.add(caseBIter,1);
    	 logger.info("Coverage CaseB "+jie.getCoverage(1));
    	 logger.info("Reads Number CaseB "+jie.getReadsNumber(1));
    	 double globalLambdaCaseB=(double)jie.getCoverage(1)/GENOMESPACE.getLength();
    	 long caseBCoverage=jie.getCoverage(1);
    	 
    	 jie.add(controlAIter,2);
    	 logger.info("Coverage ControlA "+jie.getCoverage(2));
    	 logger.info("Reads Number ControlA "+jie.getReadsNumber(2));
    	 double globalLambdaControlA=(double)jie.getCoverage(2)/GENOMESPACE.getLength();
    	 long controlACoverage=jie.getCoverage(2);
    	 
   	 
    	 jie.add(controlBIter,3);
    	 logger.info("Coverage ControlB "+jie.getCoverage(3));
    	 logger.info("Reads Number ControlB "+jie.getReadsNumber(3));
    	 double globalLambdaControlB=(double)jie.getCoverage(3)/GENOMESPACE.getLength();
    	 long controlBCoverage=jie.getCoverage(3);
	   	 
    	 double controlAToCaseAIndex=Double.valueOf(caseACoverage)/Double.valueOf(controlACoverage);
    	 double controlBToCaseBIndex=Double.valueOf(caseBCoverage)/Double.valueOf(controlBCoverage);
         
    	 FileWriter fout = new FileWriter(FOUT);
    	 FileWriter fpeak = new FileWriter(FOUT.toString()+".peak.bed");
    	 BufferedWriter bf = new BufferedWriter(fout);
    	 
    	 BedGraphMultiScoreReader c = new BedGraphMultiScoreReader(jie);
    	 LocalEnvReader b = new LocalEnvReader(c);
    	 int k=0;
    	 logger.debug("Has Next "+b.hasNext());
    	 
    	 
    	 int state=0;
    	 int lastTid=0;
    	 int lastStart=0;
    	 int lastEnd=0;
    	 double maxScore=0;
    	 int peakIndex=0;
    	 BedGraphMultiScore bufferBed=new BedGraphMultiScore(0,0,0,new int[1]);
    	 
    	 while(b.hasNext())
    	 {
    		 FocusAndEnv bgraph = b.next();
    		 
    		 k+=1;
    		 if  (k%100000==0)
    			 logger.info(String.format("%d ",k)+bgraph);
    	 
    	     BedGraphMultiScore focus=bgraph.getFocus();
    	     BedGraphMultiScore env=bgraph.getEnvs();
    	    
    	     double lambdaA = (double)env.getScore()[2]/env.length()*controlAToCaseAIndex;
    	     double lambdaB = (double)env.getScore()[3]/env.length()*controlBToCaseBIndex;
    	     double fragmentLambdaA=(double)focus.getScore()[2]*controlAToCaseAIndex;
    	     double fragmentLambdaB=(double)focus.getScore()[3]*controlAToCaseAIndex;
    	     
    	     lambdaA=Math.max(lambdaA, fragmentLambdaA);    	     
			 lambdaB=Math.max(lambdaB, fragmentLambdaB);    	     
    	     int caseAScore=focus.getScore()[0];
    	     int caseBScore=focus.getScore()[1];
    	     if ( lambdaA<globalLambdaCaseA ) lambdaA=globalLambdaCaseA;
    	     if ( lambdaB<globalLambdaCaseB ) lambdaB=globalLambdaCaseB;
    	     double score=getSkellamScore(lambdaB,lambdaA,caseBScore,caseAScore);
    	    if(Math.abs(score)>4)
    	    {
    	    	
    	    	bf.write(jie.getTid2chr().get(focus.getTid()));
    	    	bf.write("\t");
    	    	bf.write(String.format("%d",focus.getStart()));
    	    	bf.write("\t");
    	    	bf.write(String.format("%d",focus.getEnd()));
    	    	bf.write("\t");
    	    	bf.write(String.format("%f",score));
    	    	bf.write("\t");
    	    	bf.write(String.format("[  %d / %.3f  vs  %d / %.3f ]",caseAScore,lambdaA,caseBScore,lambdaB));
    	    	bf.write("\n");
    	    	
    	    	if (focus.getTid()!=lastTid || focus.getStart() - lastEnd > MAXGAP)
    	    	{
    	    		if(Math.abs(maxScore)>MINPEAKSCORE)
    	    		{
    	    			peakIndex+=1;
    	    			fpeak.write(String.format("%s\t%d\t%d\tpeak_%d\t%.2f\t%d\t%d\n", jie.getTid2chr().get(lastTid), lastStart, lastEnd , peakIndex, maxScore,bufferBed.getStart(),bufferBed.getEnd()));
    	    		    
    	    		}
    	    		lastStart=focus.getStart();
	    		    lastEnd=focus.getEnd();
	    		    lastTid=focus.getTid();
	    		    bufferBed=focus;
	    		    maxScore=score;
    	    	}
    	    	else
    	    	{
    	    		lastEnd=focus.getEnd();
    	    		if (Math.abs(score) > Math.abs(maxScore))
    	    		{
    	    			maxScore=score;
    	    			bufferBed=focus;
    	    		}
    	    	}
    	    }
    	 }
    	 
    	 
    	 if(Math.abs(maxScore)>MINPEAKSCORE)
 		{
 			peakIndex+=1;
 			fpeak.write(String.format("%s\t%d\t%d\tpeak_%d\t%.2f\t%d\t%d\n", jie.getTid2chr().get(lastTid), lastStart, lastEnd , peakIndex, maxScore,bufferBed.getStart(),bufferBed.getEnd()));
 		   
 		} 
    	 fpeak.close();
    	bf.close();
    	fout.close();
    	 }
    	 catch (IOException e)
     	{
     	//To do	
    		logger.error(e);
     	}
    	 return 0;
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
