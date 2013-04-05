package xp.test.command;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Comparator;
import java.util.Iterator;

import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.Usage;
import net.sf.samtools.util.BinaryCodec;
import net.sf.samtools.util.RuntimeEOFException;
import net.sf.samtools.util.SortingCollection;
import nextgen.core.alignment.Alignment;
import nextgen.core.annotation.Annotation;
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
 *    log          : update to bam readers.
 *    Log 2013-3-14: revise to minimize the score for three lambda. instead of using the biggest lambda for skellam distribution
 *                   revise to split + and - peaks
 *                   revise to speed up (IGNORE the low reads region)
 *                   revise to read more format? (at least shortbed index).
 *                   add qvalue table
 *                   set the correct max score for pvalue.
 *                   use blocks instead of start end ( revise in JieCodeSortCollection)
 *    To Do:  Extension 100bp ? not conflict with blocks?
 *        
 *            Control is enriched than Case?
 *            separate ScoreUSB to ExpDistribution.
 *            Windows Sliding ? ( EnvScore?  improve FocusEnv?)
 *                          
 *                          
 *    
 */

/**
 * For calculate the Qvalue
 * @author zhuxp
 *
 */
class ScoreUSB
{
	private double score;
	private int USB;
	public ScoreUSB(double score, int uSB) {
		super();
		this.score = score;
		USB = uSB;
	}
	public double getScore() {
		return score;
	}
	public void setScore(double score) {
		this.score = score;
	}
	public int getUSB() {
		return USB;
	}
	public void setUSB(int uSB) {
		USB = uSB;
	}
	
}
class ScoreUSBCodec implements SortingCollection.Codec<ScoreUSB>
{
	private final BinaryCodec binaryCodec = new BinaryCodec();
	@Override
	public ScoreUSB decode() {
		// TODO Auto-generated method stub
		double score;
		try{
			score=binaryCodec.readDouble();
		} catch ( RuntimeEOFException e)
		{
			return null;
		}
		int USB=binaryCodec.readInt();
		return new ScoreUSB(score,USB);
		
	}

	@Override
	public void encode(ScoreUSB s) {
		// TODO Auto-generated method stub
		binaryCodec.writeDouble(s.getScore());
		binaryCodec.writeInt(s.getUSB());
	}

	@Override
	public void setInputStream(InputStream is) {
		this.binaryCodec.setInputStream(is);
	}

	@Override
	public void setOutputStream(OutputStream os) {
		this.binaryCodec.setOutputStream(os);
		
	}
	public ScoreUSBCodec clone()
	{
		return new ScoreUSBCodec();
		
	}
	
	
}
class ScoreUSBComparator implements Comparator<ScoreUSB> //reverse max first
{

	@Override
	public int compare(ScoreUSB o1, ScoreUSB o2) {
		return Double.compare(o2.getScore(), o1.getScore());
	}
	
	
}




public class CallPeak5 extends CommandLineProgram{
	static Logger logger = Logger.getLogger(CallPeak5.class.getName());	
	
	private static final String PROGRAM_VERSION = "0.02";
	private static int READS_THRESHOLD=3; //less than this reads coverage in sum of case would skip 
	private static int MAXGAP=200;
    private static int NEARBY=50000;
    private static int MAXSCORE=300;
    private static int MINPEAKSCORE=7; 
    private static double MINPVALUE=Math.pow(10.0,-MAXSCORE);
    
    //report the region which at least have on fragment score larger than the MINPEAK SCORE
	@Usage

    @Option(doc="This is caseA bam", shortName="A") public File CASEA;
    @Option(doc="This is caseB bam", shortName="B") public File CASEB;
    @Option(doc="This is control A bam", shortName="X") public File CONTROLA;
    @Option(doc="This is control B bam", shortName="Y") public File CONTROLB;
    @Option(doc="caseA peak prefix", shortName="PA") public String CASEAPEAKNAME; 
    @Option(doc="caseB peak prefix", shortName="PB") public String CASEBPEAKNAME; 
    @Option(doc="Chromosome Size",shortName="G") public String CHROMSIZES;
    @Option(doc="output file",shortName="o") public String FOUT;
    @Option(doc="paired end",shortName="p") public boolean ISPAIRED;
    @Option(doc="data format", shortName="f") public String FORMAT;
	
    public static void main(String[] argv){
    	System.exit(new CallPeak5().instanceMain(argv));
        
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
    	 Iterator<? extends Annotation> caseAIter ; 
    	 Iterator<? extends Annotation> caseBIter ; 
    	 Iterator<? extends Annotation> controlAIter ; 
    	 Iterator<? extends Annotation> controlBIter ; 
    	 if(FORMAT.equalsIgnoreCase("shortbed"))
    	 {
    		 ShortBEDTabixDBI caseADBI= new ShortBEDTabixDBI(CASEA.getAbsolutePath(),GENOMESPACE);
        	 ShortBEDTabixDBI caseBDBI= new ShortBEDTabixDBI(CASEB.getAbsolutePath(),GENOMESPACE);
        	 ShortBEDTabixDBI controlADBI= new ShortBEDTabixDBI(CONTROLA.getAbsolutePath(),GENOMESPACE);
        	 ShortBEDTabixDBI controlBDBI= new ShortBEDTabixDBI(CONTROLB.getAbsolutePath(),GENOMESPACE);
        	 caseAIter=caseADBI.iterate();
        	 caseBIter=caseBDBI.iterate();
        	 controlAIter=controlADBI.iterate();
        	 controlBIter=controlBDBI.iterate();
    	 }
    	 else
    	
    	 {
    	 caseAIter = BamIteratorFactory.makeIterator(CASEA,isPaired); 
    	 caseBIter = BamIteratorFactory.makeIterator(CASEB,isPaired); 
    	 controlAIter = BamIteratorFactory.makeIterator(CONTROLA,isPaired); 
    	 controlBIter = BamIteratorFactory.makeIterator(CONTROLB,isPaired); 
    	 }
    	 
    	
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
    	 FileWriter fqvalue = new FileWriter(FOUT.toString()+".qvalue.tab.txt");
    	 SortingCollection<ScoreUSB> qValueArray = SortingCollection.newInstance(ScoreUSB.class, new ScoreUSBCodec(), new ScoreUSBComparator(), 500000); 
    	 
    	 
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
    	 int posPeakIndex=0;
    	 int negPeakIndex=0;
    	 int lastSign=0;
    	 
    	 String posPrefix=CASEAPEAKNAME; 
    	 
    	
    	 
    	 
    	 String negPrefix=CASEBPEAKNAME; 
    	 // -1 is NEG  , +1 is POS
    	 
    	 BedGraphMultiScore bufferBed=new BedGraphMultiScore(0,0,0,new int[1]);
    	 
    	 while(b.hasNext())
    	 {
    		 FocusAndEnv bgraph = b.next();
    		 
    		 k+=1;
    		 if  (k%100000==0)
    			 logger.info(String.format("%d ",k)+bgraph);
    	 
    	     BedGraphMultiScore focus=bgraph.getFocus();
    	     BedGraphMultiScore env=bgraph.getEnvs();
    	    
    	     double envLambdaA = (double)env.getScore()[2]/env.length()*controlAToCaseAIndex;
    	     double envLambdaB = (double)env.getScore()[3]/env.length()*controlBToCaseBIndex;
    	     double fragmentLambdaA=(double)focus.getScore()[2]*controlAToCaseAIndex;
    	     double fragmentLambdaB=(double)focus.getScore()[3]*controlAToCaseAIndex;
    	     
    	     //lambdaA=Math.max(lambdaA, fragmentLambdaA);    	     
			 //lambdaB=Math.max(lambdaB, fragmentLambdaB);    	     
    	     int caseAScore=focus.getScore()[0];
    	     int caseBScore=focus.getScore()[1];
    	     //if ( lambdaA<globalLambdaCaseA ) lambdaA=globalLambdaCaseA;
    	     //if ( lambdaB<globalLambdaCaseB ) lambdaB=globalLambdaCaseB;
    	     if(caseAScore+caseBScore < READS_THRESHOLD) {continue;} //TEST SPEED;
    	     
    	     
    	     double globalScore=getSkellamScore(globalLambdaCaseB,globalLambdaCaseA,caseBScore,caseAScore);
    	     double score=globalScore;
    	     String lambdaString="GLB"; 
    	     double lambdaA=globalLambdaCaseA;
    	     double lambdaB=globalLambdaCaseB;
    	     
    	     if(envLambdaB != 0  && envLambdaA != 0)
    	     {
    	    	 double envScore=getSkellamScore(envLambdaB,envLambdaA,caseBScore,caseAScore);
    	    	 if (Math.abs(score) > Math.abs(envScore))
    	    	 {
    	    		 score=envScore;
    	    		 lambdaString="ENV";
    	    		 lambdaA=envLambdaA;
    	    		 lambdaB=envLambdaB;
    	    	 }
    	     }
    	     if(fragmentLambdaA !=0  && fragmentLambdaB!=0)
    	     {
    	     double fragmentScore=getSkellamScore(fragmentLambdaB,fragmentLambdaA,caseBScore,caseAScore);
    	     	 if (Math.abs(score) > Math.abs(fragmentScore))
    	    	 {
    	    		 score=fragmentScore;
    	    		 lambdaString="FRG";
    	    		 lambdaA=fragmentLambdaA;
    	    		 lambdaB=fragmentLambdaB;
    	    	 }
    	     }
    	     
    	     
    	    if(Math.abs(score)>4)
    	    {
    	    	qValueArray.add(new ScoreUSB(Math.abs(score),focus.length()));
    	    	
    	    	bf.write(jie.getTid2chr().get(focus.getTid()));
    	    	bf.write("\t");
    	    	bf.write(String.format("%d",focus.getStart()));
    	    	bf.write("\t");
    	    	bf.write(String.format("%d",focus.getEnd()));
    	    	bf.write("\t");
    	    	bf.write(String.format("%f",score));
    	    	bf.write("\t");
    	    	bf.write(String.format("[  %s  %d / %.3f  %d / %.3f ]",lambdaString, caseAScore,lambdaA,caseBScore,lambdaB));
    	    	bf.write("\n");
    	    	
    	    	if (focus.getTid()!=lastTid || focus.getStart() - lastEnd > MAXGAP || lastSign != sign(score))  // report buffer.
    	    	{
    	    		
    	    		
    	    		if(Math.abs(maxScore)>MINPEAKSCORE)
    	    		{
    	    			
    	    			if (lastSign>0)
    	    			{
    	    			posPeakIndex+=1;	
    	    			fpeak.write(String.format("%s\t%d\t%d\t%s_peak_%d\t%.2f\t%d\t%d\n", jie.getTid2chr().get(lastTid), lastStart, lastEnd ,posPrefix,posPeakIndex, maxScore,bufferBed.getStart(),bufferBed.getEnd()));
    	    			}
    	    			else
    	    			{
    	    			negPeakIndex+=1;
    	    			fpeak.write(String.format("%s\t%d\t%d\t%s_peak_%d\t%.2f\t%d\t%d\n", jie.getTid2chr().get(lastTid), lastStart, lastEnd , negPrefix,negPeakIndex, maxScore,bufferBed.getStart(),bufferBed.getEnd()));
    	    			}
    	    		}
    	    		lastStart=focus.getStart();
	    		    lastEnd=focus.getEnd();
	    		    lastTid=focus.getTid();
	    		    lastSign=sign(score);
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
    			
    			if (lastSign>0)
    			{
    			posPeakIndex+=1;	
    			fpeak.write(String.format("%s\t%d\t%d\t%s_peak_%d\t%.2f\t%d\t%d\n", jie.getTid2chr().get(lastTid), lastStart, lastEnd ,posPrefix,posPeakIndex, maxScore,bufferBed.getStart(),bufferBed.getEnd()));
    			}
    			else
    			{
    			negPeakIndex+=1;
    			fpeak.write(String.format("%s\t%d\t%d\t%s_peak_%d\t%.2f\t%d\t%d\n", jie.getTid2chr().get(lastTid), lastStart, lastEnd , negPrefix,negPeakIndex, maxScore,bufferBed.getStart(),bufferBed.getEnd()));
    			}
    		}
    	 fpeak.close();
    	bf.close();
    	fout.close();
    	/**
    	 *  Processing Score To Qvalues;
    	 */
    	qValueArray.doneAdding();
    	Iterator<ScoreUSB> iter=qValueArray.iterator();
    	double lastScore=MAXSCORE;
    	int lastLen=0;
    	while(iter.hasNext())
    	{
    		ScoreUSB s=iter.next();
    		double nowScore=(double)Math.floor(s.getScore()*100)/100;
    		if( nowScore!=lastScore && lastLen!=0)
    		{
    		 //report	
    			double expectLen=GENOMESPACE.getLength()*Math.pow(10.0, -lastScore);
    			double qvalue=expectLen/lastLen;
    			fqvalue.write(String.format("%.2f\t%d\t%.2f\t%f\n",lastScore,lastLen,expectLen,qvalue));
    		}
    		lastScore=nowScore;
    		lastLen+=s.getUSB();
    	}
    	//for last one
    	double expectLen=GENOMESPACE.getLength()*Math.pow(10.0, -lastScore);
		double qvalue=expectLen/lastLen;
		fqvalue.write(String.format("%.2f\t%d\t%.2f\t%f\n",lastScore,lastLen,expectLen,qvalue));
    	
		fqvalue.close();
    	
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
	    	if(rightPvalue > MINPVALUE)
	    		return -Math.log10(rightPvalue);
	    	else
	    		return MAXSCORE;
	    }
	    else
	    {
	    	if(leftPvalue > MINPVALUE)
	    		return Math.log10(leftPvalue);
	    	else
	    		return -MAXSCORE;
	    }
	    
	  }
	private int sign(double x)
	{
		if (x>0) return 1;
		if (x<0) return -1;
		return 0;
	}
    
    
}
