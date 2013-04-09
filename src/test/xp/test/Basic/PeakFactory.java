
package xp.test.Basic;

import java.util.ArrayList;
import java.util.Iterator;

import org.apache.log4j.Level;
import org.apache.log4j.Logger;

import cern.colt.Arrays;

import xp.test.Converter.BedGraphMultiScoreReader;
import xp.test.Converter.LocalEnvReader;
import xp.test.Utils.JieCodeSortingCollection;

import nextgen.core.annotation.Annotation;
import nextgen.core.coordinatesystem.GenomicSpace;

/**
 * 2013-4-1
 * 
 * @author zhuxp
 * 
 * read from Iterator<? extends Annotaiton>
 *      into JieCodeSortingCollection  
 *      and then iterate the peaks
 *      
 *      Parameter: 
 *      MAXGAP (default 200);
 *      SCORE THRESHOLD ( in ScoreMachine );
 *      coverageThreshold ( speed up , don't count the score for low coverage region )
 *      TODO: AUTOMATIC define coverage threshold. (according to window size and lambda)
 *      
 * EXAMPLE:
 *       PeakFactory<SkellamScoreMachine>  peakFactory= new PeakFactory<SkellamScoreMachine>(iters,GENOMESPACE,scoreMachine,WINDOWSIZE);
 *       while(peakFactory.hasNext())
 *       {
 *        Peak peak=peakFactory.next();
 *        System.out.println(peak.toShortBed());
 *       }
 *       
 */
public class PeakFactory<F extends ScoreMachine> implements Iterator<Peak> {
	private static final Logger logger = Logger
			.getLogger(PeakFactory.class);
	private Iterator<? extends Annotation>[] annotationIters;
	private String[] annotationIterNames;
	private Long[] annotationSizes;
	
	private int sampleNumber;
	
	private GenomicSpace GENOMESPACE;
	private int MAXGAP=200;
	private JieCodeSortingCollection codes; 
	private double[] globalLambdas;
	private LocalEnvReader  reader;
	private F scoreMachine;
	private int coverageThreshold = 5;
	private int windowSize;
	
	public int getWindowSize() {
		return windowSize;
	}
	public void setWindowSize(int windowSize) {
		this.windowSize = windowSize;
	}
	public GenomicSpace getGENOMESPACE() {
		return GENOMESPACE;
	}
	public void setGENOMESPACE(GenomicSpace gENOMESPACE) {
		GENOMESPACE = gENOMESPACE;
	}
	public int getMAXGAP() {
		return MAXGAP;
	}
	public void setMAXGAP(int mAXGAP) {
		MAXGAP = mAXGAP;
	}
	public double[] getGlobalLambdas() {
		return globalLambdas;
	}
	public void setGlobalLambdas(double[] globalLambdas) {
		this.globalLambdas = globalLambdas;
	}
	public int getCoverageThreshold() {
		return coverageThreshold;
	}
	public void setCoverageThreshold(int coverageThreshold) {
		this.coverageThreshold = coverageThreshold;
	}
	private Peak bufferPeak=null;
	
	private FocusAndEnv bufferFocusEnv=null;
	private Double bufferFocusScore=null;
	private String prefix="Peak";
	
	
	private boolean[] ignoreBlocksDetails; // if true, then count start and end , instead of blocks starts and ends
	
	
	private int counter;
	/**
	 * Initialize
	 * @param iters
	 * @param GS
	 */
	
	public PeakFactory(Iterator<? extends Annotation>[] iters, GenomicSpace GS, F _scoreMachine)
	{  
		this(iters,GS,_scoreMachine,new boolean[iters.length]);
		for (int i = 0; i < ignoreBlocksDetails.length; i++) {
			ignoreBlocksDetails[i]=false;
		}
	}
	public PeakFactory(Iterator<? extends Annotation>[] iters, GenomicSpace GS, F _scoreMachine , boolean[] ignoreBlocksDetails)
	{
		this(iters,GS,_scoreMachine,0, ignoreBlocksDetails);
	}
	public PeakFactory(Iterator<? extends Annotation>[] iters, GenomicSpace GS, F _scoreMachine , int _windowSize)
	{
		this(iters,GS,_scoreMachine,_windowSize, new boolean[iters.length]);
		for (int i = 0; i < ignoreBlocksDetails.length; i++) {
			ignoreBlocksDetails[i]=false;
		}
	}
	public PeakFactory(Iterator<? extends Annotation>[] iters, GenomicSpace GS, F _scoreMachine, int _windowSize, boolean[] ignoreBlocksDetails) 
	{
	 logger.setLevel(Level.DEBUG);
	 annotationIters=iters;
	 sampleNumber=annotationIters.length;
	 GENOMESPACE=GS;	 
	 codes=new JieCodeSortingCollection(GENOMESPACE);
	 windowSize=_windowSize;
	 //codes.setWindowSize(windowSize);
	 //logger.debug("set windowsize " + windowSize);
	 //logger.debug("get windowsize " + codes.getWindowSize());
	 globalLambdas= new double[iters.length];
	 this.ignoreBlocksDetails=ignoreBlocksDetails;
	 initIters();
	 reader = new LocalEnvReader( new BedGraphMultiScoreReader(codes));
	 scoreMachine = _scoreMachine ;
	 bufferPeak=null;
	 counter=0;
	 advance();
	}
	/**
	 *  Reading Annotation Iterators Into Sorting Array.
	 */
	private void initIters()
	{
		int i=0;
		logger.debug(globalLambdas);
		codes.setWindowSize(windowSize); // default is 0, if windowSize is 0, it is coverage we count.
		for(Iterator<? extends Annotation> iter :  annotationIters)
		{
		//To Do: Report the file names?	
		codes.add(iter,i,ignoreBlocksDetails[i]); 
		globalLambdas[i]=(double)codes.getCoverage(i)/GENOMESPACE.getLength();
		i++;
		}
		
	}
	@Override
	public boolean hasNext() {
		// TODO Auto-generated method stub
		if(bufferPeak==null)
		   return false;
		return true;
	}
	@Override
	public Peak next() {
		// TODO Auto-generated method stub
		
		Peak oldBufferPeak = bufferPeak;
		bufferPeak = new Peak();
		advance();
		return oldBufferPeak;
		
	}
	@Override
	public void remove() {
		// TODO Auto-generated method stub
		
	}
	private void advance()
	{
		counter++;
		String name=String.format("%s_%d", prefix, counter);
		Peak peak = new Peak(name);
		
		ArrayList<BedGraphMultiScore> bufferGapBed= new ArrayList<BedGraphMultiScore>();
		ArrayList<Double> bufferGapScores = new ArrayList<Double>();
		
		
		int state;
		if(bufferFocusEnv!=null)
		{
			if(bufferFocusScore > scoreMachine.getThreshold())
			{
				peak.add(bufferFocusEnv.getFocus(),bufferFocusScore);
			    peak.setChr(codes.getTid2chr().get(bufferFocusEnv.getFocus().getTid()));
				bufferFocusEnv=null;
				bufferFocusScore=null;
			    state=1;
			}
			else
			{
				bufferFocusEnv=null;
				bufferFocusScore=null;
				state=0;
			}
			
		}
		else
		{
			state=0;
		}
		
		while(reader.hasNext())
		{
			
		FocusAndEnv focusAndEnv = reader.next();
		BedGraphMultiScore focus = focusAndEnv.getFocus();
		if(focus.getScoreSum() < coverageThreshold && state==0) //not in peak state and coverage is low
			continue;
		
		if(focus.getScoreSum() < coverageThreshold && state==1)
		{
			bufferGapBed.add(focus);
//			bufferGapScores.add(null);
			bufferGapScores.add(0.0);
			continue;
		}
		
		
		BedGraphMultiScore env = focusAndEnv.getEnvs();
		
		
		
		double[] envLambdas = new double[sampleNumber];
		double[] focusLambdas = new double[sampleNumber];
		for (int i = 0; i < sampleNumber; i++) {
			envLambdas[i]=(double)env.getScore()[i]/env.length();
		    focusLambdas[i]=(double)focus.getScore()[i];
		    
		}
		//logger.debug(Arrays.toString(focusLambdas) + " G: " + Arrays.toString(globalLambdas) + " E: "+ Arrays.toString(envLambdas));
		double focusScore = scoreMachine.getScore(focusLambdas,globalLambdas,envLambdas);
		if (focusScore > scoreMachine.getThreshold())  // in peak state
	    {
			state=1;
			if(peak.isEmpty())
			{
				peak.add(focus,focusScore);
				bufferGapBed.clear();
	    		bufferGapScores.clear();
			    peak.setChr(codes.getTid2chr().get(focus.getTid()));
			}
			
			else if(peak.getTid()==focus.getTid()  && (focus.getStart() - peak.getEnd()) < MAXGAP) 
	    	{
	    		for (int i = 0; i < bufferGapBed.size(); i++) {
	    			peak.add(bufferGapBed.get(i),bufferGapScores.get(i));
				}
	    		bufferGapBed.clear();
	    		bufferGapScores.clear();
	    		peak.add(focus,focusScore);
	    	}
	    	else
	    	{
	    		bufferGapBed.clear();
	    		bufferGapScores.clear();
	    		bufferFocusEnv=focusAndEnv;
	    		bufferFocusScore=focusScore;
	    		this.bufferPeak=peak;
	    		return;
	    	}
	    }
		else
		{
			bufferGapBed.add(focus);
			bufferGapScores.add(focusScore);
		}
		
		
		
		}
		if (peak.isEmpty())
			this.bufferPeak=null;
		else
			this.bufferPeak=peak;
		return;
	}
	public String getPrefix() {
		return prefix;
	}
	public void setPrefix(String prefix) {
		this.prefix = prefix;
		String name=bufferPeak.getName();
		//logger.info(name);
		name=name.replace("Peak", prefix);
		//logger.info(name);
		bufferPeak.setName(name);
	}
	
	
    
	
	
	
	

}
