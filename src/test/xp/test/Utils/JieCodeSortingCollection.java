package xp.test.Utils;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;

import org.apache.log4j.Level;
import org.apache.log4j.Logger;


import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.SortingCollection;
import nextgen.core.annotation.Annotation;
import nextgen.core.coordinatesystem.GenomicSpace;

/**
 *  Created on 2013-3-7  
 *  
 *  Log:
 *  2013-03-14 Now adding blocks instead of start and end
 */
public class JieCodeSortingCollection {
	static Logger logger = Logger.getLogger(JieCodeSortingCollection.class.getName());
	
	private HashMap<String,Integer> chr2tid;
	private ArrayList<String> tid2chr;
	private SortingCollection<JieCode> sortingArray;
	private int classIndexState=0;
	private static int MAX_NUM=500000;
	private static File tmpDir = new File(System.getProperty("java.io.tmpdir"));
	private HashMap<Integer, Long> classIndexToReadsNumber = new HashMap<Integer,Long>();
	private HashMap<Integer, Long> classIndexToCoverageLength = new HashMap<Integer,Long>();
    private int classNumber=0;
    private int windowSize=0; // extend reads ? , is half windowSize.
    private int halfWindowSize=0;
    /**
	 * 2013-4-4
	 * extends the ends of reads to halfWindowSize,
	 * now the coverage is how many reads number in the window center on this position.
	 * log: count nearby reads 
	 *         |
	 *       ----->
	 *  <----     --->
	 *   
	 *   
	 */ 
	
	
	public int getClassNumber() {
		return classNumber;
	}
	public int getWindowSize() {
		return windowSize;
	}
	public void setWindowSize(int window) {
		this.windowSize = window;
		this.halfWindowSize = window/2;
		//logger.debug("window"+window);
		//logger.debug("halfWindowSize " + this.halfWindowSize);
		//logger.debug("WindowSize " + this.windowSize);
	}
	
	public void setClassNumber(int classNumber) {
		this.classNumber = classNumber;
	}
	public JieCodeSortingCollection(HashMap<String,Integer> chr2tid) {
		super();
		this.chr2tid=chr2tid;
		this.sortingArray = SortingCollection.newInstance(JieCode.class, new JieCodec(), new JieCodeComparator(), MAX_NUM , tmpDir);
	}
	public JieCodeSortingCollection(GenomicSpace gs)
	{
		setGenomicSpace(gs);
		//logger.setLevel(Level.DEBUG);
		this.sortingArray = SortingCollection.newInstance(JieCode.class, new JieCodec(), new JieCodeComparator(), MAX_NUM , tmpDir);
	}
	public void setGenomicSpace(GenomicSpace gs)
	{
		ArrayList<String> chrs=new ArrayList(gs.getReferenceNames());
		Collections.sort(chrs);
		tid2chr=chrs;
		chr2tid=new HashMap<String, Integer>();
		for(int i=0;i<chrs.size();i++)
		{
			chr2tid.put(chrs.get(i), i);
		}
	}
	
	
	
	public HashMap<String, Integer> getChr2tid() {
		return chr2tid;
	}
	public void setChr2tid(HashMap<String, Integer> chr2tid) {
		this.chr2tid = chr2tid;
	}
	public ArrayList<String> getTid2chr() {
		return tid2chr;
	}
	
	public void add(Iterator<? extends Annotation> iter)
	{
		logger.info("reading Iterator to array ");
		this.classNumber+=1;
		int i=0;
		while (iter.hasNext())
		{
			this.add(iter.next());
			i++;
			if(i%1000000==0) logger.info("reading "+i+" reads");
		}
	}
	public void add(Iterator<? extends Annotation> iter, int classIndex)
	{
		this.classIndexState=classIndex;
		
		if (!classIndexToCoverageLength.containsKey((Integer)classIndex))
		{
			classIndexToCoverageLength.put((Integer)classIndex, Long.valueOf(0));
			
		}
		
		if (!classIndexToReadsNumber.containsKey((Integer)classIndex))
		{
			classIndexToReadsNumber.put((Integer)classIndex, Long.valueOf(0));
			
		}
		add(iter);
		logger.info("Adding Reads Number: " + this.classIndexToReadsNumber.get(classIndex));
		logger.info("Coverage NT number "+ this.classIndexToCoverageLength.get(classIndex));
	}
	
	
	public void add(JieCode a)
	{
		sortingArray.add(a);
		
	};
	public void add(Annotation b, int classIndex)
	/**
	 * if blocks > 1 then  only add blocks region.
	 * else
	 * add blocks + windowSize region ( accumulate the nearby reads to this position ) 
	 * 
	 *  
	 * 
	 */
	{
		List<? extends Annotation> alist=b.getBlocks();
		
	    if(alist.size()==1)
	    {
	    	int startPos=b.getStart();
	    	int stopPos=b.getEnd();
	    	if(halfWindowSize > 0)
	    	{	
	    	startPos-=halfWindowSize;
	    	stopPos+=halfWindowSize;
	    	if(startPos<0) startPos=0;
	    	//TODO IF > chromosome end
	    	logger.debug("windowSize "+windowSize);
	    	
	    	}
	    	 JieCode  start=new JieCode(chr2tid.get(b.getChr()),startPos,true, classIndex);
			 JieCode  stop= new JieCode(chr2tid.get(b.getChr()),stopPos,false, classIndex);
			  
			 sortingArray.add(start);
			 sortingArray.add(stop);
			 
			 int length=stopPos-startPos;
			 //logger.debug(b.getStart() + " TO " +startPos);
			 Long readsnum = classIndexToReadsNumber.get(Integer.valueOf(classIndexState));
			 Long coverage = classIndexToCoverageLength.get(Integer.valueOf(classIndexState));
			 classIndexToReadsNumber.put(Integer.valueOf(classIndexState),readsnum+1);
			 classIndexToCoverageLength.put(Integer.valueOf(classIndex), coverage+Long.valueOf(length));
	    	
	    }
	    else
	    {	
		for(Annotation a: alist)
		{
		 
		 JieCode  start=new JieCode(chr2tid.get(a.getChr()),a.getStart(),true, classIndex);
		 JieCode  stop= new JieCode(chr2tid.get(a.getChr()),a.getEnd(),false, classIndex);
		  
		 sortingArray.add(start);
		 sortingArray.add(stop);
		 int length=stop.getPos()-start.getPos();
		 Long readsnum = classIndexToReadsNumber.get(Integer.valueOf(classIndexState));
		 Long coverage = classIndexToCoverageLength.get(Integer.valueOf(classIndexState));
		 classIndexToReadsNumber.put(Integer.valueOf(classIndexState),readsnum+1);
		 classIndexToCoverageLength.put(Integer.valueOf(classIndex), coverage+Long.valueOf(length));
		 
		}
	    }
	}
	
	
	public Long getCoverage(int classIndex)
	{
		return classIndexToCoverageLength.get(Integer.valueOf(classIndex));
	}
	public Long getReadsNumber(int classIndex)
	{
		return classIndexToReadsNumber.get(Integer.valueOf(classIndex));
	}
	public void add(Annotation a)
	{
		add(a ,this.classIndexState);
	}
	
	
	public Iterator<JieCode> getIterator()	
	{
		return sortingArray.iterator();
	}
	


}
