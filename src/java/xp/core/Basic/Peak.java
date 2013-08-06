
package xp.core.Basic;

import java.util.ArrayList;
import java.util.List;

import org.apache.log4j.Logger;

import cern.colt.Arrays;

/**
 *  Created on 2013-4-1  
 */

public class Peak {
	private static final Logger logger = Logger.getLogger(Peak.class);
	private List<BedGraphMultiScore> listBedGraph;
	private List<Double>  listScores;
	//private List<FocusAndEnv> listFocusAndEnv;
	private long  start;
	private long end;
	private String chr=null;
	private int tid;
	private String name="chr";
	private double score=0.0; //default is max of scores
	
	public Peak()
	{
		listBedGraph=new ArrayList<BedGraphMultiScore>();
		listScores= new ArrayList<Double>();
		
	}
	public Peak(String name)
	{
		this();
		this.setName(name);
		
	}
	
	public String getName() {
		return name;
	}
	public void setName(String name) {
		this.name = name;
	}
	public double getScore() {
		return score;
	}
	public void setScore(double score) {
		this.score = score;
	}
	public boolean isEmpty()
	{
		if (listBedGraph.size()==0)
				return true;
		return false;
	}
	
	public int add(BedGraphMultiScore b)
	{
		if (listBedGraph.size()==0)
		{
			listBedGraph.add(b);
			start=(long)b.getStart();
			end=(long)b.getEnd();
			tid=b.getTid();
		}
		else
		{
			//assert(tid==b.getTid());
			if(tid!=b.getTid())
			{
				logger.warn("Can't add " + b + " to peak" + this +"\n"+ "Skipping adding" + b);
				return 0;
			}
			if(start > b.getStart()) start=(long)b.getStart();
			if(end < b.getEnd()) end=(long)b.getEnd();
			
			listBedGraph.add(b);
			
		}
		return 1;
	}
	
	
	public void add(BedGraphMultiScore b, double score)
	{
		int a=add(b);
		if(a==1) listScores.add(score);
		if(this.score<score) 
			this.score=score; // Score is always the max 
		
	}
	
	
	public long getStart() {
		return start;
	}
	public long getEnd() {
		return end;
	}
	public String getChr() {
		return chr;
	}
	public int getTid()
	{
		return tid;
	}
	public void setTid(int tid)
	{
		this.tid=tid;
	}
	public void setChr(String c)
	{
		chr=c;
	}
	
	
	public int length()
	{
		int len = (int) (getEnd()-getStart());
		return len;
	}
	public int getSampleSize()
	{
		return listBedGraph.get(0).getScore().length;
	}
	public double[] getAverageScores()
	{
		double[] scores = new double[getSampleSize()];
		for(BedGraphMultiScore i : listBedGraph)
		{
			for (int j = 0; j < scores.length; j++) {
				scores[j]+=(double)i.getScore()[j]*i.length();
			}
		}
		for (int i = 0; i < scores.length; i++) {
			scores[i]/=length();
		}
		return scores;
		
	}
	public int[] getSumScores()
	{
		int[] scores = new int[getSampleSize()];
		for(BedGraphMultiScore i : listBedGraph)
		{
			for (int j = 0; j < scores.length; j++) {
				scores[j]+=(double)i.getScore()[j]*i.length();
			}
		}
		return scores;	
	}
	/**
	 * 
	 * @param lambdas
	 * @return normalized coverage.
	 */
	public double[] getNormalizedScores(double[] lambdas)
	{
		double[] scores=getAverageScores();
		assert(scores.length==lambdas.length);
		for (int i = 0; i < scores.length; i++) {
			scores[i]/=lambdas[i];
		}
		return scores;	
	}
	
	public String toSimpleBED()
	{
		
		return String.format("%s\t%d\t%d\t%s\t%.2f", chr,start,end,name,score);
	}
	
	public String toPeakDetail()
	{
		String s="";
		s+=this.toSimpleBED()+"\tAVERAGE:"+Arrays.toString(getAverageScores())+"\n";
		
		for(BedGraphMultiScore i : this.listBedGraph)
		{
			s+="\t"+String.format("%d %d ",i.getStart(),i.getEnd())+Arrays.toString(i.getScore())+"\n";
		}
		return s;
		
	}
	
	
	
	
 }
