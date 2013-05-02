package xp.test.Basic;

import java.util.Arrays;

/**
 *  Created on 2013-3-8  
 */
public class BedGraphMultiScore {
    private int tid;
    private int start;
    private int end;
    private int[] score;
	public BedGraphMultiScore(int tid, int old_last_pos, int pos,
			int[] old_coverage_state) {
		// TODO Auto-generated constructor stub
	this.tid=tid;
	this.start=old_last_pos;
	this.end=pos;
	this.score=old_coverage_state;
	
	}
	public BedGraphMultiScore(BedGraphMultiScore b)
	{
	  this(b.getTid(),b.getStart(),b.getEnd(),b.getScore().clone());
	}
	
	public int getTid() {
		return tid;
	}
	public void setTid(int tid) {
		this.tid = tid;
	}
	public int getStart() {
		return start;
	}
	public void setStart(int start) {
		this.start = start;
	}
	public int getEnd() {
		return end;
	}
	public void setEnd(int end) {
		this.end = end;
	}
	public int[] getScore() {
		return score;
	}
	public void setScore(int[] score) {
		this.score = score;
	}
	@Override
	public String toString() {
		return "[ "+ tid + "\t" + start + "\t"
				+ end + "\t" + Arrays.toString(score)+" ]";
	}
	public int length() {
		// TODO Auto-generated method stub
		return getEnd()-getStart();
	}
	public int getScoreLength()
	{
		return score.length;
	}
	public int getScoreSum()
	{
		int s=0;
		for (int i = 0; i < score.length; i++) {
			s+=score[i];
		}
		return s;
	}
	

}
