package xp.core.Converter;

import java.util.Iterator;

import org.apache.log4j.Level;
import org.apache.log4j.Logger;

import xp.test.Basic.BedGraphMultiScore;
import xp.test.Utils.JieCode;
import xp.test.Utils.JieCodeSortingCollection;

import broad.core.math.Distribution;


/**
 *  Created on 2013-3-8  
 */
public class BedGraphMultiScoreReader implements Iterator<BedGraphMultiScore>{
	static Logger logger = Logger.getLogger(BedGraphMultiScoreReader.class.getName());
	
	
	
   private JieCodeSortingCollection data;
   private Iterator<JieCode> iter;
   private int last_chr=0;
   private int last_pos=0;
   private int[] coverage_state;
   private int classIndex=0;
   private int READS_THRESHOLD=3;
   private BedGraphMultiScore buffer;
   

public BedGraphMultiScoreReader(JieCodeSortingCollection data) {
	super();
	this.data = data;
	
	this.iter=data.getIterator();
	coverage_state=new int[data.getClassNumber()];
	advance();
	logger.setLevel(Level.DEBUG);
	logger.debug(buffer);
	
	
}

@Override
public boolean hasNext() {
	// TODO Auto-generated method stub
	
	if (buffer==null)
	{
		logger.debug("is null "+buffer);
		return false;
	}
	return true;
}

@Override
public BedGraphMultiScore next() {
	// TODO Auto-generated method stub
	BedGraphMultiScore oldbuffer=buffer;
	advance();
	return oldbuffer;
	
}

@Override
public void remove() {
	// TODO Auto-generated method stub
	
}

private int advance()
{
	
	while(iter.hasNext())
	   {
		 JieCode a=iter.next();
		 if (last_chr!=a.getTid())
		 {
			 //logger.info("new chromsome :"+a.getTid()+" "+tid2chr.get(a.getTid()));
			 
			 for(int j=0;j<coverage_state.length;j++)
			 {
				 coverage_state[j]=0;
			 }
			 last_chr=a.getTid();
			 last_pos=0;
			 
		 }
		 int old_last_pos=last_pos;
		 int[] old_coverage_state=coverage_state.clone();		
		 last_pos=a.getPos();		 
		 if (a.isStart())
			 
			 coverage_state[a.getClassIndex()]+=1;
		 else
			 coverage_state[a.getClassIndex()]-=1;
        
		 if (old_last_pos!=a.getPos())
			 
		 {
			 //logger.debug("in return"+a);
			 buffer=new BedGraphMultiScore(a.getTid(),old_last_pos,a.getPos(), old_coverage_state);
			 return 0;
	
		 }	 
			 
		 
		
	   }
	
	buffer=null;
	return 0;
}
   
   


}
