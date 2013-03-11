package xp.test.Converter;

import java.util.ArrayList;
import java.util.Iterator;


import org.apache.log4j.Level;
import org.apache.log4j.Logger;

import xp.test.Basic.BedGraphMultiScore;
import xp.test.Basic.FocusAndEnv;
import xp.test.Utils.ShortBEDCodec;

/**
 *  Created on 2013-3-9  
 */
public class LocalEnvReader implements Iterator<FocusAndEnv>{
	private static Logger logger = Logger.getLogger(LocalEnvReader.class.getName());
	
	
    private ArrayList<BedGraphMultiScore>  env = new ArrayList<BedGraphMultiScore>();  //local environment
    private BedGraphMultiScoreReader reader;
    private int focusIndex=0; // index of current point 
    
    private int[]  envSum ;// sum of value in Env
    private int envSize = 10000;
    private BedGraphMultiScore buffer;
    private BedGraphMultiScore bed; //corrent bed (focus)
    private FocusAndEnv bufferNext;
    
	public LocalEnvReader(BedGraphMultiScoreReader reader) {
		super();
		this.reader = reader;
		buffer=reader.next();
		bed=buffer;
		envSum=new int[buffer.getScore().length];
		bufferNext=advance();
	
		
		//logger.setLevel(Level.DEBUG);	//Set to DEBUG
		
		
	}

	
	
	
	
	
	
	public boolean hasNext() {
		// TODO Auto-generated method stub
		if(bufferNext==null)
			return false;
		return true;
	}

	@Override
	public FocusAndEnv next() {
		// TODO Auto-generated method stub
		
		FocusAndEnv oldBufferNext=bufferNext;
		bufferNext=advance();
		return oldBufferNext; //need to trim;
	}
	
	

	@Override
	public void remove() {
		// TODO Auto-generated method stub
		
	}
	
	private void updateEnvSum(BedGraphMultiScore b, boolean add)
	{
		for(int i=0;i<envSum.length;i++)
		{
			if (add)
				envSum[i]+=b.getScore()[i]*b.length();
			else
				envSum[i]-=b.getScore()[i]*b.length();
		}
	}
	private FocusAndEnv advance()
	{
	  if( env.size() > focusIndex+1)   // still reading next in env
	  {
		  focusIndex+=1;
		  bed = env.get(focusIndex);
		  
	  }
	  else
	  {
		  if( buffer!=null)        // focus is the last of env, need to move the next in buffer
		  {
			  bed=buffer;
		  }
		  else
		  {
			  
			  bed=null;
			  env.clear();
			  for(int i=0;i<envSum.length;i++)
			      envSum[i]=0;
			  logger.warn("return NULL!!!");
			  return null;        //ending  no buffer and have read all the env;
		  }
	  }
		
	  /*
	   *  First Part ADDING TO ENV
	   */
	   
	  logger.debug("advance To" + bed);
	  
	  if (buffer!=null && buffer.getTid() == bed.getTid())
	  {
		  if(buffer.getStart() < bed.getEnd() + envSize/2 )
		  {
			  env.add(buffer); updateEnvSum(buffer,true); // add buffer into env and update env sum;
			 // while(reader.hasNext())
			  while(true)
			  {
				  buffer=reader.next();
				  if (buffer==null) break;
				  if(buffer.getTid()==bed.getTid() && buffer.getStart() < bed.getEnd() + envSize/2 )
				  {
					  env.add(buffer); updateEnvSum(buffer,true);
					  logger.debug("adding to env" + buffer);
					  logger.debug("the env sum 0:" + envSum[0]);
					  
				  }
				  else
				  {
					  break;
				  }
				  
			  }
		  }
	  }
	  /*
	   * Second Part ROMOVING FROM ENV
	   */
	 
	  while(true)
	  {
		  if (env.size()==0)
				  break;
		  BedGraphMultiScore b=env.get(0);
		  if (b.getTid()!=bed.getTid() || b.getEnd() < bed.getStart() - envSize/2)
		  {
			  env.remove(b);updateEnvSum(b,false);
			  logger.debug("remove from env" + b);
			  logger.debug("the env sum 0:" + envSum[0]);
			  
			  focusIndex-=1;   //updata focusIndex;
		  }
		  else
		  {
			  break;
		  }
		  
	  }
	  /*
	   * update bufferNext;
	   */
	  int[] trimEnvs=envSum.clone();
	  int envStart=env.get(0).getStart();
	  int envEnd=env.get(env.size()-1).getEnd();
	  int trimStart=bed.getStart()-envSize/2-envStart; 
	  int trimEnd=envEnd-(bed.getEnd()+envSize/2);
	  for(int i=0;i<envSum.length;i++)
	  {
		 trimEnvs[i]-=trimStart*env.get(0).getScore()[i];
		 trimEnvs[i]-=trimEnd*env.get(env.size()-1).getScore()[i];
	  }
	  BedGraphMultiScore envs = new BedGraphMultiScore(bed.getTid(),bed.getStart()-envSize/2,bed.getEnd()+envSize/2,trimEnvs);
	  return new FocusAndEnv(bed,envs);
		
	}
	
}
