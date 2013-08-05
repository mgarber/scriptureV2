
package xp.test.Utils;

import java.util.ArrayList;


import org.apache.log4j.Logger;



/**
 *  Created on 2013-4-20  
 */

public class RNAUtils {
	@SuppressWarnings("unused")
	private static final Logger logger = Logger.getLogger(RNAUtils.class);
	/**
	 * parse secondary structure string (such as "(((...)))")
	 * and
	 * get stems 
	 * stem format Integer[4] {a,b,c,d} a-b pairs with c-d
	 * 1-index
	 * example:
	 *  input "((..))((.))"
	 *  output [1,2,5,6]  [7,8,10,11]
	 *    
	 * @param SecondaryStructureString
	 * @return stems ArrayList
	 */
	public static ArrayList<Integer []> parseStems(String ss)
	
	{
		ArrayList<Integer []> retv = new ArrayList<Integer []>();
		SecondaryStructure s=new SecondaryStructure(ss);
		retv=s.getLinks();
		return retv;
	}
	/**
	 * parse secondary structure string  (such as "(((...)))")
	 * and
	 * get loops
	 * 1-index
	 * 
	 * example:
	 *  input "((..))((.))"
	 *  output [3,4]  [9,9]
	 * @param ss
	 * @return loops ArrayList
	 */
	public static ArrayList<Integer []> parseLoops(String ss)
	{
		ArrayList<Integer []> retv = new ArrayList<Integer []>();
		int i=0;
		int state=0;//if in loop state=1;
		int lastStart=0;
		for(i=0;i<ss.length();i++)
		{
			if(ss.charAt(i)=='.')
			{
				if(state==0)
				{
					lastStart=i;
					state=1;
				}
			}
			else
			{
				if(state==1)
				{
					Integer[] loop={lastStart+1,i};
					retv.add(loop);
					state=0;
				}
				
			}
			
		}
		if (state==1) 
		{
			Integer[] loop={lastStart+1,i};
			retv.add(loop);
			state=0;
		}
		return retv;
	}
    
}







class SecondaryStructure {
	private static final Logger logger = Logger
			.getLogger(SecondaryStructure.class);
	private ArrayList<Integer> posRegister;
	private ArrayList<Integer[]> pairRegister;
	private ArrayList<Integer[]> pairs;
	private ArrayList<Integer[]> links;
	private int pos;
	private int state;
	private String structure;
	
	
	
	public ArrayList<Integer[]> getLinks() {
		return links;
	}

	public void setLinks(ArrayList<Integer[]> links) {
		this.links = links;
	}

	public String getStructure() {
		return structure;
	}

	public void setStructure(String structure) {
		this.structure = structure;
	}

	public SecondaryStructure(String structure) {
		super();
		this.structure = structure;
		posRegister=new ArrayList<Integer>();
		pairRegister=new ArrayList<Integer []>();
		pairs=new ArrayList<Integer[]>();
		links=new ArrayList<Integer[]>();
		pos=0;
		read(structure);
	}

	private void read(String S)
	{
		
		for (int i = 0; i < S.length(); i++) {
			readNext(S.charAt(i));
			
		}
		dumpPairRegister();
	}
	
	private void readNext(char s)
	{
		
		if(s == '(') readNext(1);
		if(s == ')') readNext(-1);
		if(s == '.') readNext(0);
		
		
	}
	private void readNext(int i)
	{
		pos+=1;
		if(i==1) 
		{
			posRegister.add(pos);
			state=1;
		}
		else if (i==-1)
		{   assert(posRegister.size() > 0);
			if(posRegister.size()==0)
			{
				logger.error("wrong secondary structure, ignore unpaired ) in" + pos);
				return;
			}
			int first_pos=posRegister.get(posRegister.size()-1);
			
			
			posRegister.remove(posRegister.size()-1);
			
			Integer[] newPair={first_pos,pos};
			
			if (pairRegister.size()>0)
				{
				Integer[] lastPair=pairRegister.get(pairRegister.size()-1);
				if(lastPair[0]-newPair[0]>1 || newPair[1]-lastPair[1] > 1)
				{
					dumpPairRegister();
				}
				}
			
			pairRegister.add(newPair.clone());
			pairs.add(newPair.clone());
			state=-1;
		}
		else if(i==0)
		{
			if (state!=0)
			state=0;
		}
		
	}
	private void dumpPairRegister()
	{
		if(pairRegister.size()==0)
			return;
		Integer[] end = pairRegister.get(0);
		Integer[]  start = pairRegister.get(pairRegister.size()-1);
		Integer[] l={start[0],end[0],end[1],start[1]};
		links.add(l);
		pairRegister.clear();
	}
	
}