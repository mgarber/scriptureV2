package nextgen.core.scripture;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.TreeMap;
import java.util.*;

import org.apache.log4j.Logger;
import org.broad.igv.Globals;

import broad.core.datastructures.IntervalTree;
import broad.core.math.Statistics;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;

/**
 * This class takes an alignment file 
 * @author skadri
 *
 */
public class DESimulatorAssigner {
	
	static Logger logger = Logger.getLogger(DESimulatorAssigner.class.getName());
	
	public DESimulatorAssigner(){
		
	}
	public DESimulatorAssigner(ArgumentMap argMap) throws IOException{
		
		//READ THE INPUT FILE
		BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(new File(argMap.getInput()))));
		
		//output file
		FileWriter writer=new FileWriter(argMap.getOutput());
		
		//Declare all maps
		//Region to line in pro file
		Map<String, String> lineMap = new HashMap<String,String>();
		Map<String, String> geneName = new HashMap<String,String>();
		Map<String, Integer> lengthMap = new HashMap<String,Integer>();
		Map<String, Integer> countsMap = new HashMap<String,Integer>();
		Map<String, Double> abundanceMap = new HashMap<String,Double>();
		IntervalTree<String> tree = new IntervalTree<String>();
		
		String nextLine;
		while ((nextLine = br.readLine()) != null && (nextLine.trim().length() > 0)) {
			String[] entries = nextLine.split(whitespaceDelimiter);
			String region = entries[0];
			String chr = region.split(":")[0];
			Integer start = new Integer(region.split(":")[1].split("-")[0]);
			String dummy = region.split(":")[1].split("-")[1];
			Integer end =  0;
			if(dummy.endsWith("C") || dummy.endsWith("W")){
				//logger.info("Before : "+dummy);
				end = new Integer(dummy.substring(0, dummy.length()-1));
				//logger.info("After = "+end);
			}
			else{
				logger.info("End ends with character other than C or W: "+dummy);
			}			
			String name = entries[1];
			
			//If a gene already overlaps this
			if(tree.overlappingValueIterator(start, end).hasNext() || geneName.containsKey(region)){
				//Write out the gene and do not consider for DE
				//Do not write the line
//				writer.write(nextLine);
				//Add to interval tree
				tree.put(start, end, name);
			}
			else{
				Integer length = new Integer(entries[3]);
				Double abundance = new Double(entries[4]);
				Integer molecules = new Integer(entries[5]);
				lineMap.put(region, nextLine);
				geneName.put(region, name);
				lengthMap.put(region, length);
				countsMap.put(region, molecules);
				abundanceMap.put(region, abundance);
			}
		}
		br.close();
		//NOW OUTPUT TO A LIST THE GENES THAT ARE BEING CONSIDERED.
		FileWriter bw=new FileWriter(argMap.getOutput()+".considered.list.bed");
		for(String region:geneName.keySet()){
			bw.write(region.split(":")[0]+"\t"+region.split(":")[1].split("-")[0]+"\t"+region.split(":")[1].split("-")[1].substring(0, (region.split(":")[1].split("-")[1].length())-1)+"\t"+geneName.get(region)+"\n");
		}
		bw.close();
		
		//Regions to consider
		int num = argMap.getInteger("num",1000);
		if(num<6){
			logger.info("Number of genes cannot be smaller than number of length classes. Setting num of genes to 6");
			num =6;
		}
		logger.info("\nChoosing "+num+" genes");
		
		List<Integer> lengths = new ArrayList<Integer>();
		lengths.addAll(lengthMap.values());
		Collections.sort(lengths);
		
		List<String> chosen = new ArrayList<String>();
		
		int numLenClasses = 6;
		double pct=1.0/(double)numLenClasses;
		//Number of genes to pick in each length class
		int numGenesPerClass = (int)(num/numLenClasses);
		//Take 6 length classes and add random genes from the classes
		for(int i=0;i<numLenClasses;i++){
			logger.info("Pcts: "+pct*(i)+" to "+pct*(i+1));
			double len1 = Statistics.quantile(lengths, pct*(i));
			double len2 = Statistics.quantile(lengths, pct*(i+1));
			logger.info("\nLength class "+(i+1)+" : "+len1+"-"+len2);
			List<String> geneClasses = getNonZeroGenesOfLength(lengthMap,countsMap,len1,len2);
			logger.info("Number of genes in this class: "+geneClasses.size());
			int k=0;
			Random r = new Random();
			if(numGenesPerClass<geneClasses.size()){
				logger.info("Add "+numGenesPerClass+" genes");
				while(k<numGenesPerClass){
					int pMin=0;
					int pMax=geneClasses.size();
					String g = geneClasses.get(pMin + (int)(r.nextFloat() * (pMax - pMin)));
					if(!chosen.contains(g)){
						chosen.add(g);
						k++;
					}
				}
			}
			//Add all
			else{
				logger.info("Add all");
				for(String s:geneClasses)
					chosen.add(s);
			}
		}
		
		Map<String,Integer> counts = new HashMap<String,Integer>();
//        SortedMap<String,Integer> sorted_map = new TreeMap<String,Integer>(bvc);
        counts = getAbundancesOfChosen(countsMap,chosen);
//        sorted_map.putAll(counts);
        Set<Map.Entry<String, Integer>> set = counts.entrySet();
        List<Map.Entry<String, Integer>> list = new ArrayList<Map.Entry<String, Integer>>(set);
        Collections.sort( list, new Comparator<Map.Entry<String, Integer>>()
        {
            public int compare( Map.Entry<String, Integer> o1, Map.Entry<String, Integer> o2 )
            {
                return (o2.getValue()).compareTo( o1.getValue() );
            }
        } );
        
        bw=new FileWriter(argMap.getOutput()+".unchanged.pro");
        for(String s:geneName.keySet()){
        	bw.write(s+"\t"+geneName.get(s)+"\tCDS\t"+lengthMap.get(s)+"\t"+abundanceMap.get(s)+"\t"+countsMap.get(s)+"\n");
        }	
        bw.close();
        
        bw=new FileWriter(argMap.getOutput()+".abundances.changes.txt");
        bw.write("Region\tgeneName\tBefore\tAfter\n");
        boolean odd = true;
        String prev=null;
        //Iterator<Map.Entry<String,Integer>> iter = sorted_map.entrySet().iterator();
        Iterator<Map.Entry<String,Integer>> iter = list.iterator();
        while(iter.hasNext()){
        	if(odd==true){
        		prev=iter.next().getKey();
        		odd=false;
        	}
        	else{
        		int A=counts.get(prev);
        		String s = iter.next().getKey();
        		int B=counts.get(s);
        		//A2 = A/2
        		//A = A/2;
        		//B = B+(A/2)
        		int A4 = A/4;
        		countsMap.put(prev, (int)(A4));
        		countsMap.put(s, (int)(B+(A-A4)));
        		bw.write(prev+"\t"+geneName.get(prev)+"\t"+A+"\t"+A4+"\n");
        		bw.write(s+"\t"+geneName.get(s)+"\t"+B+"\t"+(B+(A-A4))+"\n");
        		odd=true;
        	}
        }
/*        for(String s:sorted_map.keySet()){
        	if(odd==true){
        		prev=s;
        		odd=false;
        	}
        	else{
        		int A=counts.get(prev);
        		int B=counts.get(s);
        		//A2 = A/2
        		//A = A/2;
        		//B = B+(A/2)
        		int A4 = A/4;
        		countsMap.put(prev, (int)(A4));
        		countsMap.put(s, (int)(B+(A-A4)));
        		bw.write(prev+"\t"+geneName.get(prev)+"\t"+A+"\t"+A4+"\n");
        		bw.write(s+"\t"+geneName.get(s)+"\t"+B+"\t"+(B+(A-A4))+"\n");
        		odd=true;
        	}
        }*/
        bw.close();
        
        for(String s:geneName.keySet()){
        	writer.write(s+"\t"+geneName.get(s)+"\tCDS\t"+lengthMap.get(s)+"\t"+abundanceMap.get(s)+"\t"+countsMap.get(s)+"\n");
        }		
		writer.close();
	}
	
	/**
	 * Returns subset of a map
	 */
	private Map<String,Integer> getAbundancesOfChosen(Map<String,Integer> countsMap,Collection<String> chosen){
		Map<String,Integer> rtrn = new HashMap<String,Integer>();
		for(String gene:chosen){
			rtrn.put(gene, countsMap.get(gene));
		}
		return rtrn;
	}
	
	/**
	 * Returns a list of genes in the specified length range
	 * @param geneLengths
	 * @param len1
	 * @param len2
	 * @return
	 */
	private List<String> getNonZeroGenesOfLength(Map<String,Integer> geneLengths,Map<String,Integer> countsMap,double len1,double len2){
		
		List<String> genes = new ArrayList<String>();
		for(String g:geneLengths.keySet()){
			if(geneLengths.get(g)>=len1 && geneLengths.get(g)<=len2 && countsMap.get(g)>0){
				genes.add(g);
			}
		}
		return genes;
	}
	public static String whitespaceDelimiter = "\\s++"; //$NON-NLS-1$

	public static void main(String[] args)throws IOException{
		 
		Globals.setHeadless(true);
		/*
		 * @param for ArgumentMap - size, usage, default task
		 * argMap maps the command line arguments to the respective parameters
		 */
		
		ArgumentMap argMap = CLUtil.getParameters(args,usage,"score");
		
		new DESimulatorAssigner(argMap);
	}

	
	static final String usage = "Usage: DESimulatorAssigner -task score "+
			"\n**************************************************************"+
			"\n\t\tArguments"+
			"\n**************************************************************"+
			
			"\n\t\t-in <.pro file from flux simulator> "+
			"\n\t\t-num <Number of genes to make differentially expressed> "+
			"\n\t\t-out <Output file [Defaults to stdout]> "+
			"\n";

}


