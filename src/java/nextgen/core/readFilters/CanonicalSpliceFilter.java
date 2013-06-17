package nextgen.core.readFilters;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import nextgen.core.alignment.Alignment;
import nextgen.core.annotation.Annotation;

import org.apache.commons.collections15.Predicate;

import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;
import broad.pda.gene.GeneTools;

/**
 * This class will filter reads by whether they have canonical splice junctions
 * Canonical splice junctions are defined as junctions spanning an AG/GT pair
 * 
 * @author mguttman
 *
 */
public class CanonicalSpliceFilter implements Predicate<Alignment>{

	//TODO: ADD A FLAG TO ALLOW FOR NON-CANONICAL SPLICE SITES
	//private Sequence genomeSequence;
	private String genomeSequenceFile=null;
	private String currentChr;
	private Sequence chrSeq;

	public CanonicalSpliceFilter(String genomeFile){
		this.genomeSequenceFile=genomeFile;
	}
	
	public CanonicalSpliceFilter(){
		this.genomeSequenceFile=null;
	}
	
	@Override
	public boolean evaluate(Alignment read) {
		
		//TODO: CHECK
		if(genomeSequenceFile==null){
			//System.out.println("Genome sequence file is not provided");
			if(!read.getSpliceConnections().isEmpty()){return true;}
			//return true;
		}
		if(read.getSpliceConnections().isEmpty()){return false;}
		Sequence chrSeq;
		try {
			chrSeq = getChrSeq(read.getChr());
		} catch (IOException e) {
			return false;
		} //Make sure we have the chromosome sequence
		for(Annotation junction: read.getSpliceConnections()){		
			//This fucntion will allow ONLY canonical sites as well
			String orientation=GeneTools.orientationFromSpliceSites(junction, chrSeq,true);
			if(orientation.equalsIgnoreCase("*")){return false;}
			else{
				return true;
			}

		}
		
		return true;
	}

	private Sequence getChrSeq(String chr) throws IOException {
		if(currentChr==null || chrSeq==null || !chr.equalsIgnoreCase(currentChr)){
			currentChr=chr;
			chrSeq=updateSeq(chr);
		}
		return chrSeq;
	}

	private Sequence updateSeq(String chr) throws IOException {
		FastaSequenceIO fsio = new FastaSequenceIO(genomeSequenceFile);
		if(!(genomeSequenceFile==null)){
			List<String> chrIds = new ArrayList<String>(1);
			chrIds.add(chr);
			List<Sequence> seqs = fsio.extractRecords(chrIds);

			if(!seqs.isEmpty()) {
				return seqs.get(0);
			} 
		}
		return null;
	}

}
