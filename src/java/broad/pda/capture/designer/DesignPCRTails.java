package broad.pda.capture.designer;

import jaligner.Alignment;
import jaligner.SmithWatermanGotoh;
import jaligner.matrix.Matrix;
import jaligner.matrix.MatrixGenerator;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import broad.core.datastructures.Pair;
import broad.core.motif.SearchException;
import broad.core.primer3.PrimerPair;
import broad.core.primer3.PrimerUtils;
import broad.core.primer3.PcrPrimerDesigner;
import broad.core.sequence.Sequence;
import broad.pda.capture.designer.SmatchLike;

public class DesignPCRTails {

	
	int numDesigns=1;
	int minPerfectMatch=10;
	int minGoodMatch=10;
	Matrix matrix=MatrixGenerator.generate(1.0f, -1.0f);
	int numberOfAnchors=240; //double the number of actuals
	private int minMatch=7;
	Map<PrimerPair, Collection<String>> primers;
	
	public DesignPCRTails() {}
	
	public DesignPCRTails(int numTails) throws Exception{
		primers=getUniquePCRTails(numTails);
	}
	
	public Map<PrimerPair, Collection<String>> getPrimers(){return primers;}
	
	private Map<PrimerPair, Collection<String>> getUniquePCRTails(int numPCRTails) throws Exception {
		//Design PCR Tails
		Collection<PrimerPair> primers=designPCRTails(numPCRTails);
		
		System.err.println("Completed designing artificial tails");
		
		//Add anchored 5bp tails and score primer pairs
		Map<PrimerPair, Collection<String>> anchoredTails=anchorTails(primers);
		
		System.err.println("Completed adding anchors to sequence");
		
		//test uniqueness to the remainder of the primer pool (ie will any given primer pair amplify the wrong thing)
		anchoredTails=filterCrossAmplifyingPrimers(anchoredTails);//TODO Put this back
		
		System.err.println("Completed filtering cross-priming sequences");
				
		return anchoredTails;
	}

	private Map<PrimerPair, Collection<String>> anchorTails(Collection<PrimerPair> uniquePrimers) throws Exception {
		Map<PrimerPair, Collection<String>> rtrn=new TreeMap<PrimerPair, Collection<String>>();
		
		Map<PrimerPair, Collection<String>> anchoredTails=new TreeMap<PrimerPair, Collection<String>>();
		
		Collection<String> randomTails=Sequence.generateAll5Mers();
		
		for(PrimerPair primer: uniquePrimers){
			Collection<String> list=new ArrayList<String>();
			//For each possible tail test cross-priming from major primer
			for(String randomTail: randomTails){
				String leftPrimer=Sequence.get3Prime(primer.getLeftPrimer(),15)+randomTail;
				list.add(leftPrimer);
			}
			anchoredTails.put(primer, list);
		}
		
				
		Collection<String> list=new TreeSet<String>();
		
		//Align each pair to each other and filter possible primer dimer
		for(PrimerPair primer: anchoredTails.keySet()){
			Collection<String> tails=anchoredTails.get(primer);
			for(String tail: tails){
				String left=tail;
				String right=primer.getRightPrimer();
				boolean isComplementary=complementary(left, right);
				if(isComplementary){list.add(tail);}
			}
		}
		
		//Compute crude TM filter by those within range of original
		for(PrimerPair primer: anchoredTails.keySet()){
			Collection<String> tails=anchoredTails.get(primer);
			for(String tail: tails){
				if(!list.contains(tail)){
					double TM=PrimerUtils.computeTM(tail);
					double originalTM=PrimerUtils.computeTM(primer.getLeftPrimer());
					double diff=TM-originalTM;
					double percentGC=PrimerUtils.percentGC(tail);
					double originalGC=PrimerUtils.percentGC(primer.getLeftPrimer());
					if(Math.abs(diff)>1 || Math.abs(percentGC-originalGC)>5){list.add(tail);} //Making sure the GC content is identical
				}
			}
		}
		
		//remove all such primers
		//Won't compute penalties directly
		for(PrimerPair primer: anchoredTails.keySet()){
			Collection<String> tails=anchoredTails.get(primer);
			Collection<String> primerTails=new HashSet<String>();
			for(String tail: tails){
				if(!list.contains(tail)){
					primerTails.add(tail);
				}
			}
			rtrn.put(primer, primerTails);
		}
		
		return rtrn;
	}
	
	/*private Map<PrimerPair, Collection<String>> anchorTails(Collection<PrimerPair> uniquePrimers) throws Exception {
		Map<PrimerPair, Collection<String>> rtrn=new TreeMap<PrimerPair, Collection<String>>();
		
		Map<PrimerPair, Collection<String>> anchoredTails=new TreeMap<PrimerPair, Collection<String>>();
		
		Collection<String> randomTails=Sequence.generateAll5Mers();
		
		for(PrimerPair primer: uniquePrimers){
			Collection<String> list=new ArrayList<String>();
			//For each possible tail test cross-priming from major primer
			for(String randomTail: randomTails){
				String leftPrimer=Sequence.get3Prime(primer.getLeftPrimer(),15)+randomTail;
				list.add(leftPrimer);
				String rightPrimer=primer.getRightPrimer();
				jaligner.Sequence left=new jaligner.Sequence(Sequence.get3Prime(leftPrimer, this.minGoodMatch));
				jaligner.Sequence right=new jaligner.Sequence(Sequence.get3Prime(rightPrimer, this.minGoodMatch));
				for(PrimerPair other: uniquePrimers){
					if(!other.equals(primer)){
						//Remove all tails that are cross-primed
						jaligner.Sequence leftOther=new jaligner.Sequence(Sequence.get3Prime(other.getLeftPrimer(), this.minGoodMatch));
						jaligner.Sequence rightOther=new jaligner.Sequence(Sequence.get3Prime(other.getRightPrimer(), this.minGoodMatch));
						if(crossHyb(left, right, other) || crossHyb(leftOther, rightOther, new Pair<String>(leftPrimer, rightPrimer))){list.remove(leftPrimer);}
					}
					else{
						//just make sure there is no 3' perfect match
						if(crossHyb(left, right, other)){list.remove(leftPrimer);}
					}
				}
			}
			//System.err.println(randomTails.size()+" "+list.size());
			anchoredTails.put(primer, list);
		}
		
		Collection<Pair<String>> allPrimers=new ArrayList<Pair<String>>();
		for(PrimerPair primer: anchoredTails.keySet()){
			Collection<String> leftPrimers=anchoredTails.get(primer);
			for(String leftPrimer: leftPrimers){
				Pair<String> primers=new Pair<String>(leftPrimer, primer.getRightPrimer());
				allPrimers.add(primers);
			}
		}
		
		//Of remaining primers check all cross-priming subprimers
		Collection<String> list=new TreeSet<String>();
		int counter=0;
		for(Pair<String> primer1: allPrimers){
			jaligner.Sequence left=new jaligner.Sequence(Sequence.get3Prime(primer1.getValue1(), this.minGoodMatch));
			jaligner.Sequence right=new jaligner.Sequence(Sequence.get3Prime(primer1.getValue2(), this.minGoodMatch));
			for(Pair<String> primer2: allPrimers){
				if(!primer1.getValue1().equalsIgnoreCase(primer2.getValue1()) && !primer1.getValue2().equalsIgnoreCase(primer2.getValue2())){
					if(crossHyb(left, right, primer2)){list.add(primer1.getValue1());}
				}
			}
			if(counter% 10000 ==0){System.err.println(counter+" "+allPrimers.size());}
			counter++;
		}
		
		int startingSize=list.size();
		System.err.println(allPrimers.size()+" "+list.size());
		
		//Align each pair to each other and filter possible primer dimer
		for(PrimerPair primer: anchoredTails.keySet()){
			Collection<String> tails=anchoredTails.get(primer);
			for(String tail: tails){
				if(!list.contains(tail)){
					String left=tail;
					String right=primer.getRightPrimer();
					boolean isComplementary=complementary(left, right);
					if(isComplementary){list.add(tail);}
				}
			}
		}
		
		System.err.println(startingSize+" Further reduced to "+list.size());
		
		startingSize=list.size();
		//Compute crude TM filter by those within range of original
		for(PrimerPair primer: anchoredTails.keySet()){
			Collection<String> tails=anchoredTails.get(primer);
			for(String tail: tails){
				if(!list.contains(tail)){
					double TM=PrimerUtils.computeTM(tail);
					double originalTM=PrimerUtils.computeTM(primer.getLeftPrimer());
					double diff=TM-originalTM;
					double percentGC=PrimerUtils.percentGC(tail);
					double originalGC=PrimerUtils.percentGC(primer.getLeftPrimer());
					if(Math.abs(diff)>1 || Math.abs(percentGC-originalGC)>5){list.add(tail);}
				}
			}
		}
		
		System.err.println(startingSize+" Further reduced to "+list.size());
		
		//remove all such primers
		//compute primer penalty for remaining primers
		//Won't compute penalties
		counter=0;
		for(PrimerPair primer: anchoredTails.keySet()){
			long start=System.currentTimeMillis();
			Collection<String> tails=anchoredTails.get(primer);
			Collection<String> primerTails=new TreeSet<String>();
			for(String tail: tails){
				if(!list.contains(tail)){
					//String seq=tail+insertSeq+Sequence.reverseSequence(primer.getRightPrimer());
					//Collection<PrimerPair> primers=qPCRPrimerDesigner.designSyntheticPrimers(seq, tail, primer.getRightPrimer(), 1);
					//if(primers!=null && !primers.isEmpty()){
						//PrimerPair primer2=primers.iterator().next();
						//primerTails.add(primer2);
						//System.err.println("1 good primer");
					//}
					primerTails.add(tail);
				}
			}
			rtrn.put(primer, primerTails);
			counter++;
			long end=System.currentTimeMillis();
			System.err.println("Scoring primer pairs "+counter+" Took "+(end-start));
		}
		
		return rtrn;
	}*/

	public boolean complementary(String left, String right) {
		String rightRev=Sequence.reverseSequence(right);
		
		//align left and rightRev
		jaligner.Sequence leftSeq=new jaligner.Sequence(left);
		jaligner.Sequence rightSeq=new jaligner.Sequence(rightRev);
		
		Alignment align=SmithWatermanGotoh.align(leftSeq, rightSeq, matrix, 1000000, 1000000);
		
		if(align.getNumberOfMatches()>this.minMatch){return true;}
		
		return false;
	}

	private Map<PrimerPair, Collection<String>> filterCrossAmplifyingPrimers(Map<PrimerPair, Collection<String>> anchoredPrimers) throws SearchException {
		//Filter major primers that inappropriately amplify other major primer or non-self sub-primers
		Collection<PrimerPair> crossPrimed=new TreeSet<PrimerPair>();
		Collection<String> tailCrossPrimed=new TreeSet<String>();
		
		
		
		//Get all subprimers
		Collection<Pair<String>> tailPrimers=new HashSet<Pair<String>>();
		for(PrimerPair primer: anchoredPrimers.keySet()){
			if(!crossPrimed.contains(primer)){
				Collection<String> tails=anchoredPrimers.get(primer);
				for(String tail: tails){
					if(!tailCrossPrimed.contains(tail)){
						Pair<String> tailPrimer=new Pair<String>(tail, primer.getRightPrimer());
						tailPrimers.add(tailPrimer);
					}
				}
			}
		}
		
		
		for(PrimerPair primer: anchoredPrimers.keySet()){
			for(PrimerPair other: anchoredPrimers.keySet()){
				if(!primer.equals(other)){
					boolean crossPrimes=crossPrime(primer, other);
					if(crossPrimes){crossPrimed.add(other);}
				}
			}
			if(!crossPrimed.contains(primer)){
				String kmer=Sequence.get3Prime(primer.getLeftPrimer(), this.minGoodMatch);
				SmatchLike smatch=new SmatchLike(kmer, tailPrimers, primer);
				Collection<String> matches=smatch.getForwardTargets(0);
				tailCrossPrimed.addAll(matches);
			}
			
		}
		
		System.err.println("Filtered major primers that cross prime...");
		
		int counter=0;
		for(Pair<String> subprimer1: tailPrimers){
			if(!tailCrossPrimed.contains(subprimer1.getValue1())){
				String kmer=Sequence.get3Prime(subprimer1.getValue1(), this.minGoodMatch);
				SmatchLike smatch=new SmatchLike(kmer, tailPrimers, subprimer1);
				Collection<String> matches=smatch.getForwardTargets(0);
				tailCrossPrimed.addAll(matches);
				//System.err.println("Reduced from "+tailPrimers.size()+" to "+matches.size());
			}
			for(PrimerPair primer: anchoredPrimers.keySet()){
				if(!crossPrimed.contains(primer)){
					Pair<String> major=new Pair<String>(primer.getLeftPrimer(), primer.getRightPrimer());
					boolean crossPrimes=crossPrime(subprimer1, major);
					if(crossPrimes){tailCrossPrimed.add(subprimer1.getValue1());}
				}	
			}
			counter++;
			if(counter%500 ==0){System.err.println(counter);}
		}
		
		System.err.println("Reduced from "+tailPrimers.size()+" to "+tailCrossPrimed.size());
		
			
		System.err.println("Filtered subprimers that cross prime...");
		
		Map<PrimerPair, Collection<String>> rtrn=new TreeMap<PrimerPair, Collection<String>>();
		
		for(PrimerPair primer: anchoredPrimers.keySet()){
			if(!crossPrimed.contains(primer)){
				Collection<String> tails=anchoredPrimers.get(primer);
				tails.removeAll(tailCrossPrimed);
				//tails.removeAll(uncertainTails);
				rtrn.put(primer, tails);
			}
		}
		
		return rtrn;
	}

	
	private boolean crossPrime(Pair<String> primer1, Pair<String> primer2) throws SearchException {
		String kmer=Sequence.get3Prime(primer1.getValue1(), this.minPerfectMatch);
		SmatchLike smatch=new SmatchLike(kmer, primer2.getValue1(), kmer.length());
		Collection<String> matches=smatch.getForwardTargets(0);
		if(matches==null || matches.isEmpty()){return false;}
		
		//Test if primer1 primes primer2
		jaligner.Sequence left=new jaligner.Sequence(Sequence.get3Prime(primer1.getValue1(), this.minGoodMatch));
		jaligner.Sequence right=new jaligner.Sequence(Sequence.get3Prime(primer1.getValue2(), this.minGoodMatch));
		
		jaligner.Sequence leftRandom=new jaligner.Sequence(primer2.getValue1());
		jaligner.Sequence rightRandom=new jaligner.Sequence(primer2.getValue2());
		Alignment leftAlignment=SmithWatermanGotoh.align(left, leftRandom, matrix, 1000000, 1000000);
		Alignment rightAlignment=SmithWatermanGotoh.align(right, rightRandom, matrix, 1000000, 1000000);
		
		return crossHyb(leftAlignment, left.getSequence().toCharArray().length)&& crossHyb(rightAlignment, right.getSequence().toCharArray().length);
	}
	
	private boolean crossPrime(PrimerPair primer1, Pair<String> primer2) throws SearchException {
		String kmer=Sequence.get3Prime(primer1.getLeftPrimer(), this.minPerfectMatch);
		SmatchLike smatch=new SmatchLike(kmer, primer2.getValue1(), kmer.length());
		Collection<String> matches=smatch.getForwardTargets(0);
		if(matches==null || matches.isEmpty()){return false;}
		
		//Test if primer1 primes primer2
		jaligner.Sequence left=new jaligner.Sequence(Sequence.get3Prime(primer1.getLeftPrimer(), this.minGoodMatch));
		jaligner.Sequence right=new jaligner.Sequence(Sequence.get3Prime(primer1.getRightPrimer(), this.minGoodMatch));
		
		jaligner.Sequence leftRandom=new jaligner.Sequence(primer2.getValue1());
		jaligner.Sequence rightRandom=new jaligner.Sequence(primer2.getValue2());
		Alignment leftAlignment=SmithWatermanGotoh.align(left, leftRandom, matrix, 1000000, 1000000);
		Alignment rightAlignment=SmithWatermanGotoh.align(right, rightRandom, matrix, 1000000, 1000000);
		
		return crossHyb(leftAlignment, left.getSequence().toCharArray().length)&& crossHyb(rightAlignment, right.getSequence().toCharArray().length);
	}
	
	public boolean crossPrime(PrimerPair primer1, PrimerPair primer2) throws SearchException {
		String kmer=Sequence.get3Prime(primer1.getLeftPrimer(), this.minPerfectMatch);
		SmatchLike smatch=new SmatchLike(kmer, primer2.getLeftPrimer(), kmer.length());
		Collection<String> matches=smatch.getForwardTargets(0);
		if(matches==null || matches.isEmpty()){return false;}
		
		//Test if primer1 primes primer2
		jaligner.Sequence left=new jaligner.Sequence(Sequence.get3Prime(primer1.getLeftPrimer(), this.minGoodMatch));
		jaligner.Sequence right=new jaligner.Sequence(Sequence.get3Prime(primer1.getRightPrimer(), this.minGoodMatch));
		
		jaligner.Sequence leftRandom=new jaligner.Sequence(primer2.getLeftPrimer());
		jaligner.Sequence rightRandom=new jaligner.Sequence(primer2.getRightPrimer());
		Alignment leftAlignment=SmithWatermanGotoh.align(left, leftRandom, matrix, 1000000, 1000000);
		Alignment rightAlignment=SmithWatermanGotoh.align(right, rightRandom, matrix, 1000000, 1000000);
		
		return crossHyb(leftAlignment, left.getSequence().toCharArray().length)&& crossHyb(rightAlignment, right.getSequence().toCharArray().length);
	}

	//Cross-hybs are defined by large number of matches in the 3' end of the primer
	private boolean crossHyb(Alignment leftAlignment, int length) {
		//Look at the total number of mismatches
		boolean mismatchCutoff=mismatchCutoff(leftAlignment, length);
		
		if(!mismatchCutoff){return false;}
		
		jaligner.Sequence seq1=new jaligner.Sequence(Sequence.get3Prime(leftAlignment.getOriginalSequence1().getSequence(), this.minPerfectMatch));
		jaligner.Sequence seq2=leftAlignment.getOriginalSequence2();
		Alignment perfect=SmithWatermanGotoh.align(seq1, seq2, matrix, 1000000, 1000000);
		boolean perfect3Prime=perfectMatch(perfect);
		if(!perfect3Prime){return false;}
		
		//for(int i=0; i<perfect.getLength(); i++){System.err.println(perfect.getSequence1()[i]+" "+perfect.getSequence2()[i]);}
		
		return true;
	}

	private boolean perfectMatch(Alignment perfect) {
		return (perfect.getNumberOfMatches()==this.minPerfectMatch);
	}

	private boolean mismatchCutoff(Alignment leftAlignment, int length) {
		int minMatch=(length/2);
		return (leftAlignment.getLength()>minMatch && leftAlignment.getNumberOfMatches()>minMatch);
	}

	private Collection<PrimerPair> designPCRTails(int numPCRTails) throws IOException {
		Collection<PrimerPair> rtrn=new TreeSet<PrimerPair>();
		
		for(int i=0; i<numPCRTails; i++){
			//System.err.println("Designing tail "+i);
			//generate random sequences
			String seq=Sequence.generateRandomSequence(500);
			
			//have primer3 pick good primers (regardless of length distance)
			Collection<PrimerPair> primers=PcrPrimerDesigner.designSyntheticPrimers(seq, numDesigns);
			if(primers!=null && !primers.isEmpty()){
				PrimerPair primer=primers.iterator().next();
				rtrn.add(primer);
			}
		}
	
		return rtrn;
	}
	
	public static void main(String[] args) throws Exception{
		String leftPrimer=args[0];
		String rightPrimer=args[1];
		
		String seq=(leftPrimer)+Sequence.generateRandomSequence(500)+Sequence.reverseSequence(rightPrimer);
		System.err.println(seq);
		Collection<PrimerPair> primers=PcrPrimerDesigner.designSyntheticPrimers(seq, leftPrimer, rightPrimer, 1);
		if(primers!=null && !primers.isEmpty()){
			PrimerPair primer=primers.iterator().next();
			System.err.println(primer.getPrimerPairPenalty()+" "+primer.getLeftPrimerTM()+" "+primer.getRightPrimerTM());
		}
		
		System.err.println(Sequence.reverseSequence("AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"));
	}
	
	static String usage=" args[0]=num designs \n args[1]=save";
	
}
