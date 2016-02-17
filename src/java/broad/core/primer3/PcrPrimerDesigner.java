package broad.core.primer3;


import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

import org.apache.log4j.Logger;


import nextgen.core.annotation.Gene;

import broad.core.annotation.GenomicAnnotation;
import broad.core.parser.CommandLineParser;
import broad.core.parser.StringParser;
import broad.core.primer3.Primer3SequenceInputTags.SequenceRegionCoordinates;
import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;


public class PcrPrimerDesigner  {
	public static int n=1;
	public static final int MIN_PRIMER_DES_SPACE = 150; 
	public static final int MIN_PROD_SIZE = 100;
	public static final String USAGE = "Usage: PrimerDesigner TASK=<task_num> <task_args>\n" +
	"\tTasks:\n" +
	"\t\t1. prime gaps: OUTDIR=output directory IN=<sequence AGP> file MAXPRODSIZE=maximum product size BUFFER=size of buffer between primers and target\n" +
	"\t\t2. Generate short amplicon set. This is useful to say desing qPCR exeperiments \n"+ 
		"\t\t\tIN=<sequence AGP>\n"+
		"\t\t\tOUTDIR=<output directory>\n"+
		"\t\t\tDIST=<Distance between amplicons> OR NUM<number of amplicons>\n"+
		"\t\t\tAMPSIZE=<Desired amplicon size>\n"+
		"\t\t\tTARGETSTART=<If supplied, the start of region to cover>\n"+
		"\t\t\tTARGETEND=<If supplied, the end of region to cover>\n" +
		"\t\t\tREPEATOUT=<RepeatMasker output file for the sequence>\n" +
	"\t\t3. Generate gene set amplification set -in <gene list, a gene symbol per line> -minsize <minimum amplicon size> -cdsonly <include this flag if only interested in cds only amplification>\n" +
		"\n\t\t-maxsize <max amplicon size> [-numOfDesigns <number of different designs to ouput (default 1)> -buffer <intronic bases to include (default 0)>] " +
	"\n\t\t4. Generate flanking primers -in <Multifasta file default is standard in> -outdir <directory where to write output files> -optimalDistFromEnds <Optimal distance from sequence begining and end" +
	"\n\t\t\t -maxDistFromEnds <Maximum distance from ends, used when no primer could be found within optimal distance> [-outprefix <a prefix to prepend to output, if not specified no prefix will be prepended to file names>]";
	
	private ArrayList<SequenceRegionCoordinates> exludedAnnotations = new ArrayList<SequenceRegionCoordinates>();
	private static Logger logger = Logger.getLogger(PcrPrimerDesigner.class.getName());

	public PcrPrimerDesigner() {}
	
	/*public PrimedRegion findPrimers(Sequence sequence, 
			GenomicAnnotation targetAnnotation, 
			int numberOfPrimers, 
			int maxSize) 
		throws IllegalArgumentException, IllegalAccessException, IOException {
		
		return findPrimers(sequence, targetAnnotation, numberOfPrimers, maxSize, null, null);
		
	}*/
	
	/*public PrimedRegion findPrimers(Sequence sequence, 
			GenomicAnnotation targetAnnotation, 
			int numberOfPrimers, 
			int maxSize,
			String leftPrimer,
			String rightPrimer) 
		throws IllegalArgumentException, IllegalAccessException, IOException {
		PrimedRegion pr = new PrimedRegion(targetAnnotation);
		
		Primer3SequenceInputTags p3sit = new Primer3SequenceInputTags(sequence);
		p3sit.setPrimerSequenceId(targetAnnotation.getName());
		SequenceRegionCoordinates target = new SequenceRegionCoordinates((int) targetAnnotation.getStart() - buffer,(int) targetAnnotation.getEnd() + buffer);
		logger.info("Designing primers for Sequence " + sequence.getId() + " \nBASES:" + sequence.getSequenceBases() +
				" target " + target);
		p3sit.addTarget(target);		
		p3sit.addExcludedRegions(exludedAnnotations);
		if(leftPrimer != null) {
			p3sit.setPrimerLeftInput(leftPrimer);
		} 
		
		if(rightPrimer != null) {
			p3sit.setPrimerRightInput(rightPrimer);
		}
		
		for(int i = 1; i <= numberOfPrimers; i++) {
			int confId	 = 0;
			PrimerPair lastPair = pr.getLastPrimerPair();
			if(lastPair != null && lastPair.hasPrimers() ) {
				if(leftPrimer == null || leftPrimer.length() == 0) {
					int leftPrimerMiddle = (int) Math.ceil((2*lastPair.getLeftPrimerPosition() + lastPair.getLeftPrimer().length())/2);
					p3sit.addExcludedRegion(new SequenceRegionCoordinates(leftPrimerMiddle - 1, leftPrimerMiddle) );
				}
				
				if(rightPrimer == null || rightPrimer.length() == 0) {
					int rightPrimerMiddle = (int) Math.ceil((2*lastPair.getRightPrimerPosition() - lastPair.getRightPrimer().length())/2);
					p3sit.addExcludedRegion(new SequenceRegionCoordinates(rightPrimerMiddle - 1, rightPrimerMiddle) );
				}
			}

			PrimerPair pp = new PrimerPair();
			while(!pp.hasPrimers() && confId < configurations.length) {
				Primer3Configuration conf = configurations[confId++].copy();
				conf.minProductSize = (int) Math.min(conf.minProductSize, target.length());
				conf.maxProductSize = conf.minProductSize;
				while(!pp.hasPrimers()&& conf.maxProductSize < maxSize) {
					conf.maxProductSize = conf.maxProductSize + (int)(conf.minProductSize/10);
					pp = p3io.findPrimerPair(p3sit, conf);
					logger.info("Primer has more primers? "+pp.hasPrimers() + " # confs " + configurations.length + " doing configuration " + confId + " minProductSize " + conf.minProductSize + " maxProdSize " + conf.maxProductSize);
					StringWriter sw = new StringWriter();
					BufferedWriter bw = new BufferedWriter(sw);
					p3io.writeRecord(bw,p3sit, conf);
					bw.flush();
				}
			}
			pp.setPrimerPairId(targetAnnotation.getName() + "_" + i);
			pr.addPrimerPair(pp);
			StringWriter sw = new StringWriter();
			BufferedWriter bw = new BufferedWriter(sw);
			if(pp.hasPrimers()) {				
				pp.setProductSequence(sequence.getSequenceBases().substring(pp.getLeftPrimerPosition() , pp.getRightPrimerPosition()  + 1));
				p3io.writeRecord(bw, p3sit, pp.getConfigurationUsed());
				logger.info(sw);
			}

		}
		return pr;
	}*/
	

	@SuppressWarnings("unused")
	private void addExcludedAnnotations(List<? extends GenomicAnnotation> excludeList) {
		Iterator<? extends GenomicAnnotation> it = excludeList.iterator();
		while(it.hasNext()) {
			GenomicAnnotation annot = it.next();
			this.exludedAnnotations.add(new SequenceRegionCoordinates(annot.getStart(), annot.getEnd()));
		}
	}

	
	
	
	/*public static PrimerPair designPrimers(String sequenceDirectory, Alignments align)throws Exception{
		
		String sequenceFile=sequenceDirectory+"/"+align.getChr().replaceAll("chr", "").trim()+"/"+align.getChr()+".agp";
		Chromosome chr = new Chromosome(sequenceFile);
		chr.loadSequence();
		
		SequenceRegion target = new SequenceRegion("chr"+chr.getSymbol());
		target.setRegionStart(align.getStart());
		target.setRegionEnd(align.getEnd());
		
		Primer3Configuration best = Primer3ConfigurationFactory.getqRTPCRConfiguration();
		logger.info("Max probe size " + best.maxProductSize + " min " + best.minProductSize+" "+target.toString());
		
		Primer3Configuration [] configurations = { best};				
		Primer3IO p3io = new Primer3IO();
		p3io.startPrimer3Communications();
			
		String suffix = "_" + align.toFileName();
		//BufferedWriter primerOutput = new BufferedWriter(new FileWriter(outDir + "/primers"+suffix+".txt"));
		//BufferedWriter primerConfig   = new BufferedWriter(new FileWriter(outDir + "/primersInput"+suffix+".primer3"));
		//String primerFastaFile   =outDir + "/primers"+suffix+".fa";
		//BufferedWriter primerGFF = new BufferedWriter(new FileWriter(outDir + "/primers"+suffix+".gff"));
		
				
		qPCRPrimerDesigner pd = new qPCRPrimerDesigner(configurations, p3io, 0);
		pd.addExcludedAnnotations(chr.getRepeatReader().getRepeatList());
		ArrayList<Sequence> primerSequences = new ArrayList<Sequence>();
			
		
		SequenceRegion probeRegion = new SequenceRegion("");
		chr.getRegion(target);
		probeRegion.setName("probe_" + target.toString());
		Primer3SequenceInputTags p3sit = new Primer3SequenceInputTags(target);
		p3sit.setPrimerSequenceId(target.getId());
		PrimerPair pp = p3io.findPrimerPair(p3sit, best);
		
			
		p3io.endPrimer3Communications();
		//logger.info(pp.getPrimerSequences().get(0));
			
		
		
		return pp;
	}*/
	
	

	/*public static PrimerPair designPrimers(Chromosome chr, Alignments align, boolean repeatMask)throws Exception{
		
		SequenceRegion target = new SequenceRegion("chr"+chr.getSymbol());
		target.setRegionStart(align.getStart());
		target.setRegionEnd(align.getEnd());
		
		Primer3Configuration best = Primer3ConfigurationFactory.getqRTPCRConfiguration();
		logger.info("Max probe size " + best.maxProductSize + " min " + best.minProductSize+" "+target.toString());
		
		Primer3Configuration [] configurations = { best};				
		Primer3IO p3io = new Primer3IO();
		p3io.startPrimer3Communications();
			
		String suffix = "_" + align.toFileName();
		//BufferedWriter primerOutput = new BufferedWriter(new FileWriter(outDir + "/primers"+suffix+".txt"));
		//BufferedWriter primerConfig   = new BufferedWriter(new FileWriter(outDir + "/primersInput"+suffix+".primer3"));
		//String primerFastaFile   =outDir + "/primers"+suffix+".fa";
		//BufferedWriter primerGFF = new BufferedWriter(new FileWriter(outDir + "/primers"+suffix+".gff"));
		
				
		qPCRPrimerDesigner pd = new qPCRPrimerDesigner(configurations, p3io, 0);
		pd.addExcludedAnnotations(chr.getRepeatReader().getRepeatList());
		ArrayList<Sequence> primerSequences = new ArrayList<Sequence>();
			
		
		SequenceRegion probeRegion = new SequenceRegion("");
		chr.getRegion(target, repeatMask);
		probeRegion.setName("probe_" + target.toString());
		Primer3SequenceInputTags p3sit = new Primer3SequenceInputTags(target);
		p3sit.setPrimerSequenceId(target.getId());
		
		
		//add targets p3sit.
		
		try{PrimerPair pp = p3io.findPrimerPair(p3sit, best);
		
			
		p3io.endPrimer3Communications();
		//logger.info(pp.getPrimerSequences().get(0));
			
		
		
		return pp;
		}catch(NullPointerException ex){return null;}
	}*/
	
	
	/*public static PrimerPair designPrimers(Chromosome chr, RefSeqGene gene, boolean repeatMask)throws Exception{
		
		
		Primer3Configuration best = Primer3ConfigurationFactory.getqRTPCRConfiguration();
		
		Primer3Configuration [] configurations = {best};				
		Primer3IO p3io = new Primer3IO();
		p3io.startPrimer3Communications();
			
		
				
		qPCRPrimerDesigner pd = new qPCRPrimerDesigner(configurations, p3io, 0);
		ArrayList<Sequence> primerSequences = new ArrayList<Sequence>();
			
		Sequence mRNA=new Sequence(gene.getName());
		mRNA.setSequenceBases(gene.getSequence(chr, repeatMask, false));
		
		Primer3SequenceInputTags p3sit = new Primer3SequenceInputTags(mRNA);	
		p3sit.setPrimerSequenceId(gene.getName());
		
		ArrayList<Integer> splicePositions=gene.getSpliceJunctionCoordinates();
		//logger.info(targetList);
		ArrayList<SequenceRegionCoordinates> targetList=new ArrayList();
		for(Integer spliceJunction: splicePositions){
			SequenceRegionCoordinates region=new SequenceRegionCoordinates(spliceJunction-n, spliceJunction+n);
			targetList.add(region);
		}
		p3sit.addTargets(targetList);
		//add targets p3sit.
		
		
		try{
			PrimerPair pp = p3io.findPrimerPair(p3sit, best);
			p3io.endPrimer3Communications();
			return pp;
		}catch(NullPointerException ex){return null;}
	}*/
	
	
	/*public static Set<PrimerPair> designPrimers(Chromosome chr, RefSeqGene gene, boolean repeatMask, int numDesigns)throws Exception{
		
		
		Primer3Configuration best = Primer3ConfigurationFactory.getqRTPCRConfiguration();
		
		Primer3Configuration [] configurations = {best};				
		Primer3IO p3io = new Primer3IO();
		p3io.startPrimer3Communications();
			
		
				
		qPCRPrimerDesigner pd = new qPCRPrimerDesigner(configurations, p3io, 0);
		ArrayList<Sequence> primerSequences = new ArrayList<Sequence>();
			
		Sequence mRNA=new Sequence(gene.getName());
		mRNA.setSequenceBases(gene.getSequence(chr, repeatMask, false));
		logger.info(mRNA.getSequenceBases());
		
		Primer3SequenceInputTags p3sit = new Primer3SequenceInputTags(mRNA);	
		p3sit.setPrimerSequenceId(gene.getName());
		
		ArrayList<Integer> splicePositions=gene.getSpliceJunctionCoordinates();
		//logger.info(targetList);
		ArrayList<SequenceRegionCoordinates> targetList=new ArrayList();
		for(Integer spliceJunction: splicePositions){
			SequenceRegionCoordinates region=new SequenceRegionCoordinates(spliceJunction-n, spliceJunction+n);
			targetList.add(region);
		}
		p3sit.addTargets(targetList);
		//add targets p3sit.
		
		//try{
			Set<PrimerPair> primers=new HashSet();
			for(int i=0; i<numDesigns; i++){
				PrimerPair pp = p3io.findPrimerPair(p3sit, best);
				if(pp!=null &&pp.getLeftPrimer()!=null && pp.getRightPrimer()!=null){
					p3sit.addExcludedRegion(new SequenceRegionCoordinates(pp.getLeftPrimerPosition(), pp.getLeftPrimerPosition()+pp.getLeftPrimer().toCharArray().length));
					p3sit.addExcludedRegion(new SequenceRegionCoordinates(pp.getRightPrimerPosition()-pp.getRightPrimer().toCharArray().length, pp.getRightPrimerPosition()));
					logger.info(gene.getAlignment());
					primers.add(pp);
				}
				//return pp;
			}
			p3io.endPrimer3Communications();
			return primers;
		//}catch(NullPointerException ex){return null;}
	}*/
	
	//TODO should add all junctions as sequence regions
	public static Set<PrimerPair> designJenPrimers(Sequence chr, Gene gene, boolean repeatMask, int numDesigns, String pathPrimer3core)throws Exception{
		
		Primer3Configuration best = Primer3ConfigurationFactory.getJenRTPCRConfiguration();
		
		Primer3IO p3io = new Primer3IO(pathPrimer3core);
		p3io.startPrimer3Communications();
				
		Sequence mRNA=new Sequence(gene.getName());
		mRNA.setSequenceBases(gene.getSequence(chr, repeatMask, false));
		//logger.info(mRNA.getSequenceBases());
		
		Primer3SequenceInputTags p3sit = new Primer3SequenceInputTags(mRNA);	
		p3sit.setPrimerSequenceId(gene.getName());
		
		/**Add exon junction sequence targets**/
		ArrayList<Integer> splicePositions=gene.getSpliceJunctionCoordinates();
		logger.info(gene.getName()+" "+gene.getNumExons());
		//logger.info(targetList);
		ArrayList<SequenceRegionCoordinates> targetList=new ArrayList<SequenceRegionCoordinates>();
		for(Integer spliceJunction: splicePositions){
			SequenceRegionCoordinates region=new SequenceRegionCoordinates(spliceJunction.intValue()-1, spliceJunction.intValue()+1);
			targetList.add(region);
		}
		p3sit.addTargets(targetList);
		//add targets p3sit.
		/***/
		
		//try{
			Set<PrimerPair> primers=new HashSet<PrimerPair>();
			
				Collection<PrimerPair> pp = p3io.findPrimerPair(p3sit, best);
				
					primers.addAll(pp);
				
				//return pp;
			
			p3io.endPrimer3Communications();
			return primers;
		//}catch(NullPointerException ex){return null;}
	}
	

	
	/*public static Set<PrimerPair> designqPCRAcrossAndWithJunctions(Chromosome chr, RefSeqGene gene, boolean repeatMask, int numDesigns)throws Exception{
		
		
		Primer3Configuration best = Primer3ConfigurationFactory.getJenRTPCRConfiguration();
		
		Primer3Configuration [] configurations = {best};				
		Primer3IO p3io = new Primer3IO();
		p3io.startPrimer3Communications();
			
		
				
		qPCRPrimerDesigner pd = new qPCRPrimerDesigner(configurations, p3io, 0);
		ArrayList<Sequence> primerSequences = new ArrayList<Sequence>();
			
		Sequence mRNA=new Sequence(gene.getName());
		mRNA.setSequenceBases(gene.getSequence(chr, repeatMask, false));
		logger.info(mRNA.getSequenceBases());
		
		Primer3SequenceInputTags p3sit = new Primer3SequenceInputTags(mRNA);	
		p3sit.setPrimerSequenceId(gene.getName());
		
		ArrayList<Integer> splicePositions=gene.getSpliceJunctionCoordinates();
		//logger.info(targetList);
		ArrayList<SequenceRegionCoordinates> targetList=new ArrayList();
		for(Integer spliceJunction: splicePositions){
			SequenceRegionCoordinates region=new SequenceRegionCoordinates(spliceJunction-n, spliceJunction+n);
			targetList.add(region);
		}
		p3sit.addTargets(targetList);
		//add targets p3sit.
		
		
		//try{
			Set<PrimerPair> primers=new HashSet();
			for(int i=0; i<numDesigns; i++){
				PrimerPair pp = p3io.findPrimerPair(p3sit, best);
				if(pp!=null &&pp.getLeftPrimer()!=null && pp.getRightPrimer()!=null){
					p3sit.addExcludedRegion(new SequenceRegionCoordinates(pp.getLeftPrimerPosition(), pp.getLeftPrimerPosition()+pp.getLeftPrimer().toCharArray().length));
					p3sit.addExcludedRegion(new SequenceRegionCoordinates(pp.getRightPrimerPosition()-pp.getRightPrimer().toCharArray().length, pp.getRightPrimerPosition()));
					logger.info(gene.getAlignment());
					primers.add(pp);
				}
				//return pp;
			}
			p3io.endPrimer3Communications();
			return primers;
		//}catch(NullPointerException ex){return null;}
	}*/
	
	/*public static Set<PrimerPair> designqPCRWithinExons(Chromosome chr, RefSeqGene gene, boolean repeatMask, int numDesigns)throws Exception{
		
		
		Primer3Configuration best = Primer3ConfigurationFactory.getJenRTPCRConfiguration();
		
		Primer3Configuration [] configurations = {best};				
		Primer3IO p3io = new Primer3IO();
		p3io.startPrimer3Communications();
			
		
				
		qPCRPrimerDesigner pd = new qPCRPrimerDesigner(configurations, p3io, 0);
		ArrayList<Sequence> primerSequences = new ArrayList<Sequence>();
			
		Sequence mRNA=new Sequence(gene.getName());
		mRNA.setSequenceBases(gene.getSequence(chr, repeatMask, false));
		logger.info(mRNA.getSequenceBases());
		
		Primer3SequenceInputTags p3sit = new Primer3SequenceInputTags(mRNA);	
		p3sit.setPrimerSequenceId(gene.getName());
		
		//ArrayList<SequenceRegionCoordinates> targetList=gene.getSpliceJunctionCoordinates(chr, repeatMask, 0);
		//logger.info(targetList);
		//p3sit.addTargets(targetList);
		//add targets p3sit.
		
		//try{
			Set<PrimerPair> primers=new HashSet();
			for(int i=0; i<numDesigns; i++){
				PrimerPair pp = p3io.findPrimerPair(p3sit, best);
				if(pp!=null &&pp.getLeftPrimer()!=null && pp.getRightPrimer()!=null){
					p3sit.addExcludedRegion(new SequenceRegionCoordinates(pp.getLeftPrimerPosition(), pp.getLeftPrimerPosition()+pp.getLeftPrimer().toCharArray().length));
					p3sit.addExcludedRegion(new SequenceRegionCoordinates(pp.getRightPrimerPosition()-pp.getRightPrimer().toCharArray().length, pp.getRightPrimerPosition()));
					logger.info(gene.getAlignment());
					primers.add(pp);
				}
				//return pp;
			}
			p3io.endPrimer3Communications();
			return primers;
		//}catch(NullPointerException ex){return null;}
	}*/
	
	/*public static Set<PrimerPair> designPrimers(Chromosome chr, Alignments gene, boolean repeatMask, int numDesigns)throws Exception{
		
		
		Primer3Configuration best = Primer3ConfigurationFactory.getqRTPCRConfiguration();
		
		Primer3Configuration [] configurations = {best};				
		Primer3IO p3io = new Primer3IO();
		p3io.startPrimer3Communications();
			
		
				
		qPCRPrimerDesigner pd = new qPCRPrimerDesigner(configurations, p3io, 0);
		ArrayList<Sequence> primerSequences = new ArrayList<Sequence>();
			
		Sequence mRNA=new Sequence(gene.getName());
		mRNA.setSequenceBases(gene.getSequence(chr, repeatMask));
		logger.info(mRNA.getSequenceBases());
		
		Primer3SequenceInputTags p3sit = new Primer3SequenceInputTags(mRNA);	
		p3sit.setPrimerSequenceId(gene.getName());
		
		
		
		//try{
			Set<PrimerPair> primers=new HashSet();
			for(int i=0; i<numDesigns; i++){
				PrimerPair pp = p3io.findPrimerPair(p3sit, best);
				if(pp!=null &&pp.getLeftPrimer()!=null && pp.getRightPrimer()!=null){
					p3sit.addExcludedRegion(new SequenceRegionCoordinates(pp.getLeftPrimerPosition(), pp.getLeftPrimerPosition()+pp.getLeftPrimer().toCharArray().length));
					p3sit.addExcludedRegion(new SequenceRegionCoordinates(pp.getRightPrimerPosition()-pp.getRightPrimer().toCharArray().length, pp.getRightPrimerPosition()));
					logger.info(gene);
					primers.add(pp);
				}
				//return pp;
			}
			p3io.endPrimer3Communications();
			return primers;
		//}catch(NullPointerException ex){return null;}
	}*/
	
	/*public static Set<PrimerPair> designChIPPrimers(Chromosome chr, Alignments gene, boolean repeatMask, int numDesigns)throws Exception{
		
		
		Primer3Configuration best = Primer3ConfigurationFactory.getChIPqPCRConfiguration();
		
		Primer3Configuration [] configurations = {best};				
		Primer3IO p3io = new Primer3IO();
		p3io.startPrimer3Communications();
			
		
				
		qPCRPrimerDesigner pd = new qPCRPrimerDesigner(configurations, p3io, 0);
		ArrayList<Sequence> primerSequences = new ArrayList<Sequence>();
			
		Sequence mRNA=new Sequence(gene.getName());
		mRNA.setSequenceBases(gene.getSequence(chr, repeatMask));
		logger.info(mRNA.getSequenceBases());
		
		Primer3SequenceInputTags p3sit = new Primer3SequenceInputTags(mRNA);	
		p3sit.setPrimerSequenceId(gene.getName());
		
		
		
		//try{
			Set<PrimerPair> primers=new HashSet();
			for(int i=0; i<numDesigns; i++){
				PrimerPair pp = p3io.findPrimerPair(p3sit, best);
				if(pp!=null &&pp.getLeftPrimer()!=null && pp.getRightPrimer()!=null){
					p3sit.addExcludedRegion(new SequenceRegionCoordinates(pp.getLeftPrimerPosition(), pp.getLeftPrimerPosition()+pp.getLeftPrimer().toCharArray().length));
					p3sit.addExcludedRegion(new SequenceRegionCoordinates(pp.getRightPrimerPosition()-pp.getRightPrimer().toCharArray().length, pp.getRightPrimerPosition()));
					logger.info(gene);
					primers.add(pp);
				}
				//return pp;
			}
			p3io.endPrimer3Communications();
			return primers;
		//}catch(NullPointerException ex){return null;}
	}*/
	
	
	/*public static Set<PrimerPair> designPrimers(Sequence mRNA, int numDesigns)throws Exception{
		
		
		Primer3Configuration best = Primer3ConfigurationFactory.getqRTPCRConfiguration();
		
		Primer3Configuration [] configurations = {best};				
		Primer3IO p3io = new Primer3IO();
		p3io.startPrimer3Communications();
			
		
				
		qPCRPrimerDesigner pd = new qPCRPrimerDesigner(configurations, p3io, 0);
		ArrayList<Sequence> primerSequences = new ArrayList<Sequence>();
			
		Primer3SequenceInputTags p3sit = new Primer3SequenceInputTags(mRNA);	
		p3sit.setPrimerSequenceId(mRNA.getId());
		
		
		
		//try{
			Set<PrimerPair> primers=new HashSet();
			for(int i=0; i<numDesigns; i++){
				PrimerPair pp = p3io.findPrimerPair(p3sit, best);
				if(pp!=null &&pp.getLeftPrimer()!=null && pp.getRightPrimer()!=null){
					p3sit.addExcludedRegion(new SequenceRegionCoordinates(pp.getLeftPrimerPosition(), pp.getLeftPrimerPosition()+pp.getLeftPrimer().toCharArray().length));
					p3sit.addExcludedRegion(new SequenceRegionCoordinates(pp.getRightPrimerPosition()-pp.getRightPrimer().toCharArray().length, pp.getRightPrimerPosition()));
					logger.info(mRNA.getId());
					primers.add(pp);
				}
				//return pp;
			}
			p3io.endPrimer3Communications();
			return primers;
		//}catch(NullPointerException ex){return null;}
	}
	*/
	
	

	public static class PrimedRegion {
		ArrayList<PrimerPair> primers = new ArrayList<PrimerPair>();
		GenomicAnnotation targetAnnotation;
		
		public PrimedRegion(GenomicAnnotation annotation) {
			this.targetAnnotation = annotation;
		}
		
		public Collection<? extends Sequence> getPrimerSequences() {
			ArrayList<Sequence> primerSequences = new ArrayList<Sequence>(); 
			for(int i = 0; i < primers.size(); i++) {
				primerSequences.addAll(primers.get(i).getPrimerSequences());
			}
			return primerSequences;
		}

		public void addPrimerPair(PrimerPair pp) {
			primers.add(pp);
		}
		
		public PrimerPair getLastPrimerPair() {
			return primers.size() > 0 ? primers.get(primers.size() - 1) : null;
		}
		
		public List<PrimerPair> getPrimers() { return primers; }
		
		public void write(BufferedWriter goodPrimerWriter,  BufferedWriter goodPrimerConfig, Primer3IO p3io, boolean printExpectedProduct) throws IOException {
			Iterator<PrimerPair> it = primers.iterator();
			PrimerPair pp = null;
			while(it.hasNext()) {
				pp = it.next();
				logger.info(pp);
				if(pp.hasPrimers()) {
					goodPrimerWriter.write(pp.toString(printExpectedProduct));
					goodPrimerWriter.newLine();
					pp.writePrimer3Record(goodPrimerConfig, p3io);
				} 
			}
		}
		
		public void write(BufferedWriter goodPrimerWriter,  BufferedWriter goodPrimerConfig, Primer3IO p3io) throws IOException {
			write(goodPrimerWriter, goodPrimerConfig, p3io, false);
		}
		
	}


	public static Collection<PrimerPair> designIntronPrimers(Sequence chr, Gene gene, boolean repeatMask, String sequencePrimer, String sequencePrimerRevComp, int numDesigns, int min3, int min5, String pathPrimer3core) throws Exception {
		Primer3Configuration best = Primer3ConfigurationFactory.getJenRTPCRConfiguration();
		
		if(numDesigns>0){best.setPrimerNumReturn(numDesigns);}
		if(min3>0){best.setPrimerMin3PrimeOverlapOfJunction(min3);}
		if(min5>0){best.setPrimerMin5PrimeOverlapOfJunction(min5);}
		
		Primer3IO p3io = new Primer3IO(pathPrimer3core);
		p3io.startPrimer3Communications();
				
		Sequence mRNA=new Sequence(gene.getName());
		mRNA.setSequenceBases(gene.getSequence(chr, repeatMask, false));
		ArrayList<Integer> splicePositions=gene.getSpliceJunctionCoordinates();
		Primer3SequenceInputTags p3sit = new Primer3SequenceInputTags(mRNA);	
		p3sit.setPrimerSequenceId(gene.getName());
		p3sit.addJunctions(splicePositions);
		if(sequencePrimer!=null && !sequencePrimer.isEmpty()){p3sit.setPrimerLeftInput(sequencePrimer);}
		if(sequencePrimerRevComp!=null && !sequencePrimerRevComp.isEmpty()){p3sit.setPrimerRightInput(sequencePrimerRevComp);}
		
		Collection<PrimerPair> pp = p3io.findPrimerPair(p3sit, best);
		

		p3io.endPrimer3Communications();
		
		TreeSet<PrimerPair> rtrn=new TreeSet<PrimerPair>();
		for(PrimerPair pair: pp){
			if(pair.getLeftPrimer()!=null && pair.getRightPrimer()!=null){rtrn.add(pair);}
		}
		
		return rtrn;
	}
	
	
	public static Collection<PrimerPair> designRACEPrimers(Sequence mRNA, String leftPrimer, String rightPrimer, String pathPrimer3core) throws Exception {
		Primer3Configuration best = Primer3ConfigurationFactory.getJenRTPCRConfiguration();
				
		Primer3IO p3io = new Primer3IO(pathPrimer3core);
		p3io.startPrimer3Communications();
			
		//ArrayList<Integer> splicePositions=gene.getSpliceJunctionCoordinates();
		Primer3SequenceInputTags p3sit = new Primer3SequenceInputTags(mRNA);	
		p3sit.setPrimerSequenceId(mRNA.getId());
		//p3sit.addJunctions(splicePositions);
		if(leftPrimer!=null && !leftPrimer.isEmpty()){p3sit.setPrimerLeftInput(leftPrimer);}
		if(rightPrimer!=null && !rightPrimer.isEmpty()){p3sit.setPrimerRightInput(rightPrimer);}
		
		Collection<PrimerPair> pp = p3io.findPrimerPair(p3sit, best);
		

		p3io.endPrimer3Communications();
		
		
		
		return pp;
	}

	public static Collection<PrimerPair> designPCRPrimers(Sequence chr, Gene gene, boolean repeatMask,	int numDesigns, boolean crossJunction, String pathPrimer3core) throws Exception {
		Primer3Configuration best = Primer3ConfigurationFactory.getJenRTPCRConfiguration();
		
		if(numDesigns>0){best.setPrimerNumReturn(numDesigns);}
		if(!crossJunction){repeatMask=true;}		
		
		Primer3IO p3io = new Primer3IO(pathPrimer3core);
		p3io.startPrimer3Communications();
			
		
				
		Sequence mRNA=new Sequence(gene.getName());
		mRNA.setSequenceBases(gene.getSequence(chr, repeatMask, false));
		ArrayList<Integer> splicePositions=gene.getSpliceJunctionCoordinates();
		Primer3SequenceInputTags p3sit = new Primer3SequenceInputTags(mRNA);	
		p3sit.setPrimerSequenceId(gene.getName());
		
		if(crossJunction){
			ArrayList<SequenceRegionCoordinates> targetList=new ArrayList<SequenceRegionCoordinates>();
			for(Integer spliceJunction: splicePositions){
				SequenceRegionCoordinates region=new SequenceRegionCoordinates(spliceJunction.intValue()-1, spliceJunction.intValue()+1);
				targetList.add(region);
			}
			p3sit.addTargets(targetList);
		}
		
		
		Collection<PrimerPair> pp = p3io.findPrimerPair(p3sit, best);
		

		p3io.endPrimer3Communications();
		
		TreeSet<PrimerPair> rtrn=new TreeSet<PrimerPair>();
		for(PrimerPair pair: pp){
			if(pair.getLeftPrimer()!=null && pair.getRightPrimer()!=null){rtrn.add(pair);}
		}
		
		return rtrn;
	}
	
	public static Collection<PrimerPair> designPCRPrimers(Sequence chr, Gene gene, boolean repeatMask,	int numDesigns,  SequenceRegionCoordinates target, String pathPrimer3core) throws Exception {
		Primer3Configuration best = Primer3ConfigurationFactory.getJenRTPCRConfiguration();
		
		if(numDesigns>0){best.setPrimerNumReturn(numDesigns);}	
		
		Primer3IO p3io = new Primer3IO(pathPrimer3core);
		p3io.startPrimer3Communications();		
				
		Sequence mRNA=new Sequence(gene.getName());
		mRNA.setSequenceBases(gene.getSequence(chr, repeatMask, false));
		Primer3SequenceInputTags p3sit = new Primer3SequenceInputTags(mRNA);	
		p3sit.setPrimerSequenceId(gene.getName());
		p3sit.addTarget(target);
		
		
		Collection<PrimerPair> pp = p3io.findPrimerPair(p3sit, best);
		

		p3io.endPrimer3Communications();
		
		TreeSet<PrimerPair> rtrn=new TreeSet<PrimerPair>();
		for(PrimerPair pair: pp){
			if(pair.getLeftPrimer()!=null && pair.getRightPrimer()!=null){rtrn.add(pair);}
		}
		
		return rtrn;
	}
	
	public static Collection<PrimerPair> designCloningPrimers(Sequence geneSequence, int numDesigns, SequenceRegionCoordinates target, String pathPrimer3core) throws Exception {
		Primer3Configuration best = Primer3ConfigurationFactory.getLongRangePCRConfiguration();
		
		//if(numDesigns>0){best.setPrimerNumReturn(numDesigns);}
		
		
		Primer3IO p3io = new Primer3IO(pathPrimer3core);
		p3io.startPrimer3Communications();
			
				
		Primer3SequenceInputTags p3sit = new Primer3SequenceInputTags(geneSequence);	
		p3sit.setPrimerSequenceId(geneSequence.getId());
		if(target!=null){p3sit.addTarget(target);}
		//Collection<PrimerPair> pp = p3io.findPrimerPair(p3sit, best);
		
		Collection<PrimerPair> primers=new TreeSet<PrimerPair>();
		for(int i=0; i<numDesigns; i++){
			Collection<PrimerPair> pp = p3io.findPrimerPair(p3sit, best);
			if(pp.iterator().hasNext()){
			PrimerPair primer=pp.iterator().next();
			if(primer!=null && primer.getLeftPrimer()!=null && primer.getRightPrimer()!=null){
				p3sit.addExcludedRegion(new SequenceRegionCoordinates(primer.getLeftPrimerPosition(), primer.getLeftPrimerPosition()+primer.getLeftPrimer().toCharArray().length));
				p3sit.addExcludedRegion(new SequenceRegionCoordinates(primer.getRightPrimerPosition()-primer.getRightPrimer().toCharArray().length, primer.getRightPrimerPosition()));
				//logger.info(mRNA.getId());
				primers.add(primer);
			}
			}
			//return pp;
		}
		
		p3io.endPrimer3Communications();
		
		Collection<PrimerPair> rtrn=new TreeSet<PrimerPair>();
		for(PrimerPair pair: primers){
			if(pair.getLeftPrimer()!=null && pair.getRightPrimer()!=null){rtrn.add(pair);}
		}
		
		return rtrn;
	}
	
	public static Collection<PrimerPair> designCloningPrimers(Sequence geneSequence, int numDesigns, String pathPrimer3core) throws Exception {
		return designCloningPrimers(geneSequence, numDesigns, null, pathPrimer3core);
	}

	public static Collection<PrimerPair> designSyntheticPrimers(String seq, int numDesigns, String pathPrimer3core, double optimalMeltingTemp) throws IOException {
		Primer3Configuration best = Primer3ConfigurationFactory.getSyntheticConfiguration(optimalMeltingTemp);
		
		if(numDesigns>0){best.setPrimerNumReturn(numDesigns);}
				
		Primer3IO p3io = new Primer3IO(pathPrimer3core);
		p3io.startPrimer3Communications();
			
		
				
		Sequence mRNA=new Sequence("Gene");
		mRNA.setSequenceBases(seq);
		
		Primer3SequenceInputTags p3sit = new Primer3SequenceInputTags(mRNA);	
		p3sit.setPrimerSequenceId(mRNA.getId());
		
				
		Collection<PrimerPair> pp = p3io.findPrimerPair(p3sit, best);
		

		p3io.endPrimer3Communications();
		
		TreeSet<PrimerPair> rtrn=new TreeSet<PrimerPair>();
		for(PrimerPair pair: pp){
			if(pair.getLeftPrimer()!=null && pair.getRightPrimer()!=null){rtrn.add(pair);}
		}
		
		return rtrn;
	}
	
	/**
	 * Design a primer pair that flanks a region of interest
	 * One primer is in each flanking region
	 * Neither primer overlaps the actual region
	 * @param config primer3 configuration
	 * @param sequence The sequence containing the region of interest
	 * @param pathPrimer3core primer3core executable
	 * @param regionStart Start position of region of interest
	 * @param regionEnd Position after last position of region of interest
	 * @param flankingRegionSize Size of flanking regions to search for primers
	 * @return Primer pair with one primer in each flank or null if none exist
	 * @throws IOException
	 */
	public static PrimerPair designPrimerPairFlankingWindow(Primer3Configuration config, Sequence sequence, String pathPrimer3core, int regionStart, int regionEnd, int flankingRegionSize) throws IOException {

		Sequence leftFlank = sequence.getSubSequence("", regionStart - flankingRegionSize, regionStart);
		Sequence rightFlank = sequence.getSubSequence("", regionEnd, regionEnd + flankingRegionSize);
		
		// Make a string of Ns to put between the two flanking regions
		int innerLength = Math.max(regionEnd - regionStart, flankingRegionSize);
		char[] ns = new char[innerLength];
		for(int i = 0; i < ns.length; i++) ns[i] = 'N';
		String nStr = new String(ns);
		
		String modifiedSequence = leftFlank.getSequenceBases() + nStr + rightFlank.getSequenceBases(); // Sequence to design primers against
		
		// Change the product size in the primer3 config
		config.minProductSize = innerLength;
		config.maxProductSize = modifiedSequence.length();
		
		// Change number of primers to return
		config.maxNumPrimersToReturn = 1;
		
		// Get the primer
		return designBestPrimer(config, modifiedSequence, pathPrimer3core);
		
	}
	
	/**
	 * Design primers using primer3
	 * @param config Primer3 configuration
	 * @param seq Sequence
	 * @param pathPrimer3core primer3core executable
	 * @return Collection of primers
	 * @throws IOException
	 */
	public static Collection<PrimerPair> designPrimers(Primer3Configuration config, String seq, String pathPrimer3core) throws IOException {
		Sequence sequence = new Sequence("");
		sequence.setSequenceBases(seq);
		return designPrimers(config, sequence, pathPrimer3core);
	}
	
	/**
	 * Design primers that amplify the entire sequence using primer3 and get the one with the lowest penalty
	 * @param seq Sequence
	 * @param pathPrimer3core primer3core executable
	 * @return Primer with lowest penalty whose product is the entire sequence, or null if there are no valid primers
	 * @throws IOException
	 */
	public static PrimerPair designBestPrimersFullSequence(Sequence seq, String pathPrimer3core) throws IOException {
		Primer3Configuration config = Primer3ConfigurationFactory.getSpecificLengthPCRConfiguration(seq.getLength());
		return designBestPrimer(config, seq, pathPrimer3core);
	}
	
	/**
	 * Design primers using primer3
	 * @param config Primer3 configuration
	 * @param seq Sequence
	 * @param pathPrimer3core primer3core executable
	 * @return Collection of primers
	 * @throws IOException
	 */
	public static Collection<PrimerPair> designPrimers(Primer3Configuration config, Sequence seq, String pathPrimer3core) throws IOException {
		Primer3IO p3io = new Primer3IO(pathPrimer3core);
		p3io.startPrimer3Communications();
		Primer3SequenceInputTags p3sit = new Primer3SequenceInputTags(seq);	
		p3sit.setPrimerSequenceId(seq.getId());
		Collection<PrimerPair> pp = p3io.findPrimerPair(p3sit, config);
		p3io.endPrimer3Communications();
		TreeSet<PrimerPair> rtrn=new TreeSet<PrimerPair>();
		for(PrimerPair pair: pp){
			if(pair.getLeftPrimer()!=null && pair.getRightPrimer()!=null){rtrn.add(pair);}
		}
		return rtrn;
	}
	/**
	 * Design primers using primer3 and get the one with the lowest penalty
	 * @param config Primer3 configuration
	 * @param seq Sequence
	 * @param pathPrimer3core primer3core executable
	 * @return Primer with lowest penalty, or null if there are no primers
	 * @throws IOException
	 */
	public static PrimerPair designBestPrimer(Primer3Configuration config, String seq, String pathPrimer3core) throws IOException {
		Sequence sequence = new Sequence("");
		sequence.setSequenceBases(seq);
		return designBestPrimer(config, sequence, pathPrimer3core);
	}
	
	/**
	 * Design primers using primer3 and get the one with the lowest penalty
	 * @param config Primer3 configuration
	 * @param seq Sequence
	 * @param pathPrimer3core primer3core executable
	 * @return Primer with lowest penalty, or null if there are no primers
	 * @throws IOException
	 */
	public static PrimerPair designBestPrimer(Primer3Configuration config, Sequence seq, String pathPrimer3core) throws IOException {
		Collection<PrimerPair> allPrimers = designPrimers(config, seq, pathPrimer3core);
		if(allPrimers.isEmpty()) {
			return null;
		}
		float minPenalty = Float.MAX_VALUE;
		PrimerPair bestPrimer = null;
		for(PrimerPair primer : allPrimers) {
			if(primer.getPrimerPairPenalty() < minPenalty) {
				minPenalty = primer.getPrimerPairPenalty();
				bestPrimer = primer;
			}
		}
		return bestPrimer;
	}
	
	public static Collection<PrimerPair> designSyntheticPrimers(String seq, int numDesigns, int primerSize, int seqSize, String pathPrimer3core, double optimalMeltingTemp) throws IOException {
		Primer3Configuration best = Primer3ConfigurationFactory.getSyntheticConfiguration(optimalMeltingTemp);
		best.minPrimerSize = primerSize;
		best.maxPrimerSize = primerSize;
		best.optimalPrimerSize = primerSize;
        best.minProductSize = primerSize + 1;
        best.maxProductSize = seqSize;
		
		if(numDesigns>0){best.setPrimerNumReturn(numDesigns);}
				
		Primer3IO p3io = new Primer3IO(pathPrimer3core);
		p3io.startPrimer3Communications();
			
		
				
		Sequence mRNA=new Sequence("Gene");
		mRNA.setSequenceBases(seq);
		
		Primer3SequenceInputTags p3sit = new Primer3SequenceInputTags(mRNA);	
		p3sit.setPrimerSequenceId(mRNA.getId());
		
				
		Collection<PrimerPair> pp = p3io.findPrimerPair(p3sit, best);
		

		p3io.endPrimer3Communications();
		
		TreeSet<PrimerPair> rtrn=new TreeSet<PrimerPair>();
		for(PrimerPair pair: pp){
			if(pair.getLeftPrimer()!=null && pair.getRightPrimer()!=null){rtrn.add(pair);}
		}
		
		return rtrn;
	}


	public static Collection<PrimerPair> designSyntheticPrimers(String seq, String sequencePrimer, String sequencePrimerRevComp, int numDesigns, String pathPrimer3core, double optimalMeltingTemp) throws Exception {
		Primer3Configuration best = Primer3ConfigurationFactory.getSyntheticConfiguration2(optimalMeltingTemp);
		
		if(numDesigns>0){best.setPrimerNumReturn(numDesigns);}
				
		Primer3IO p3io = new Primer3IO(pathPrimer3core);
		p3io.startPrimer3Communications();
			
		Sequence mRNA=new Sequence("Primer");
		mRNA.setSequenceBases(seq);
		Primer3SequenceInputTags p3sit = new Primer3SequenceInputTags(mRNA);	
		p3sit.setPrimerSequenceId(mRNA.getId());
		if(sequencePrimer!=null && !sequencePrimer.isEmpty()){p3sit.setPrimerLeftInput(sequencePrimer);}
		if(sequencePrimerRevComp!=null && !sequencePrimerRevComp.isEmpty()){p3sit.setPrimerRightInput(sequencePrimerRevComp);}
		
		Collection<PrimerPair> pp = p3io.findPrimerPair(p3sit, best);
		

		p3io.endPrimer3Communications();
		
		TreeSet<PrimerPair> rtrn=new TreeSet<PrimerPair>();
		for(PrimerPair pair: pp){
			if(pair.getLeftPrimer()!=null && pair.getRightPrimer()!=null){rtrn.add(pair);}
		}
		
		return rtrn;
	}
	
	private static final String SYNTHETIC_CONFIG_NAME = "synthetic";
	private static final String QPCR_CONFIG_NAME = "qPCR";
	private static final String RAP_QPCR_CONFIG_NAME = "RAPqPCR";
	private static final String DELETION_PLASMID_CONFIG_NAME = "deletion_plasmid";
	
	private static final String[] CONFIG_NAMES = {SYNTHETIC_CONFIG_NAME,QPCR_CONFIG_NAME,RAP_QPCR_CONFIG_NAME,DELETION_PLASMID_CONFIG_NAME};
	
	private static String configNamesList() {
		StringBuilder sb = new StringBuilder(CONFIG_NAMES[0]);
		for(int i = 1; i < CONFIG_NAMES.length; i++) {
			sb.append(", " + CONFIG_NAMES[i]);
		}
		return sb.toString();
	}
	
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-s","Fasta file of sequences to design primers against",true);
		p.addStringArg("-r","Optional file of region coordinates to design primers against (format: sequence_name start_pos end_pos)",false,null);
		p.addStringArg("-c","Primer3 configuration name (" + configNamesList() + ")",true);
		p.addBooleanArg("-rc", "Design primers against antisense strand", true);
		p.addStringArg("-o", "Outfile", true);
		p.addStringArg("-p3c", "Primer3core executable", true);
		
		p.parse(args);
		
		String inputfile = p.getStringArg("-s");
		String regions = p.getStringArg("-r");
		String config = p.getStringArg("-c");
		String outfile = p.getStringArg("-o");
		boolean rc = p.getBooleanArg("-rc");
		String primer3core = p.getStringArg("-p3c");
		
		
		boolean configOk = false;
		for(int i=0; i<CONFIG_NAMES.length; i++) {
			if(config.equals(CONFIG_NAMES[i])) configOk = true;
		}
		
		// validate configuration name
		if(!configOk) {
			logger.info("\nValid Primer3 configuration names:");
			for(int i=0; i<CONFIG_NAMES.length; i++) {
				logger.info(CONFIG_NAMES[i]);
			}
			logger.info("");
			throw new IllegalArgumentException();
		}
		
		// the sequences to design primers for
		ArrayList<Sequence> seqs = new ArrayList<Sequence>();
		
		// the input sequences
		FastaSequenceIO fsio = new FastaSequenceIO(inputfile);
		List<Sequence> inputseqs = fsio.loadAll();
		
		// if regions are not specified, just keep the whole sequences
		if(regions == null) {
			for(Sequence seq : inputseqs) {
				if(rc) seqs.add(seq.getAntisense());
				else seqs.add(seq);
			}
		}
		
		// if regions are specified, extract the regions
		if(regions != null) {
			
			FileReader r = new FileReader(regions);
			BufferedReader b = new BufferedReader(r);
			
			StringParser s = new StringParser();
			
			while(b.ready()) {
				
				String line = b.readLine();
				s.parse(line);
				if(s.getFieldCount() != 3 && s.getFieldCount() != 4) {
					b.close();
					throw new IllegalArgumentException("Region line is not valid:\n" + line + "\nFormat: sequence_name start_pos end_pos <optional ID>");
				}
				String seqname = s.asString(0);
				int start = s.asInt(1);
				int end = s.asInt(2);
				
				Collection<String> seqnames = new ArrayList<String>();
				seqnames.add(seqname);
				
				Collection<Sequence> theseqs = fsio.extractRecords(seqnames);
				Iterator<Sequence> iter = theseqs.iterator();
				Sequence theseq = iter.next();
				
				Sequence subseq = theseq.getSubSequence(theseq.getId() + ":" + start + "-" + end, start, end);
				if(s.getFieldCount() == 4) {
					subseq.setId(s.asString(3));
				}
				if(rc) seqs.add(subseq.getAntisense());
				else seqs.add(subseq);
				
			}
			
			b.close();
			
		}
		
		// set up primer3
		Primer3Configuration primer3config = new Primer3Configuration();
		
		if(config.equals(SYNTHETIC_CONFIG_NAME)) primer3config = Primer3ConfigurationFactory.getSyntheticConfiguration(60);
		if(config.equals(RAP_QPCR_CONFIG_NAME)) primer3config = Primer3ConfigurationFactory.getRAPqPCRConfiguration();
		if(config.equals(QPCR_CONFIG_NAME)) primer3config = Primer3ConfigurationFactory.getQpcrConfiguration();
		// If plasmids, will get separate configuration for each plasmid
		if(config.equals(DELETION_PLASMID_CONFIG_NAME)) primer3config = null;
		
		// output file
		FileWriter writer = new FileWriter(outfile);
		
		String header = "primer_ID\t";
		header += "left_primer\t";
		header += "right_primer\t";
		header += "left_primer_TM\t";
		header += "right_primer_TM\t";
		header += "primer_pair_penalty";
		writer.write(header + "\n");

		for(Sequence seq : seqs) {
			
			PrimerPair primer = PcrPrimerDesigner.designBestPrimer(primer3config, seq, primer3core);
			String lineToWrite = seq.getId() + "\t";
			if(primer == null) {
				lineToWrite += "NO_PRIMERS";
				writer.write(lineToWrite + "\n");
				continue;
			}
			lineToWrite += primer.getLeftPrimer().toUpperCase() + "\t";
			lineToWrite += primer.getRightPrimer().toUpperCase() + "\t";
			lineToWrite += primer.getLeftPrimerTM() + "\t";
			lineToWrite += primer.getRightPrimerTM() + "\t";
			lineToWrite += primer.getPrimerPairPenalty();
			writer.write(lineToWrite + "\n");

			
		}
		
		writer.close();
		
		logger.info("");
		logger.info("All done.");
		
	}
	
	
	
}
