package broad.core.primer3;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.*;

import nextgen.core.annotation.Gene;

import broad.core.annotation.GenomicAnnotation;
import broad.core.error.ParseException;
import broad.core.parser.CommandLineParser;
import broad.core.parser.StringParser;
import broad.core.primer3.Primer3SequenceInputTags.SequenceRegionCoordinates;
import broad.core.sequence.Extractor;
import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;
import broad.core.sequence.SequenceRegion;
import broad.pda.datastructures.Alignments;


public class qPCRPrimerDesigner  {
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
	
	private Primer3Configuration [] configurations;
	private Primer3IO p3io;
	private ArrayList<SequenceRegionCoordinates> exludedAnnotations = new ArrayList<SequenceRegionCoordinates>();
	private boolean storeExpectedProd;

	private int buffer;

	public qPCRPrimerDesigner(Primer3Configuration [] confs, Primer3IO p3io, int buffer) {
		this.configurations = confs;
		this.p3io = p3io; 
		this.buffer = buffer;
	}
	
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
		System.out.println("Designing primers for Sequence " + sequence.getId() + " \nBASES:" + sequence.getSequenceBases() +
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
					System.out.println("Primer has more primers? "+pp.hasPrimers() + " # confs " + configurations.length + " doing configuration " + confId + " minProductSize " + conf.minProductSize + " maxProdSize " + conf.maxProductSize);
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
				System.out.println(sw);
			}

		}
		return pr;
	}*/
	

	private void addExcludedAnnotations(List<? extends GenomicAnnotation> excludeList) {
		Iterator<? extends GenomicAnnotation> it = excludeList.iterator();
		while(it.hasNext()) {
			GenomicAnnotation annot = it.next();
			this.exludedAnnotations.add(new SequenceRegionCoordinates((int)annot.getStart(), (int)annot.getEnd()));
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
		System.out.println("Max probe size " + best.maxProductSize + " min " + best.minProductSize+" "+target.toString());
		
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
		//System.err.println(pp.getPrimerSequences().get(0));
			
		
		
		return pp;
	}*/
	
	

	/*public static PrimerPair designPrimers(Chromosome chr, Alignments align, boolean repeatMask)throws Exception{
		
		SequenceRegion target = new SequenceRegion("chr"+chr.getSymbol());
		target.setRegionStart(align.getStart());
		target.setRegionEnd(align.getEnd());
		
		Primer3Configuration best = Primer3ConfigurationFactory.getqRTPCRConfiguration();
		System.out.println("Max probe size " + best.maxProductSize + " min " + best.minProductSize+" "+target.toString());
		
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
		//System.err.println(pp.getPrimerSequences().get(0));
			
		
		
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
		//System.err.println(targetList);
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
		System.err.println(mRNA.getSequenceBases());
		
		Primer3SequenceInputTags p3sit = new Primer3SequenceInputTags(mRNA);	
		p3sit.setPrimerSequenceId(gene.getName());
		
		ArrayList<Integer> splicePositions=gene.getSpliceJunctionCoordinates();
		//System.err.println(targetList);
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
					System.err.println(gene.getAlignment());
					primers.add(pp);
				}
				//return pp;
			}
			p3io.endPrimer3Communications();
			return primers;
		//}catch(NullPointerException ex){return null;}
	}*/
	
	//TODO should add all junctions as sequence regions
	public static Set<PrimerPair> designJenPrimers(Sequence chr, Gene gene, boolean repeatMask, int numDesigns)throws Exception{
		
		
		Primer3Configuration best = Primer3ConfigurationFactory.getJenRTPCRConfiguration();
		
		Primer3Configuration [] configurations = {best};				
		Primer3IO p3io = new Primer3IO();
		p3io.startPrimer3Communications();
			
		
				
		qPCRPrimerDesigner pd = new qPCRPrimerDesigner(configurations, p3io, 0);
		ArrayList<Sequence> primerSequences = new ArrayList<Sequence>();
			
		Sequence mRNA=new Sequence(gene.getName());
		mRNA.setSequenceBases(gene.getSequence(chr, repeatMask, false));
		//System.err.println(mRNA.getSequenceBases());
		
		Primer3SequenceInputTags p3sit = new Primer3SequenceInputTags(mRNA);	
		p3sit.setPrimerSequenceId(gene.getName());
		
		/**Add exon junction sequence targets**/
		ArrayList<Integer> splicePositions=gene.getSpliceJunctionCoordinates();
		System.err.println(gene.getName()+" "+gene.getNumExons());
		//System.err.println(targetList);
		ArrayList<SequenceRegionCoordinates> targetList=new ArrayList();
		for(Integer spliceJunction: splicePositions){
			SequenceRegionCoordinates region=new SequenceRegionCoordinates(spliceJunction-1, spliceJunction+1);
			targetList.add(region);
		}
		p3sit.addTargets(targetList);
		//add targets p3sit.
		/***/
		
		//try{
			Set<PrimerPair> primers=new HashSet();
			
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
		System.err.println(mRNA.getSequenceBases());
		
		Primer3SequenceInputTags p3sit = new Primer3SequenceInputTags(mRNA);	
		p3sit.setPrimerSequenceId(gene.getName());
		
		ArrayList<Integer> splicePositions=gene.getSpliceJunctionCoordinates();
		//System.err.println(targetList);
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
					System.err.println(gene.getAlignment());
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
		System.err.println(mRNA.getSequenceBases());
		
		Primer3SequenceInputTags p3sit = new Primer3SequenceInputTags(mRNA);	
		p3sit.setPrimerSequenceId(gene.getName());
		
		//ArrayList<SequenceRegionCoordinates> targetList=gene.getSpliceJunctionCoordinates(chr, repeatMask, 0);
		//System.err.println(targetList);
		//p3sit.addTargets(targetList);
		//add targets p3sit.
		
		//try{
			Set<PrimerPair> primers=new HashSet();
			for(int i=0; i<numDesigns; i++){
				PrimerPair pp = p3io.findPrimerPair(p3sit, best);
				if(pp!=null &&pp.getLeftPrimer()!=null && pp.getRightPrimer()!=null){
					p3sit.addExcludedRegion(new SequenceRegionCoordinates(pp.getLeftPrimerPosition(), pp.getLeftPrimerPosition()+pp.getLeftPrimer().toCharArray().length));
					p3sit.addExcludedRegion(new SequenceRegionCoordinates(pp.getRightPrimerPosition()-pp.getRightPrimer().toCharArray().length, pp.getRightPrimerPosition()));
					System.err.println(gene.getAlignment());
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
		System.err.println(mRNA.getSequenceBases());
		
		Primer3SequenceInputTags p3sit = new Primer3SequenceInputTags(mRNA);	
		p3sit.setPrimerSequenceId(gene.getName());
		
		
		
		//try{
			Set<PrimerPair> primers=new HashSet();
			for(int i=0; i<numDesigns; i++){
				PrimerPair pp = p3io.findPrimerPair(p3sit, best);
				if(pp!=null &&pp.getLeftPrimer()!=null && pp.getRightPrimer()!=null){
					p3sit.addExcludedRegion(new SequenceRegionCoordinates(pp.getLeftPrimerPosition(), pp.getLeftPrimerPosition()+pp.getLeftPrimer().toCharArray().length));
					p3sit.addExcludedRegion(new SequenceRegionCoordinates(pp.getRightPrimerPosition()-pp.getRightPrimer().toCharArray().length, pp.getRightPrimerPosition()));
					System.err.println(gene);
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
		System.err.println(mRNA.getSequenceBases());
		
		Primer3SequenceInputTags p3sit = new Primer3SequenceInputTags(mRNA);	
		p3sit.setPrimerSequenceId(gene.getName());
		
		
		
		//try{
			Set<PrimerPair> primers=new HashSet();
			for(int i=0; i<numDesigns; i++){
				PrimerPair pp = p3io.findPrimerPair(p3sit, best);
				if(pp!=null &&pp.getLeftPrimer()!=null && pp.getRightPrimer()!=null){
					p3sit.addExcludedRegion(new SequenceRegionCoordinates(pp.getLeftPrimerPosition(), pp.getLeftPrimerPosition()+pp.getLeftPrimer().toCharArray().length));
					p3sit.addExcludedRegion(new SequenceRegionCoordinates(pp.getRightPrimerPosition()-pp.getRightPrimer().toCharArray().length, pp.getRightPrimerPosition()));
					System.err.println(gene);
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
					System.err.println(mRNA.getId());
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
				System.out.println(pp);
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


	public static Collection<PrimerPair> designIntronPrimers(Sequence chr, Gene gene, boolean repeatMask, String sequencePrimer, String sequencePrimerRevComp, int numDesigns, int min3, int min5) throws Exception {
		Primer3Configuration best = Primer3ConfigurationFactory.getJenRTPCRConfiguration();
		
		if(numDesigns>0){best.setPrimerNumReturn(numDesigns);}
		if(min3>0){best.setPrimerMin3PrimeOverlapOfJunction(min3);}
		if(min5>0){best.setPrimerMin5PrimeOverlapOfJunction(min5);}
		
		Primer3Configuration [] configurations = {best};				
		Primer3IO p3io = new Primer3IO();
		p3io.startPrimer3Communications();
			
		
				
		qPCRPrimerDesigner pd = new qPCRPrimerDesigner(configurations, p3io, 0);
		ArrayList<Sequence> primerSequences = new ArrayList<Sequence>();
			
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
		
		TreeSet<PrimerPair> rtrn=new TreeSet();
		for(PrimerPair pair: pp){
			if(pair.getLeftPrimer()!=null && pair.getRightPrimer()!=null){rtrn.add(pair);}
		}
		
		return rtrn;
	}
	
	
	public static Collection<PrimerPair> designRACEPrimers(Sequence mRNA, String leftPrimer, String rightPrimer, int numDesigns) throws Exception {
		Primer3Configuration best = Primer3ConfigurationFactory.getJenRTPCRConfiguration();
				
		Primer3Configuration [] configurations = {best};				
		Primer3IO p3io = new Primer3IO();
		p3io.startPrimer3Communications();
			
		
				
		qPCRPrimerDesigner pd = new qPCRPrimerDesigner(configurations, p3io, 0);
		ArrayList<Sequence> primerSequences = new ArrayList<Sequence>();
			
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

	public static Collection<PrimerPair> designPCRPrimers(Sequence chr, Gene gene, boolean repeatMask,	int numDesigns, boolean crossJunction) throws Exception {
		Primer3Configuration best = Primer3ConfigurationFactory.getJenRTPCRConfiguration();
		
		if(numDesigns>0){best.setPrimerNumReturn(numDesigns);}
		if(!crossJunction){repeatMask=true;}		
		
		Primer3IO p3io = new Primer3IO();
		p3io.startPrimer3Communications();
			
		
				
		Sequence mRNA=new Sequence(gene.getName());
		mRNA.setSequenceBases(gene.getSequence(chr, repeatMask, false));
		ArrayList<Integer> splicePositions=gene.getSpliceJunctionCoordinates();
		Primer3SequenceInputTags p3sit = new Primer3SequenceInputTags(mRNA);	
		p3sit.setPrimerSequenceId(gene.getName());
		
		if(crossJunction){
			ArrayList<SequenceRegionCoordinates> targetList=new ArrayList<SequenceRegionCoordinates>();
			for(Integer spliceJunction: splicePositions){
				SequenceRegionCoordinates region=new SequenceRegionCoordinates(spliceJunction-1, spliceJunction+1);
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
	
	public static Collection<PrimerPair> designPCRPrimers(Sequence chr, Gene gene, boolean repeatMask,	int numDesigns,  SequenceRegionCoordinates target) throws Exception {
		Primer3Configuration best = Primer3ConfigurationFactory.getJenRTPCRConfiguration();
		
		if(numDesigns>0){best.setPrimerNumReturn(numDesigns);}	
		
		Primer3IO p3io = new Primer3IO();
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
	
	public static Collection<PrimerPair> designCloningPrimers(Sequence geneSequence, int numDesigns, SequenceRegionCoordinates target) throws Exception {
		Primer3Configuration best = Primer3ConfigurationFactory.getLongRangePCRConfiguration();
		
		//if(numDesigns>0){best.setPrimerNumReturn(numDesigns);}
		
		
		Primer3Configuration [] configurations = {best};				
		Primer3IO p3io = new Primer3IO();
		p3io.startPrimer3Communications();
			
		qPCRPrimerDesigner pd = new qPCRPrimerDesigner(configurations, p3io, 0);
				
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
				//System.err.println(mRNA.getId());
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
	
	public static Collection<PrimerPair> designCloningPrimers(Sequence geneSequence, int numDesigns) throws Exception {
		return designCloningPrimers(geneSequence, numDesigns, null);
	}

	public static Collection<PrimerPair> designSyntheticPrimers(String seq, int numDesigns) throws IOException {
		Primer3Configuration best = Primer3ConfigurationFactory.getSyntheticConfiguration();
		
		if(numDesigns>0){best.setPrimerNumReturn(numDesigns);}
				
		Primer3IO p3io = new Primer3IO();
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
	
	public static Collection<PrimerPair> designSyntheticPrimers(String seq, int numDesigns, int primerSize, int seqSize) throws IOException {
		Primer3Configuration best = Primer3ConfigurationFactory.getSyntheticConfiguration();
		best.minPrimerSize = primerSize;
		best.maxPrimerSize = primerSize;
		best.optimalPrimerSize = primerSize;
        best.minProductSize = primerSize + 1;
        best.maxProductSize = seqSize;
		
		if(numDesigns>0){best.setPrimerNumReturn(numDesigns);}
				
		Primer3IO p3io = new Primer3IO();
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


	public static Collection<PrimerPair> designSyntheticPrimers(String seq, String sequencePrimer, String sequencePrimerRevComp, int numDesigns) throws Exception {
		Primer3Configuration best = Primer3ConfigurationFactory.getSyntheticConfiguration2();
		
		if(numDesigns>0){best.setPrimerNumReturn(numDesigns);}
				
		Primer3Configuration [] configurations = {best};				
		Primer3IO p3io = new Primer3IO();
		p3io.startPrimer3Communications();
			
		qPCRPrimerDesigner pd = new qPCRPrimerDesigner(configurations, p3io, 0);
		
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
	
	private static final String[] CONFIG_NAMES = {"synthetic","qPCR","RAPqPCR"};
	
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-s","Input sequences to design primers against",true);
		p.addStringArg("-r","Optional file of region coordinates to design primers against (format: sequence_name start_pos end_pos)",false,null);
		p.addStringArg("-c","Primer3 configuration name",true);
		p.addBooleanArg("-rc", "Design primers against antisense strand", true);
		p.addStringArg("-o", "outfile", true);
		
		p.parse(args);
		
		String inputfile = p.getStringArg("-s");
		String regions = p.getStringArg("-r");
		String config = p.getStringArg("-c");
		String outfile = p.getStringArg("-o");
		boolean rc = p.getBooleanArg("-rc").booleanValue();
		
		boolean configOk = false;
		for(int i=0; i<CONFIG_NAMES.length; i++) {
			if(config.equals(CONFIG_NAMES[i])) configOk = true;
		}
		
		// validate configuration name
		if(!configOk) {
			System.err.println("\nValid Primer3 configurarion names:");
			for(int i=0; i<CONFIG_NAMES.length; i++) {
				System.err.println(CONFIG_NAMES[i]);
			}
			System.err.println();
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
				if(s.getFieldCount() != 3) {
					throw new IllegalArgumentException("Region line is not valid:\n" + line + "\nFormat: sequence_name start_pos end_pos");
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
				if(rc) seqs.add(subseq.getAntisense());
				else seqs.add(subseq);
				
			}
			
		}
		
		// set up primer3
		Primer3Configuration primer3config = new Primer3Configuration();
		
		if(config.equals("synthetic")) primer3config = Primer3ConfigurationFactory.getSyntheticConfiguration();
		if(config.equals("RAPqPCR")) primer3config = Primer3ConfigurationFactory.getRAPqPCRConfiguration();
		if(config.equals("qPCR")) primer3config = Primer3ConfigurationFactory.getQpcrConfiguration();

		Primer3IO p3io = new Primer3IO();
		p3io.startPrimer3Communications();
		
		// output file
		FileWriter writer = new FileWriter(outfile);
		writer.write("Sequence\t" + PrimerPair.getPrimerPairInformationFieldNames() + "\n");
		
		for(Sequence seq : seqs) {
			
			Primer3SequenceInputTags p3sit = new Primer3SequenceInputTags(seq);	
			p3sit.setPrimerSequenceId(seq.getId());
			Collection<PrimerPair> pp = p3io.findPrimerPair(p3sit, primer3config);
			if(pp == null || pp.isEmpty()) {
				writer.write(seq.getId() + "\tNO_PRIMERS\n");
				continue;
			}
			
			for(PrimerPair pair: pp){
				if(pair.getLeftPrimer()!=null && pair.getRightPrimer()!=null) {
					writer.write(seq.getId() + "\t" + pair.getPrimerPairInformation() + "\n");
				}
			}
			
		}
		
		p3io.endPrimer3Communications();
		writer.close();
		
	}
	
	
	
}
