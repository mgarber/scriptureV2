package nextgen.editing;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Pattern;

import org.apache.log4j.Logger;

import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.Annotation.Strand;
import nextgen.core.annotation.Gene;

import broad.core.motif.SearchException;
import broad.core.motif.SequenceMotif;
import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;
import broad.core.sequence.SequenceRegion;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;
import broad.pda.annotation.BEDFileParser;



public class CRISPRDesigner {

	static Logger logger = Logger.getLogger(CRISPRDesigner.class.getName());
	public static final SequenceMotif CRISPR_TARGET_MOTIF = new SequenceMotif( Pattern.compile("G[A,C,G,T]{20}GG"));

	public List<CRISPRTarget> design (Sequence seq, int number) {
		List<SequenceRegion> targetRegions = CRISPR_TARGET_MOTIF.match(seq);
		List<CRISPRTarget> targets = new ArrayList<CRISPRTarget>(number);

		int targetNum = 0;
		for (SequenceRegion targetRegion : targetRegions) {

			CRISPRTarget target = new CRISPRTarget(targetRegion);
			targets.add(target);
			targetNum++;
			if(targetNum >= number) {
				break;
			}
		}

		return targets;
	}

	public static class CRISPRTarget {
		String sequence;
		int    start;
		int    genomicStart;
		Strand orientation;
		private Gene gene;

		public CRISPRTarget(SequenceRegion targetRegion) {
			sequence = targetRegion.getSequenceBases();
			start = targetRegion.getStart();
		}

		public String toString() {
			StringBuilder sb = new StringBuilder(String.valueOf(start));
			sb.append("\t").append(sequence)
				.append("\t").append(genomicStart)
				.append("\t").append(getDistanceToTarget())
				.append("\t").append(orientation)
				.append("\t").append(gene.getOrientation());
			return sb.toString();
		}
		
		void setGenomicLocation(Annotation containingRegion) {
			orientation = containingRegion.getOrientation();
			genomicStart = containingRegion.getStart() + start;// ( containingRegion.isNegativeStrand() ?  (containingRegion.length() - start) : start); 
		}

		public int getDistanceToTarget() {
			return gene.isNegativeStrand() ?   genomicStart -  gene.getOrientedStart() : gene.getOrientedStart() - genomicStart;
			
		}

		public void setGene(Gene gene) {
			this.gene = gene;
			
		}

	}
	
	private static void adjustPositions(List<CRISPRTarget> targets, SequenceRegion promoterRegion, Gene gene) {
		for (CRISPRTarget target : targets) {
			target.setGenomicLocation(promoterRegion);
			target.setGene(gene);
		}
	}

	private static Sequence getChromosomeSequence(String chr, File sequenceDir) throws IOException {
		String chrNum = chr.replace("chr", "");
		String filePath = sequenceDir + "/" + chrNum + "/" + chr + ".fa";
		FastaSequenceIO fsio = new FastaSequenceIO(filePath);
		return fsio.loadAll().get(0);
	}

	private static void writeDesign(ArgumentMap argMap,
			LinkedHashMap<String, List<CRISPRTarget>> result)
			throws IOException {
		BufferedWriter bw = argMap.getOutputWriter();
		for (String sid : result.keySet()) {
			bw.write(sid);
			bw.newLine();
			List<CRISPRTarget> targets = result.get(sid);
			for(CRISPRTarget t : targets) {
				bw.write(t.toString());
				bw.newLine();
			}
		}
		bw.close();
	}

	public static String USAGE = "Usage: CRISPRDesigner TASK=<task> <task_args>\n" +
			"\tTasks:\n" +
			"\t\tDesign. Get putative CRISPR matches to sequence: \n\t\t-in <Annotation file in BED format> \n\t\t-num <Number of desired targets per sequence> "+
			"\n\t\t-sequenceDir <Directory of the genomic sequence. Assumes each chromosome in its own directory <sequenceDIr>/N/chrN.fa> " +
			"\n\t\t-promoterStart <In bases before the TSS> \n\t\t-promoterEnd <In bases past the TSS>" +
			"\t\tDesign2. Get putative CRISPR matches to sequence: \n\t\t-in <FASTA Sequence> \n\t\t-num <Number of desired targets>\n" +
			"\n";

	/**
	 * @param args
	 * @throws SearchException 
	 * @throws NumberFormatException 
	 * @throws IOException 
	 */
	public static void main(String[] args) throws NumberFormatException, SearchException, IOException {
		ArgumentMap argMap = CLUtil.getParameters(args, USAGE, "Design2");
		CRISPRDesigner designer = new CRISPRDesigner();
		if ("Design".equals(argMap.getTask())) {
			FastaSequenceIO fsio = new FastaSequenceIO(argMap.getInput());
			int numToDesign = argMap.getInteger("num");
			List<Sequence> all = fsio.loadAll();

			LinkedHashMap<String, List<CRISPRTarget>> result = new LinkedHashMap<String, List<CRISPRTarget>>();
			for (Sequence seq : all) {
				List<CRISPRTarget> targets = designer.design(seq, numToDesign);
				result.put(seq.getId(), targets);
			}

			writeDesign(argMap, result);

		} else if ("Design2".equalsIgnoreCase(argMap.getTask()) ) {
			String annotationFile = argMap.getMandatory("genes");
			File sequenceDir = new File(argMap.getMandatory("sequenceDir"));
			Map<String, Collection<Gene>> genes = BEDFileParser.loadDataByChr(new File(annotationFile));
			int tssMinus = argMap.getInteger("promoterStart");
			int tssPlus = argMap.getInteger("promoterEnd");
			int numToDesign = argMap.getInteger("num");
			LinkedHashMap<String, List<CRISPRTarget>> result = new LinkedHashMap<String, List<CRISPRTarget>>();

			for (String chr : genes.keySet()) {
				Sequence chromosomeSequence  = getChromosomeSequence(chr, sequenceDir);
				Collection<Gene> chrGenes = genes.get(chr);
				for (Gene g : chrGenes) {
					int start = g.getOrientedStart();
					SequenceRegion promoterRegion = new SequenceRegion(chromosomeSequence.getId(), chr, g.isNegativeStrand() ? start - tssPlus  : start - tssMinus, g.isNegativeStrand() ? start + tssMinus : start + tssPlus);
					promoterRegion.setOrientation(g.getOrientation());
					chromosomeSequence.getRegion(promoterRegion);


					List<CRISPRTarget> targets = designer.design(promoterRegion, numToDesign);
					adjustPositions(targets, promoterRegion, g);

					if(targets.size() < numToDesign ) {
						promoterRegion.reverse();
						List<CRISPRTarget> reverseTargets = designer.design(promoterRegion, numToDesign - targets.size()); 
						promoterRegion.setOrientation(!promoterRegion.isNegativeStrand() ?   "+" : "-");
						adjustPositions(reverseTargets, promoterRegion, g);
						targets.addAll(reverseTargets);
					}

					result.put(g.getName(), targets);
				}

			}
			writeDesign(argMap, result);


		}else {
			System.err.println(USAGE);
		}

	}




}
