package nextgen.editing.crispr;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Pattern;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.Annotation.Strand;
import nextgen.core.annotation.Gene;

import org.apache.log4j.Logger;

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
		public String fastaSeqId;
		List<SAMRecord> matches;

		public CRISPRTarget(SequenceRegion targetRegion) {
			sequence = targetRegion.getSequenceBases();
			start = targetRegion.getStart();
			matches = new ArrayList<SAMRecord>();
		}

		public String toString() {
			StringBuilder sb = new StringBuilder(String.valueOf(start));
			sb.append("\t").append(sequence)
				.append("\t").append(genomicStart)
				.append("\t").append(getDistanceToTarget())
				.append("\t").append(orientation)
				.append("\t").append(gene != null ? gene.getOrientation() : ".")
				.append("\t").append(matches.size());
			
			for(SAMRecord r : matches) {
				sb.append("\t")
					.append(r.getReferenceName())
					.append("_")
					.append(r.getAlignmentStart())
					.append("(")
					.append(r.getReadNegativeStrandFlag() ? "-" : "+")
					.append(")");
			}
			return sb.toString();
		}
		
		void setGenomicLocation(Annotation containingRegion) {
			orientation = containingRegion.getOrientation();
			//genomicStart = containingRegion.getStart() + start;
			genomicStart =   containingRegion.getStart() + ( containingRegion.isNegativeStrand() ? (containingRegion.length() -( sequence.length() + start)) : start);
			logger.debug("Promoter region: " + containingRegion.toShortBED() + "("+containingRegion.getOrientation() + ") match start: " + start + " genomicStart: " + genomicStart);
		}

		public int getDistanceToTarget() {
			if(gene != null)  {
				return gene.isNegativeStrand() ?   	genomicStart -  gene.getOrientedStart() : gene.getOrientedStart() - genomicStart;
			} else {
				return 0;
			}
			
		}

		public void setGene(Gene gene) {
			this.gene = gene;
			
		}

		public void setFastaSeqId(String fastaId) {
			this.fastaSeqId = fastaId;
			
		}

		public String getFastaSeqId() {
			return fastaSeqId;
		}

		public boolean isSequence(String seq) {
			return seq.equals(this.sequence);
		}

		public void addMatch(SAMRecord aln) {
			matches.add(aln);
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
			"\n\tDesign. Get putative CRISPR matches to sequence: \n\t\t-genes <Annotation file in BED format> \n\t\t-num <Number of desired targets per sequence> "+
			"\n\t\t-sequenceDir <Directory of the genomic sequence. Assumes each chromosome in its own directory <sequenceDIr>/N/chrN.fa> " +
			"\n\t\t-promoterStart <In bases before the TSS> \n\t\t-promoterEnd <In bases past the TSS>" +
			"\n\t\tIf  you wish to invoke Bowtie to test for other possible matches by adding the following paramters: " +
			"\n\t\t-bowtieBuild <e.g. full path to the bowtie build> -bowtieExcutable <path to the Bowtie executable > " +
			"\n\tDesignFromFasta. Get putative CRISPR matches to sequence: \n\t\t-in <FASTA Sequences> \n\t\t-num <Number of desired targets>" +
			"\n\t\tIf  you wish to invoke Bowtie to test for other possible matches by adding the following paramters: " +
			"\n\t\t-bowtieBuild <e.g. full path to the bowtie build> -bowtieExcutable <path to the Bowtie executable > " +
			"\n";

	/**
	 * @param args
	 * @throws SearchException 
	 * @throws NumberFormatException 
	 * @throws IOException 
	 * @throws InterruptedException 
	 */
	public static void main(String[] args) throws NumberFormatException, SearchException, IOException, InterruptedException {
		ArgumentMap argMap = CLUtil.getParameters(args, USAGE, "Design");
		CRISPRDesigner designer = new CRISPRDesigner();
		if ("DesignFromFasta".equalsIgnoreCase(argMap.getTask())) {
			FastaSequenceIO fsio = new FastaSequenceIO(argMap.getInput());
			int numToDesign = argMap.getInteger("num");
			List<Sequence> all = fsio.loadAll();

			LinkedHashMap<String, List<CRISPRTarget>> result = new LinkedHashMap<String, List<CRISPRTarget>>();
			for (Sequence seq : all) {
				List<CRISPRTarget> targets = designer.design(seq, numToDesign);
				if (targets.size() < numToDesign ) {
					seq.reverse();
					List<CRISPRTarget> reverseTargets = designer.design(seq, numToDesign - targets.size()); 
					targets.addAll(reverseTargets);
				}
				result.put(seq.getId(), targets);
			}

			if(argMap.isPresent("bowtieBuild")) {
				String bowtieBuild = argMap.getMandatory("bowtieBuild");
				String bowtie     =  argMap.getMandatory("bowtieExecutable");
				checkOfTargets(result, bowtieBuild, bowtie);
				
			}

			writeDesign(argMap, result);

		} else if ("Design".equalsIgnoreCase(argMap.getTask()) ) {
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

					if (targets.size() < numToDesign ) {
						promoterRegion.reverse();
						List<CRISPRTarget> reverseTargets = designer.design(promoterRegion, numToDesign - targets.size()); 
						promoterRegion.setOrientation(promoterRegion.isNegativeStrand() ?  Strand.POSITIVE : Strand.NEGATIVE);
						adjustPositions(reverseTargets, promoterRegion, g);
						targets.addAll(reverseTargets);
					}

					result.put(g.getName(), targets);
				}

			}
			
			if(argMap.isPresent("bowtieBuild")) {
				String bowtieBuild = argMap.getMandatory("bowtieBuild");
				String bowtie     =  argMap.getMandatory("bowtieExecutable");
				checkOfTargets(result, bowtieBuild, bowtie);
				
			}

			writeDesign(argMap, result);
		}else {
			System.err.println(USAGE);
		}

	}

	protected static void checkOfTargets(
			LinkedHashMap<String, List<CRISPRTarget>> result,
			String bowtieBuild, String bowtie) throws IOException,
			InterruptedException {
		String fastaFile = writeTargetsAsFasta(result );
		String alignmentFile = runBowtie (bowtie, bowtieBuild, fastaFile);
		if( alignmentFile != null ) {
			updateResultsWithAlignment(alignmentFile, result);
		}
	}

	private static void updateResultsWithAlignment(String alignmentFile, LinkedHashMap<String, List<CRISPRTarget>> result) {
		File alignmentFileFile = new File(alignmentFile);
		if(!alignmentFileFile.exists()) {
			logger.info("Alignmnent file " + alignmentFile + " did not exists. Bowtie must not have succeeded");
			return;
		} 
		
		SAMFileReader samReader = new SAMFileReader(alignmentFileFile);
		SAMRecordIterator rIt = samReader.iterator();
		while (rIt.hasNext()) {
			SAMRecord aln = rIt.next();
			String readName = aln.getReadName();
			String targetGene = readName.split("___")[0];
			String seq = aln.getReadString();
			
			List<CRISPRTarget> geneTargets = result.get(targetGene);
			
			for (CRISPRTarget t : geneTargets) {
				if (t.isSequence(seq) || t.isSequence(Sequence.complement(seq))) {
					t.addMatch(aln);
					break;
				}
			}

		}
		rIt.close();
		samReader.close();
		
	}

	private static String runBowtie(String bowtie, String bowtieBuild, String fastaFile) throws IOException, InterruptedException {
		String samTmpFile = fastaFile+".bowtie2.sam";
		String bowtieCmd = bowtie + " --local -f -k 10  --very-sensitive-local -L 9 -N 1 " +  bowtieBuild + " -U  " + fastaFile + "  " +  samTmpFile ;
		logger.debug("starting bowtie: " + bowtieCmd);
		Process p = Runtime.getRuntime().exec(bowtieCmd);
		
		logger.debug("Waiting for bowtie to finish");
		int time = p.waitFor();
		logger.debug("Bowtie finished " + time + " " + p.exitValue());
		

		return p.exitValue() == 0 ? samTmpFile : null;
	}

	private static String writeTargetsAsFasta(LinkedHashMap<String, List<CRISPRTarget>> result) throws IOException {
		String out = "tmp."+System.currentTimeMillis()+".fa" ;
		FastaSequenceIO fsio = new FastaSequenceIO(out);
		for (String targetGene : result.keySet()) {
			List<CRISPRTarget> constructs = result.get(targetGene);
			int n = 0;
			for(CRISPRTarget c : constructs) {
				String fastaId = targetGene + "___" + n + "_" + String.valueOf(c.getDistanceToTarget() ).replace("-", "neg");
				c.setFastaSeqId(fastaId);
				Sequence cS = new Sequence(c.getFastaSeqId() );
				cS.setSequenceBases(c.sequence.toUpperCase());
				n++;
				fsio.append(cS, out);
			}
		}

		return out;
	}




}
