package broad.pda.seq.rap;

import java.util.*;
import java.util.regex.Pattern;
import java.io.*;
import java.nio.file.Files;

import org.apache.commons.collections15.Predicate;
import org.apache.commons.collections15.CollectionUtils;

import broad.core.motif.SequenceMotif;
import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;
import broad.core.sequence.SequenceRegion;
import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.io.IoUtil;
import net.sf.picard.util.Log;
import net.sf.samtools.util.CloseableIterator;
import nextgen.core.annotation.*;
import nextgen.editing.crispr.*;
import nextgen.editing.crispr.score.*;
import nextgen.editing.crispr.predicate.*;
import nextgen.core.capture.filter.PolyBaseFilter;

/**
 * @author engreitz
 * Jesse's CRISPR designer
 */
public class CRISPRDesigner extends CommandLineProgram {
    private static final Log log = Log.getInstance(CRISPRDesigner.class);
    
    public static final SequenceMotif CRISPR_TARGET_MOTIF = new SequenceMotif( Pattern.compile("G[A,C,G,T]{20}GG"));
	
	@Option(doc="Bowtie build")
	public String BOWTIE_BUILD = "/seq/lincRNA/data/mm9.nonrandom";
	
	@Option(doc="Genome fasta file")
	public File GENOME_FASTA = new File("/seq/lincRNA/data/mm9.nonrandom.fa");
	
	@Option(doc="Bed file with annotations to design CRISPRs to")
	public File TARGETS;
	
	@Option(doc="Bitpacked file containing off target sites", optional=true)
	public File OFF_TARGETS = null;

	@Option(doc="Output directory to create")
	public File OUTPUT_DIR;
	
	@Option(doc="Skip initial step of generating guides ... assumes that file exists and reads from there")
	public Boolean SKIP_GENERATION = false;
	
	@Option(doc="Minimum distance between guides in a pair")
	public Integer MIN_DISTANCE = 20;
	
	@Option(doc="Maximum distance between guides in a pair")
	public Integer MAX_DISTANCE = 200;
	
	@Option(doc="Skip scoring ... do this to skip the off target and efficacy scoring steps")
	public Boolean SKIP_SCORING = false;
	

	/**
	 * Stock main method.
	 *
	 * @param args main arguments
	 */
	public static void main(final String[] args) {
		System.exit(new CRISPRDesigner().instanceMain(args));
	}
	

    @Override
    protected int doWork() {
        IoUtil.assertFileIsReadable(TARGETS);
        IoUtil.assertFileIsReadable(GENOME_FASTA);
        if (OFF_TARGETS != null) IoUtil.assertFileIsReadable(OFF_TARGETS);

        if (!OUTPUT_DIR.exists()) OUTPUT_DIR.mkdir();

        File allGuidesFile = new File(OUTPUT_DIR.getPath() + "/allGuides.bed");
		File filteredFile = new File(OUTPUT_DIR.getPath() + "/filteredGuides.bed");
		if (filteredFile.exists() & !SKIP_SCORING) {
			log.error("filteredGuides.bed file already exists but SKIP_SCORING is set to false ... remove file or set SKIP_SCORING to true");
			return 1;
		}
		
		File filteredPairs = new File(OUTPUT_DIR.getPath() + "/filteredPairs.bed");
        
        try {
        	
        	AnnotationList<GuideRNA> guides = null;
        	if (!SKIP_GENERATION) {
        		// Iterate through buffered query annotations to handle large files gracefully
        		CloseableIterator<Annotation> itr = AnnotationFileReader.read(TARGETS, Annotation.class, new BasicAnnotation.Factory());
        		Map<String,Sequence> chrs = FastaSequenceIO.getChrSequencesFromFasta(GENOME_FASTA.getAbsolutePath());
				guides = getAllGuides(itr, chrs);
				chrs = null;  // release memory, in case Java gets around to it
				writeGuides(guides, allGuidesFile);
        	} else {
        		guides = AnnotationFileReader.load(allGuidesFile, GuideRNA.class, new GuideRNA.Factory());
        	}
			
        	if (!SKIP_SCORING) {
        		filterGuides(guides);
        		scoreGuides(guides);
        		writeGuides(guides, filteredFile);
        	} else {
        		guides = AnnotationFileReader.load(filteredFile, GuideRNA.class, new GuideRNA.Factory());
        	}
			
			AnnotationList<GuideRNAPair> guidePairs = getPairs(guides, MIN_DISTANCE, MAX_DISTANCE);
			writeGuidePairs(guidePairs, new File(OUTPUT_DIR.getPath() + "/allGuidePairs.bed"));
			
			guidePairs = filterGuidePairs(guidePairs);
			writeGuidePairs(guidePairs, filteredPairs);
			
			
        } catch (Exception e) {
        	log.error(e);
        }

        return 0;
    }
    
    
    public AnnotationList<GuideRNA> getAllGuides(CloseableIterator<Annotation> itr, Map<String,Sequence> chrs) {
    	AnnotationList<GuideRNA> allGuides = new AnnotationList<GuideRNA>();
    	while (itr.hasNext()) {
    		Gene target = new Gene(itr.next());
    		Collection<GuideRNA> guides = GuideRNA.findAll(chrs.get(target.getChr()), target.getStart(), target.getEnd(), target);
    		allGuides.addAll(guides);
    	}
    	return allGuides;
    }
    
    
    private void writeGuides(AnnotationList<GuideRNA> guides, File output) throws IOException {
    	BufferedWriter writer = new BufferedWriter(new FileWriter(output));
  		for (GuideRNA guide : guides) {
			writer.append(guide.toBedWithSequence() + "\n");
		}
  		writer.close();
    }
    
    
    private void filterGuides(AnnotationList<GuideRNA> guides) {
    	List<GuideRNA> guideList = guides.toList();
		CollectionUtils.filter(guideList, new UContentPredicate());
		CollectionUtils.filter(guideList, new PolyBaseFilter("ACGTN",5,5));
		guides = new AnnotationList<GuideRNA>();
		guides.addAll(guideList);
    }
    
    
    private void scoreGuides(AnnotationList<GuideRNA> guides) throws IOException, InterruptedException {
		GuideEfficacyScore efficacyScorer = new GuideEfficacyScore(guides.toList());
		GuideOffTargetScore offTargetScorer = null;
		if (OFF_TARGETS != null) offTargetScorer = new GuideOffTargetScore(OFF_TARGETS);
		
		for (GuideRNA guide : guides) {
			double total = 0.0;
			if (OFF_TARGETS != null)  {	
				// off target score is 0-100 where >50 is best
				double offTargetScore = offTargetScorer.getScore(guide);
				total += Math.floor(offTargetScore);
			}
			
			// efficacy score is 0-1, 1 is best
			double efficacyScore = efficacyScorer.getScore(guide);
			total += efficacyScore;
			
			guide.setScore(total);
		}
    }
    
    private AnnotationList<GuideRNAPair> filterGuidePairs(AnnotationList<GuideRNAPair> guides) {
    	AnnotationList<GuideRNAPair> guidePairs = new AnnotationList<GuideRNAPair>();  // necessary b/c AnnotationList does not support concurrent modification
    	GuideFilter filter = new GuideFilter();
    	for (GuideRNAPair pair : guides) {
    		if (filter.evaluate(pair.getLeft()) &&filter.evaluate(pair.getRight())) {
    			guidePairs.add(pair);
    		}
    	}
    	return(guidePairs);
    }
    
    
    private void writeGuidePairs(AnnotationList<GuideRNAPair> guides, File output) throws IOException {
    	BufferedWriter writer = new BufferedWriter(new FileWriter(output));
    	for (GuideRNAPair pair : guides) {
    		writer.append(pair.toBedWithSequence() + "\n");
    	}
    	writer.close();
    }
    
    
    public AnnotationList<GuideRNAPair> getPairs(AnnotationList<GuideRNA> guides, int minSpan, int maxSpan) {
    	AnnotationList<GuideRNAPair> pairs = new AnnotationList<GuideRNAPair>();
    	for (GuideRNA guide : guides) {
    		CloseableIterator<GuideRNA> matchItr = guides.getOverlappingAnnotations(new BasicAnnotation(guide.getReferenceName(), guide.getEnd() + minSpan, guide.getEnd() + maxSpan));
    		while (matchItr.hasNext()) {
    			pairs.add(new GuideRNAPair(guide, matchItr.next()));
    		}
    	}
    	return pairs;
    }
    
    
    
    private class GuideRNAPair extends BasicAnnotation {
    	private GuideRNA left, right;
    	
    	public GuideRNAPair(GuideRNA left, GuideRNA right) {
    		super(left);
    		addBlocks(right);
    		setName(left.getName() + "_" + right.getName());
    		this.left = left;
    		this.right = right;
    	}

		public GuideRNA getLeft() {
			return left;
		}

		public GuideRNA getRight() {
			return right;
		}
		
		public String toBedWithSequence() {
			return toBED() + "\t" + left.getSequenceWithPAM().getSequenceBases() + "\t" + right.getSequenceWithPAM().getSequenceBases() + "\t" + left.getScore() + "\t" + right.getScore();
		}
    }
    
    
    private class GuideFilter implements GuideRNAPredicate {

		@Override
		public boolean evaluate(GuideRNA guide) {
			double offTargetScore = Math.floor(guide.getScore());
			double efficacyScore = guide.getScore() - offTargetScore;
			return (offTargetScore > 30) && (efficacyScore > 0.1) && new UContentPredicate().evaluate(guide);
		}

		@Override
		public String getPredicateName() {
			return "scoreFilter";
		}

		@Override
		public String getShortFailureMessage(GuideRNA g) throws IOException,
				InterruptedException {
			return "did not pass score filter";
		}
    	
    }

}
