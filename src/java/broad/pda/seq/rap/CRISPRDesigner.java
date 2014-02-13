package broad.pda.seq.rap;

import java.util.*;
import java.util.regex.Pattern;
import java.io.*;


import org.apache.commons.collections15.Predicate;
import org.apache.commons.collections15.CollectionUtils;

import com.mysql.jdbc.StringUtils;

import broad.core.motif.SequenceMotif;
import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;
import broad.core.sequence.SequenceRegion;
import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.io.IoUtil;
import net.sf.picard.util.Log;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
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
	
	@Option(doc="BED file containing off target sites (23-mer sequence in name column)")
	public File OFF_TARGETS;

	@Option(doc="Output directory to create")
	public File OUTPUT_DIR;
	
	@Option(doc="Skip initial step of generating guides ... assumes that file exists and reads from there")
	public Boolean SKIP_GENERATION = false;

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
        IoUtil.assertFileIsReadable(OFF_TARGETS);

        if (!OUTPUT_DIR.exists()) OUTPUT_DIR.mkdir();
        
        File allGuidesFile = new File(OUTPUT_DIR.getAbsolutePath() + "/allGuides.bed");

        
        try {
        	
        	List<GuideRNA> guides = null;
        	if (!SKIP_GENERATION) {
        		// Iterate through buffered query annotations to handle large files gracefully
        		CloseableIterator<Annotation> itr = AnnotationFileReader.read(TARGETS, Annotation.class, new BasicAnnotation.Factory());
        		Map<String,Sequence> chrs = FastaSequenceIO.getChrSequencesFromFasta(GENOME_FASTA.getAbsolutePath());
				guides = getAllGuides(itr, chrs);
				chrs = null;  // release memory, in case Java gets around to it
				writeGuides(guides, allGuidesFile);
        	} else {
        		guides = AnnotationFileReader.load(allGuidesFile, GuideRNA.class, new GuideRNA.Factory()).toList();
        	}
			
			CollectionUtils.filter(guides, new GuideRNASeedU());
			CollectionUtils.filter(guides, new PolyBaseFilter("ACGTN",5,5));
			
			GuideEfficacyScore score = new GuideEfficacyScore(guides);
			GuideOffTargetScore offTargetScore = new GuideOffTargetScore(OFF_TARGETS);
			
			for (GuideRNA guide : guides) {
				// off target score is 0-100 where >50 is best
				// efficacy score is 0-1, 1 is best
				guide.setScore(Math.floor(offTargetScore.getScore(guide)) + score.getScore(guide));
			}
			
			writeGuides(guides, new File(OUTPUT_DIR.getAbsolutePath() + "/filteredGuides.bed"));

			
			
        } catch (Exception e) {
        	log.error(e);
        }

        return 0;
    }
    
    
    public List<GuideRNA> getAllGuides(CloseableIterator<Annotation> itr, Map<String,Sequence> chrs) {
    	List<GuideRNA> allGuides = new ArrayList<GuideRNA>();
    	while (itr.hasNext()) {
    		Gene target = new Gene(itr.next());
    		Collection<GuideRNA> guides = GuideRNA.findAll(chrs.get(target.getChr()), target.getStart(), target.getEnd(), target);
    		allGuides.addAll(guides);
    	}
    	return allGuides;
    }
    
    
    public void writeGuides(List<GuideRNA> guides, File output) throws IOException {
    	BufferedWriter writer = new BufferedWriter(new FileWriter(output));
  		for (GuideRNA guide : guides) {
			writer.append(guide.toBedWithSequence() + "\n");
		}
  		writer.close();
    }
    

}
