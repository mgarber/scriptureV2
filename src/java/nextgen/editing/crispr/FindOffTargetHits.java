package nextgen.editing.crispr;

import java.io.File;
import java.io.BufferedWriter;
import java.io.FileWriter;

import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.io.IoUtil;
import net.sf.picard.util.Log;
import net.sf.samtools.util.CloseableIterator;
import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.AnnotationFileReader;
import nextgen.core.annotation.AnnotationList;
import nextgen.core.annotation.BasicAnnotation;
import nextgen.editing.crispr.score.*;

import nextgen.editing.crispr.score.GuideOffTargetScore.OffTargetHits;

/**
 * Find all off target hits for a given guide RNA
 * @author engreitz
 */
public class FindOffTargetHits extends CommandLineProgram {
    private static final Log log = Log.getInstance(FindOffTargetHits.class);
    
	@Option(doc="Bed file containing guides")
	public File GUIDES;

	@Option(doc="Bitpacked file containing off target sites")
	public File OFF_TARGET_BITS;
	
	@Option(doc="BED file containing off target sites (23-mer sequence in name column)", optional=true)
	public File OFF_TARGETS = null;
	
	@Option(doc="Output file")
	public File OUTPUT;
	
	/**
	 * Stock main method.
	 *
	 * @param args main arguments
	 */
	public static void main(final String[] args) {
		System.exit(new FindOffTargetHits().instanceMain(args));
	}


    @Override
    protected int doWork() {
        IoUtil.assertFileIsReadable(GUIDES);
        IoUtil.assertFileIsReadable(OFF_TARGET_BITS);
        IoUtil.assertFileIsReadable(OFF_TARGETS);
        IoUtil.assertFileIsWritable(OUTPUT);

        try {
        	
        	AnnotationList<GuideRNA> guides = AnnotationFileReader.load(GUIDES, GuideRNA.class, new GuideRNA.Factory());
        	GuideOffTargetScore scorer = new GuideOffTargetScore(OFF_TARGET_BITS, OFF_TARGETS);
        	BufferedWriter writer = new BufferedWriter(new FileWriter(OUTPUT));
        	
			for (GuideRNA guide : guides) {
				OffTargetHits hits = scorer.getOffTargetHits(guide.getSequenceString());
				//log.info("Found " + hits.offTargetIndices.size() + " off target hits.  Score = " + hits.score + ". Annotating ...");
				scorer.annotateOffTargetHits(hits);
				//log.info("length of offTargetHits is " + hits.offTargetHits.size());
				
				for (GuideRNA offTargetHit : hits.offTargetHits) {
					writer.write(guide.getSequenceWithPAM().getSequenceBases() + "\t" + offTargetHit.toBedWithSequence() +"\n");
				}
			}
			writer.close();
			
        } catch (Exception e) {
        	log.error(e);
        }

        return 0;
    }
}
