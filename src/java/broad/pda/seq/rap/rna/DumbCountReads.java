package broad.pda.seq.rap.rna;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.*;

import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.util.Log;
import net.sf.samtools.*;
import nextgen.core.alignment.Alignment;
import nextgen.core.alignment.SingleEndAlignment;
import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.AnnotationFileReader;
import nextgen.core.annotation.AnnotationList;
import nextgen.core.annotation.BasicAnnotation;
import nextgen.core.annotation.filter.FullyContainedFilter;
import nextgen.core.model.score.WindowProcessor;
import nextgen.core.model.score.WindowScore;
import nextgen.core.model.score.WindowScoreIterator;


public class DumbCountReads extends CommandLineProgram {
    private static final Log log = Log.getInstance(DumbCountReads.class);

    @Usage
    public String USAGE = "Counts reads overlapping the provided BED files.";
    
    @Option(doc="Input SAM or BAM files", shortName="I")
    public List<File> INPUT;
    
    @Option(doc="Target BED file")
    public File ANNOTATION_FILE;
    
    @Option(doc="Buffer")
    public Integer BUFFER = 1000;
    
    @Option(doc="Output BED file", shortName="O")
    public File OUTPUT;
    
	/**
	 * Stock main method.
	 *
	 * @param args main arguments
	 */
	public static void main(final String[] args) {
		System.exit(new DumbCountReads().instanceMain(args));
	}
	

	@Override
	protected int doWork() {
		
		try {
			if (OUTPUT.exists()) {
				OUTPUT.delete();
			}
        	
			
			List<SAMFileReader> readers = new ArrayList<SAMFileReader>();
			for (File input : INPUT) {
				readers.add(new SAMFileReader(input));
			}

			AnnotationList<Annotation> targets = AnnotationFileReader.load(ANNOTATION_FILE, Annotation.class, new BasicAnnotation.Factory());
			log.info("Loaded " + targets.size() + " annotations.");
			
			BufferedWriter bw = new BufferedWriter(new FileWriter(OUTPUT,true));
			bw.write("Target");
			for (File input : INPUT) {
				bw.write("\t");
				bw.write(input.getPath());
			}
			bw.write("\n");
			
			for (Annotation region : targets) {
				log.info("Starting: " + region.toUCSC());
				
				bw.write(region.getName());
				
				for (SAMFileReader reader : readers) {
					int count = 0;
					int antisense = 0;
					SAMRecordIterator itr = reader.query(region.getReferenceName(), region.getStart()-BUFFER, region.getEnd()+BUFFER, false);
					while (itr.hasNext()) {
						SAMRecord rec = itr.next();
						Alignment read = new SingleEndAlignment(rec);
				    	if (region.getStrand() == read.getStrand() && !rec.getFirstOfPairFlag() ||
				    		region.getStrand() != read.getStrand() && rec.getFirstOfPairFlag()) {
				    			++count;
				    	} else {
				    		++antisense;
				    	}
					}
					itr.close();
					bw.write("\t" + count + ";" + antisense);
				}
				
				bw.write("\n");
			}
			
			bw.close();
			
		} catch (Exception e) {
			log.error(e);
		}
		
		return 0;
	}
}
