package broad.pda.methylation;

import java.io.File;
import java.util.Iterator;

import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.io.IoUtil;
import net.sf.picard.util.Log;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;


import org.apache.commons.lang3.StringUtils;

public class FilterSamForMethylationAlignment extends CommandLineProgram {
    private static final Log log = Log.getInstance(FilterSamForMethylationAlignment.class);
	
    @Usage
    public String USAGE = "Filters a SAM/BAM file for methylation alignment.";
    
    @Option(doc = "SAM or BAM file", optional=false, shortName = "I")
	public File INPUT;
    
    @Option(doc = "Output SAM or BAM file", shortName = "O")
    public File OUTPUT = null;
    
    @Option(doc = "Filtered reads", shortName = "FR")
    public File FILTERED_READS = null;
   
    
    @Override
    protected int doWork() {
    	try {
			IoUtil.assertFileIsReadable(INPUT);
			IoUtil.assertFileIsWritable(OUTPUT);
			IoUtil.assertFileIsWritable(FILTERED_READS);
			
			SAMFileReader reader = new SAMFileReader(INPUT);
	    	final SAMFileHeader header = reader.getFileHeader();
	    	SAMFileWriter writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(header, true, OUTPUT);
	    	SAMFileWriter filterWriter = new SAMFileWriterFactory().makeSAMOrBAMWriter(header, true, FILTERED_READS);
	    	
	    	Iterator<SAMRecord> itr = reader.iterator();
			int countFilteredN = 0;
			int countFilteredGC = 0;
			int counter = 0;
			
			while (itr.hasNext()) {
				SAMRecord rec = itr.next();
				
				boolean filter = false;
				
				// Count number of N bases and filter those with 2+
				if (StringUtils.countMatches(rec.getReadString(), "N") > 1) {
					filter = true;
					countFilteredN++;
				}
				
				// Check to see if the read might be unconverted.
				if (MethylationUtils.guessUnconverted(rec.getReadString())) {
					filter = true;
					countFilteredGC++;
				}

				if (!filter) writer.addAlignment(rec);
				else filterWriter.addAlignment(rec);
				
				counter++;
				if (counter % 1000000 == 0) log.info(counter + " reads completed.");
			}

			log.info(countFilteredN + " reads were filtered due to too many N bases.");
			log.info(countFilteredGC + " reads were filtered due to appearing to be unconverted.");
			log.info((counter - countFilteredN - countFilteredGC) + " reads were retained.");
			
			writer.close();
			filterWriter.close();
    		
    	} catch (Exception e) {
    		log.error(e);
    		return 1;
    	}
    	return 0;
    }
    
    
	@Override
	protected String[] customCommandLineValidation() {
		if ((OUTPUT != null && INPUT.equals(OUTPUT))) {
			return new String[]{"INPUT file and OUTPUT file must differ!"};
		}
		return super.customCommandLineValidation();
	}
	
	/**
	 * Stock main method.
	 *
	 * @param args main arguments
	 */
	public static void main(final String[] args) {
		System.exit(new FilterSamForMethylationAlignment().instanceMain(args));
	}
}
