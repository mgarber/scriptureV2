package broad.pda.methylation;

import java.io.File;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.Iterator;
import java.util.Map;
import java.util.LinkedHashMap;

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

import broad.core.math.EmpiricalDistribution;

public class FilterSamByBaseComposition extends CommandLineProgram {
    private static final Log log = Log.getInstance(FilterSamByBaseComposition.class);
	
    @Usage
    public String USAGE = "Filters a SAM/BAM file by read base composition and/or outputs metrics.";
    
    @Option(doc = "SAM or BAM file", optional=false, shortName = "I")
	public File INPUT;
    
    @Option(doc = "Output SAM or BAM file", optional=true, shortName = "O")
    public File OUTPUT = null;
    
    @Option(doc = "Filtered reads", optional=true, shortName = "FR")
    public File FILTERED_READS = null;
    
    @Option(doc = "Metrics file", optional=false, shortName="M")
    public File METRICS = null;
    
    @Option(doc = "Whether to filter", optional=true)
    public boolean FILTER = true;
    
    @Option(doc = "Reads with GC content above this will be filtered", optional=true)
    public double GC_CUTOFF = 0.5;
    
    
    @Override
    protected int doWork() {
    	try {
			IoUtil.assertFileIsReadable(INPUT);
			
			Map<String, EmpiricalDistribution> ed = getBaseCompositionMetrics(INPUT, OUTPUT, FILTERED_READS, FILTER);
			if (METRICS != null) {
				log.info("Writing metrics to " + METRICS.getPath());
				for (String s : ed.keySet())
					ed.get(s).write(METRICS.getAbsolutePath() + "." + s);
			}
    		
    	} catch (Exception e) {
    		log.error(e);
    		return 1;
    	}
    	return 0;
    }
    
    
    public Map<String,EmpiricalDistribution> getBaseCompositionMetrics(final File inSamOrBam, final File outSamOrBam, final File filteredReadsSamOrBam, boolean filter) {
    	SAMFileReader reader = new SAMFileReader(inSamOrBam);
    	final SAMFileHeader header = reader.getFileHeader();
    	SAMFileWriter writer = null, filteredReadsWriter = null;
    	if (outSamOrBam != null) writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(header, true, outSamOrBam);
    	if (filteredReadsSamOrBam != null) filteredReadsWriter = new SAMFileWriterFactory().makeSAMOrBAMWriter(header, true, filteredReadsSamOrBam);
   
    	Map<String,EmpiricalDistribution> ed = new LinkedHashMap<String,EmpiricalDistribution>();
    	EmpiricalDistribution pct = new EmpiricalDistribution(100, 0.0, 1.0);
    	EmpiricalDistribution ratio = new EmpiricalDistribution(100, -3.0, 3, true);
    	ed.put("gcPct", pct);
    	ed.put("gcRatio", ratio);
    	
		Iterator<SAMRecord> itr = reader.iterator();
		int countFiltered = 0;
		while (itr.hasNext()) {
			SAMRecord rec = itr.next();
			double gcPct = MethylationUtils.getGcPct(rec.getReadString());
			double gcRatio = MethylationUtils.getGcRatio(rec.getReadString());
			ed.get("gcPct").add(gcPct);
			ed.get("gcRatio").add(Math.log10(gcRatio));
			
			if (gcPct > GC_CUTOFF) {
				countFiltered++;
				if (filter && filteredReadsWriter != null) filteredReadsWriter.addAlignment(rec);
			} else {
				if (writer != null) writer.addAlignment(rec);
			}
		}
		
		if (filter) {
			log.info(countFiltered + " reads were filtered.");
		} else {
			log.info(countFiltered + " reads would have been filtered.");
		}
		
		if (writer != null) writer.close();
		if (filteredReadsWriter != null) filteredReadsWriter.close();
		return ed;
    }
    
    
	@Override
	protected String[] customCommandLineValidation() {
		if ((OUTPUT != null && INPUT.equals(OUTPUT)) || 
		    (FILTERED_READS != null && INPUT.equals(FILTERED_READS)) ||
		    (OUTPUT != null && FILTERED_READS != null && OUTPUT.equals(FILTERED_READS))) {
			return new String[]{"INPUT file, OUTPUT file, and FILTERED_READS file must differ!"};
		}
		if (FILTERED_READS != null && !FILTER) {
			return new String[]{"Cannot output FILTERED_READS if FILTER is false."};
		}
		return super.customCommandLineValidation();
	}
	
	/**
	 * Stock main method.
	 *
	 * @param args main arguments
	 */
	public static void main(final String[] args) {
		System.exit(new FilterSamByBaseComposition().instanceMain(args));
	}
}
