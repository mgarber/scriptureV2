package nextgen.core.programs;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

import org.apache.log4j.Logger;


import broad.core.annotation.ShortBED;


import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.samtools.TabixWriter;
import net.sf.samtools.TabixWriter.TabixException;
import net.sf.samtools.util.BlockCompressedInputStream;
import net.sf.samtools.util.BlockCompressedOutputStream;
import net.sf.samtools.util.SortingCollection;
import nextgen.core.alignment.Alignment;
import nextgen.core.readers.BamToPairedEndIterator;
import nextgen.core.utils.ShortBEDCodec;
import nextgen.core.utils.ShortBEDComparator;

/**
 * 
 * TO DO:
 * 
 * 3. move to nextgen.core.programs
 * 
 * @author zhuxp
 *
 */
public class BamToPairedEndSortedShortBED extends CommandLineProgram {
	
	private static final String PROGRAM_VERSION = "0.01";
	private static final int MAX_NUM = 500000;
	private static final String SUFFIX=".PairedEnd.ShortBED.gz";
	@Usage
    public String USAGE =  "java nextgen.core.programs.BamToPairedEndSortedShortBED I=file.bam";
	static Logger logger = Logger.getLogger(BamToPairedEndSortedShortBED.class.getName());
    @Option(doc="INPUT STANDARD BAM FILE", shortName=StandardOptionDefinitions.INPUT_SHORT_NAME) public File INPUT;
    public static void main(String[] argv){
    	System.exit(new BamToPairedEndSortedShortBED().instanceMain(argv));
        
    }
    
    protected int doWork()
    {
    	BamToPairedEndIterator iter = new BamToPairedEndIterator(INPUT);
    	String OUTPUT=INPUT.toString()+SUFFIX;
    	if (Files.exists(Paths.get(OUTPUT)))
    			{
    		      logger.warn("overwriting exists file "+OUTPUT);
    		      try {
					Files.deleteIfExists(Paths.get(OUTPUT));
				      } catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				      }
    			}
    	BlockCompressedOutputStream out = new BlockCompressedOutputStream(OUTPUT);
    	
    	File tmpDir = new File(System.getProperty("java.io.tmpdir"));
		//File tmpDir = new File("/Users/zhuxp/tmp");
    	logger.info("BEGINNING");
    	SortingCollection<ShortBED> sortingArray = SortingCollection.newInstance(ShortBED.class, new ShortBEDCodec(), new ShortBEDComparator(), MAX_NUM , tmpDir);
        int i=0;
		while (iter.hasNext()) {
			i++;
			if (i%100000==0) 
				logger.info("reading "+i );
		    Alignment  record= iter.next();
		    ShortBED r=new ShortBED(record.getName(),record.getChr(),record.getStart(),record.getEnd());
		    sortingArray.add(r);
		}	
		iter.close();
		logger.info("done reading");
		sortingArray.doneAdding();
		logger.info("done adding");
		logger.info(tmpDir.getAbsoluteFile());
		i=0;
		for (ShortBED bed : sortingArray)
		{
			i++;
			if (i%100000==0) 
				logger.info("writing "+i );
		    	//System.out.println(bed.toShortBED());
			try {
				String line=bed.toShortBED()+"\n";
				out.write(line.getBytes());
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
			
		}
		try {
			out.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		/*
		 *  Create Tabix Index
		 */
		
		Path fp = Paths.get(OUTPUT);
		
		TabixWriter tabixWriter;
		try {
			tabixWriter = new TabixWriter(fp,TabixWriter.BED_CONF);
			try {
				Files.deleteIfExists(Paths.get(OUTPUT+".tbi"));
			      } catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			      }
			tabixWriter.createIndex();
		} catch (IOException | TabixException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		
		
		return 0;
    }
}

