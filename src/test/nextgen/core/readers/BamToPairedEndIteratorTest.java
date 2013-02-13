package nextgen.core.readers;

import java.io.File;
import java.io.IOException;

import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.samtools.util.CloseableIterator;
import nextgen.core.alignment.Alignment;
import nextgen.core.alignment.PairedEndAlignment;
import nextgen.core.readers.BamToPairedEndIteratorTest;

public class BamToPairedEndIteratorTest extends CommandLineProgram {
	
	private static final String PROGRAM_VERSION = "0.01";
	@Usage
    public String USAGE =  "java nextgen.core.readers.BamToPairedEndIteratorTest I=file.bam";

    @Option(doc="This is test paired end bam iterater and convert to simple bed format.", shortName=StandardOptionDefinitions.INPUT_SHORT_NAME) public File INPUT;
	
    public static void main(String[] argv){
    	System.exit(new BamToPairedEndIteratorTest().instanceMain(argv));
        
    }
    
    protected int doWork()
    {
    	BamToPairedEndIterator iter = new BamToPairedEndIterator(INPUT);
		//SAMFileHeader header = reader.getFileHeader();
    	System.out.println("BEGINNING");
		
        int i=0;
		while (iter.hasNext()) {
		    Alignment  record= iter.next();
			//System.out.println("aha");
		    System.out.println(record);
		    i+=1;
		  if ( record instanceof PairedEndAlignment)
		  {
			    System.out.println("First Mate:\t"+((PairedEndAlignment) record).getFirstMate());
			    System.out.println("Second Mate:\t"+((PairedEndAlignment) record).getSecondMate());
		  }
			    System.out.println();
		   if (i%10000==0) {System.err.format("Processing %d reads      BufferSize:%d         \r",i,iter.getBufferSize());}
		}	
		iter.close();
		
		System.out.print("ENDING");
		return 0;
    }
}

