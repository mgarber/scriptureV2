package nextgen.core.programs;

import java.io.File;
import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.util.Log;

import nextgen.core.writers.PairedEndWriter;

public class ConvertBamToPairedEnd extends CommandLineProgram {
    private static final Log log = Log.getInstance(ConvertBamToPairedEnd.class);
    
    @Option(doc="Input BAM file (default format)", shortName="I")
    String INPUT;
    
    @Option(doc="Maximum insert length allowed")
    Integer MAX_INSERT = 5000000;
    
	/**
	 * Stock main method.
	 *
	 * @param args main arguments
	 */
	public static void main(final String[] args) {
		System.exit(new ConvertBamToPairedEnd().instanceMain(args));
	}
	

	@Override
	protected int doWork() {
		
		try {
			PairedEndWriter writer = new PairedEndWriter(new File(INPUT));
			writer.setMaxAllowableInsert(MAX_INSERT);
			writer.convertInputToPairedEnd();
			
		} catch (Exception e) {
			log.error(e);
		}
		
		return 0;
	}

	
	
}
