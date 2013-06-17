package broad.pda.seq.rap;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.Iterator;
import java.util.List;

import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.util.Log;

import nextgen.core.annotation.*;
import nextgen.core.model.score.*;


public class SlideAndCalculate extends GenomeScoringProgramFromBed {
    private static final Log log = Log.getInstance(SlideAndCalculate.class);
	
    @Usage
    public String USAGE = "Slides across the genome and counts reads (or calculates ratios for two models).  This is currently written for Genomic Space analyses but should be modified to allow other coordinate systems.";
   
	@Option(doc="Window size")
	public int WINDOW;
	
	@Option(doc="Overlap between windows")
	public int OVERLAP;

	@Option(doc="Output file", shortName="O")
	public File OUTPUT;

	@Option(doc="Whether to force paired end behavior", optional=true)
	public boolean PAIRED_END=false;
	

	/**
	 * Stock main method.
	 *
	 * @param args main arguments
	 */
	public static void main(final String[] args) {
		System.exit(new SlideAndCalculate().instanceMain(args));
	}

	@Override
	protected int doWork() {
		
		try {
			if (OUTPUT.exists()) {
				OUTPUT.delete();
			}
			
			List<Annotation> regions = getRegions();

			BufferedWriter bw = new BufferedWriter(new FileWriter(OUTPUT,true));
			WindowProcessor<? extends WindowScore> processor = getWindowProcessor();
					
			for (Annotation region : regions) {
				/*		log.info("Starting: " + region.toUCSC());
				Iterator<? extends Annotation> windowIterator = getCoordinateSpace().getWindowIterator(region, WINDOW, OVERLAP);
				WindowScoreIterator<? extends WindowScore> itr = new WindowScoreIterator(windowIterator, processor, region);
				WindowScore curr = processor.processWindow(region);
				bw.write(curr.toString());
				bw.newLine();
				while (itr.hasNext()) {
				
					WindowScore curr = itr.next();
					
				}
				itr.close();*/
				WindowScore curr = processor.processWindow(region);
				bw.write(curr.toString());
				bw.newLine();
			}
			
			bw.close();
			
		} catch (Exception e) {
			log.error(e);
		}
		
		return 0;
	}
	
}
