package broad.pda.seq.rap;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;

import broad.core.annotation.ShortBEDReader;

import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.io.IoUtil;
import net.sf.picard.util.Log;
import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.AnnotationFileReader;
import nextgen.core.annotation.AnnotationList;
import nextgen.core.annotation.BasicAnnotation;
import nextgen.core.annotation.filter.FullyContainedFilter;
import nextgen.core.model.score.WindowProcessor;
import nextgen.core.model.score.WindowScore;
import nextgen.core.model.score.WindowScoreIterator;

public class PlotAggregateRegions extends GenomeScoringProgram {
    private static final Log log = Log.getInstance(PlotAggregateRegions.class);
    
    @Usage
    public String USAGE = "Plot scores as as a function of position from / in a target.";
    
    @Option(doc="Output file", shortName="O")
    public File OUTPUT;
    
	@Option(doc="File containing genomic annotations (bed or bedGraph format)")
	public File ANNOTATION_FILE;
	
	@Option(doc="Window size")
	public int WINDOW;

	@Option(doc="Overlap between windows", optional=true)
	public int OVERLAP=0;
	
	@Option(doc="Length of the inner region", optional=true)
	public int INNER_LENGTH = 0;
	
	@Option(doc="Length of the outer region", optional=true)
	public int OUTER_LENGTH = 0;
	
	@Option(doc="Length of the middle region", optional=true)
	public int MIDDLE_LENGTH = 0;
	
	@Option(doc="Set to true to generate one-sided diagrams")
	public boolean SYMMETRIC=false;
	

	/**
	 * Stock main method.
	 *
	 * @param args main arguments
	 */
	public static void main(final String[] args) {
		System.exit(new PlotAggregateRegions().instanceMain(args));
	}
	
	@Override
	protected int doWork() {
		try {

			List<Annotation> regions = getRegions();
			AnnotationList<Annotation> targets = AnnotationFileReader.load(ANNOTATION_FILE, Annotation.class, new BasicAnnotation.Factory(), getCoordinateSpace(), new FullyContainedFilter(regions));
			log.info("Loaded " + targets.size() + " annotations.");
			
			BufferedWriter bw = new BufferedWriter(new FileWriter(OUTPUT,true));
			WindowProcessor<? extends WindowScore> processor = getWindowProcessor();
			
			for (Annotation region : regions) {
				log.info("Starting: " + region.toUCSC());
				
				PlotRegions subregions = generateSubregions(region);
				
				
			}
			
			bw.close();
			
			
		} catch (Exception e) {
			log.error(e);
		}
		
		return 0;
	}
	
	public PlotRegions generateSubregions(Annotation region) throws IOException {
		PlotRegions subregions;
		
		Annotation outerLeft = new BasicAnnotation(region.getReferenceName(), region.getStart() - OUTER_LENGTH, region.getStart(), region.getOrientation());
		Annotation outerRight = new BasicAnnotation(region.getReferenceName(), region.getEnd(), region.getEnd() + OUTER_LENGTH, region.getOrientation());
		Annotation innerLeft = new BasicAnnotation(region.getReferenceName(), region.getStart(), region.getStart() + INNER_LENGTH, region.getOrientation());
		Annotation innerRight = new BasicAnnotation(region.getReferenceName(), region.getEnd() - INNER_LENGTH, region.getEnd(), region.getOrientation());
		// TODO:  Adjust so that these regions don't extend beyond the length of a short gene
		
		int midway = (int) Math.floor((region.getEnd() + region.getStart())/2.0);
		int midAdjust = (int) Math.floor(MIDDLE_LENGTH/2.0);
		Annotation middle = new BasicAnnotation(region.getReferenceName(), midway - midAdjust, midway - midAdjust + MIDDLE_LENGTH, region.getOrientation());
		
		if (region.isNegativeStrand()) {
			subregions = new PlotRegions(outerRight, innerRight, middle, innerLeft, outerLeft);
		} else {
			subregions = new PlotRegions(outerLeft, innerLeft, middle, innerRight, outerRight);
		} 
		
		if (SYMMETRIC) {
			reverseOrientation(subregions.endOuter);
			reverseOrientation(subregions.endInner);
		}
		
		return subregions;
	}
	
	
	private void reverseOrientation(Annotation a) {
		if (a.isNegativeStrand()) {
			a.setOrientation(Annotation.Strand.POSITIVE);
		} else {
			a.setOrientation(Annotation.Strand.NEGATIVE);
		}
	}
	
	
	private void scanAndPrint(Annotation subregion, String name, BufferedWriter bw) {
		
	}
	
	
	private class PlotRegions {
		public Annotation beginOuter, beginInner, middle, endInner, endOuter;
		public PlotRegions(Annotation beginOuter, Annotation beginInner, Annotation middle, Annotation endInner, Annotation endOuter) {
			this.beginOuter = beginOuter;
			this.beginInner = beginInner;
			this.middle = middle;
			this.endInner = endInner;
			this.endOuter = endOuter;
		}
	}
}
