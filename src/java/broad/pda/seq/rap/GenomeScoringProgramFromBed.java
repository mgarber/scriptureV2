package broad.pda.seq.rap;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import nextgen.core.model.score.*;
import nextgen.core.annotation.*;
import nextgen.core.coordinatesystem.GenomicSpace;

public abstract class GenomeScoringProgramFromBed extends CommandLineProgram{
	
	@Option(doc="Input bed file", shortName="I")
	public File TARGET;
	
	@Option(doc="Control bed file", shortName="C")
	public File CONTROL = null;
	
	@Option(doc="Scoring function to use")
	public String SCORE = "all";

	@Option(doc="File specifying chromosome sizes.")
	public String SIZES = "/seq/lincRNA/data/mm9/sizes";
	
	@Option(doc="File containing masked regions.", optional=true)
	public String MASK_FILE = "/seq/mguttman/ChIPData/MaskFiles/MM9Segments/all.mask.mouse.n36.d2.bin";
	
	@Option(doc="Percent masked allowable per sliding window", optional=true)
	public double PCT_MASKED_ALLOWED = 50.0;
	
	@Option(doc="Region to process (e.g. chr1, chr1:5000-50230).  Default (null) processes the entire genome", optional=true)
	public String REGION = null;

	protected GenomicSpace coordinateSpace;
	
	protected String[] customCommandLineValidation() {
		loadCoordinateSpace();
		return super.customCommandLineValidation();
	}
	
	public WindowProcessor<? extends WindowScore> getWindowProcessor() {
		AnnotationCollection<? extends Annotation> target = AnnotationFileReader.load
				(TARGET, broad.core.annotation.BEDGraph.class, new broad.core.annotation.AnnotationFactoryFactory.BEDGraphFactory(), getCoordinateSpace());
		AnnotationCollection<? extends Annotation> control = AnnotationFileReader.load
				(CONTROL, broad.core.annotation.BEDGraph.class, new broad.core.annotation.AnnotationFactoryFactory.BEDGraphFactory(), getCoordinateSpace());
		return getWindowProcessor(target,control);
	}

	//TODO: Do we need this?
	private void loadCoordinateSpace() {
		coordinateSpace = new GenomicSpace(SIZES, MASK_FILE, PCT_MASKED_ALLOWED);
	}
	
	public List<Annotation> getRegions() {
		List<Annotation> regions = new ArrayList<Annotation>();

		if (REGION != null) {
			if (coordinateSpace.hasChromosome(REGION)) {
				regions.add(coordinateSpace.getReferenceAnnotation(REGION));
			} else {
				try {
					regions.add(new BasicAnnotation(REGION));
				} catch (RuntimeException e) {
					throw new IllegalArgumentException("REGION is improperly formatted");
				}
			}
		} else {
			regions.addAll(coordinateSpace.getReferenceAnnotations());
		}
		return regions;
	}
	
	protected GenomicSpace getCoordinateSpace() { return coordinateSpace; }
	
	public WindowProcessor<? extends WindowScore> getWindowProcessor(AnnotationCollection<? extends Annotation> target,AnnotationCollection<? extends Annotation> control) {
		WindowProcessor<? extends WindowScore> processor;
		
		if (SCORE.equalsIgnoreCase("ratio")) {
			processor = new NewRatioScore.Processor(target, control);
		} else{
			throw new IllegalArgumentException("Could not find scoring class for " + SCORE);
		}
		
		return processor;
	}
}
