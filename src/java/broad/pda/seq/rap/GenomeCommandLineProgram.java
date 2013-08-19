package broad.pda.seq.rap;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.io.File;
import java.util.List;

import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.io.IoUtil;
import net.sf.picard.util.Log;

import nextgen.core.alignment.AbstractPairedEndAlignment.TranscriptionRead;
import nextgen.core.coordinatesystem.GenomicSpace;
import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.AnnotationList;
import nextgen.core.annotation.BasicAnnotation;
import nextgen.core.model.AlignmentModel;
import nextgen.core.readFilters.*;


public abstract class GenomeCommandLineProgram extends CommandLineProgram {
    private static final Log log = Log.getInstance(GenomeCommandLineProgram.class);
    
	@Option(doc="File specifying chromosome sizes.")
	public String SIZES = "/seq/lincRNA/data/mm9/sizes";
	
	@Option(doc="File containing masked regions.", optional=true)
	public String MASK_FILE = "/seq/mguttman/ChIPData/MaskFiles/MM9Segments/all.mask.mouse.n36.d2.bin";
	
	@Option(doc="Percent masked allowable per sliding window", optional=true)
	public double PCT_MASKED_ALLOWED = 50.0;
	
	@Option(doc="Region to process (e.g. chr1, chr1:5000-50230).  Default (null) processes the entire genome", optional=true)
	public String REGION = null;
	
	@Option(doc="Maximum paired-end read fragment length to consider.", optional=true)
	public Integer MAX_FRAGMENT_LENGTH = 10000;
	
	@Option(doc="Minimum mapping quality for reads") 
	public Integer MIN_MAPPING_QUALITY = 30;
	
	protected GenomicSpace coordinateSpace;

	
	
	@Override
	protected String[] customCommandLineValidation() {
		loadCoordinateSpace();
		return super.customCommandLineValidation();
	}
	
	protected void loadCoordinateSpace() {
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
	
	public Map<String,Annotation> getRegionsMap() {
		Map<String,Annotation> regions = new HashMap<String,Annotation>();

		if (REGION != null) {
			if (coordinateSpace.hasChromosome(REGION)) {
				regions.put(REGION,coordinateSpace.getReferenceAnnotation(REGION));
			} else {
				try {
					regions.put(REGION,new BasicAnnotation(REGION));
				} catch (RuntimeException e) {
					throw new IllegalArgumentException("REGION is improperly formatted");
				}
			}
		} else {
			for(String chr:coordinateSpace.getReferenceNames()){
				regions.put(chr,coordinateSpace.getReferenceAnnotation(chr));
			}
		}
		return regions;
	}
	
	public AnnotationList<Annotation> getRegionSet() {
		AnnotationList<Annotation> annotations = new AnnotationList<Annotation>(coordinateSpace, getRegions());
		return annotations;
	}


	protected GenomicSpace getCoordinateSpace() { return coordinateSpace; }
	
	
	public AlignmentModel loadAlignmentModel(File bamFile) {
		return loadAlignmentModel(bamFile, true);
	}

	public AlignmentModel loadAlignmentModel(File bamFile, boolean pairedEnd) {
		// TODO: Extend to handle BED or other annotation files
		
		IoUtil.assertFileIsReadable(bamFile);
		AlignmentModel model = new AlignmentModel(bamFile.getAbsolutePath(), coordinateSpace, pairedEnd);
		model.addFilter(new GenomicSpanFilter(MAX_FRAGMENT_LENGTH));
		model.addFilter(new ChimeraFilter());
		model.addFilter(new ProperPairFilter());
		model.addFilter(new MappingQualityFilter(MIN_MAPPING_QUALITY));
		// TODO need to modify PairedEndWriter to save information about the other read
		return model;
	}
}
