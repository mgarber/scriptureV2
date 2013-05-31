package nextgen.core.model.score;

import broad.pda.seq.segmentation.AlignmentDataModelStats;
import net.sf.samtools.util.CloseableIterator;
import nextgen.core.alignment.Alignment;
import nextgen.core.annotation.Annotation;
import nextgen.core.coordinatesystem.CoordinateSpace;
import nextgen.core.model.AlignmentModel;

public class ScanStatisticScore extends CountScore {

	/**
	 * Global statistics are stored in the CountScore class
	 *
	 * P-Value statistics are stored here
	 * "Region" refers to an arbitrary region being scanned, like a gene, or a chromosome, or some subset of the whole coordinate space
	 */
	private double scanPvalue;
	private double regionLength;
	private double fullyContainedNumberOfReads;
	private double globalLength;
	private CoordinateSpace coordSpace;

	
	public ScanStatisticScore(AlignmentModel model, Annotation annotation) {
		super(model, annotation);
		coordSpace = model.getCoordinateSpace();
		setGlobalLength(model.getGlobalLength());
		try {
			setScanPvalue(AlignmentDataModelStats.calculatePVal(new Double(getCount()).intValue(), model.getGlobalLambda(), model.getCoordinateSpace().getSize(annotation), getGlobalLength()));
		} catch(Exception e) {
			logger.info("Could not set scan P value for annotation " + annotation.getName());
		}
		getAnnotation().setScore(getScanPvalue());
		
		// by default, set the "region" stats to the "global" stats
		setRegionLength(model.getGlobalLength());
		setRegionTotal(getTotal());
		
	}
	
	public ScanStatisticScore(AlignmentModel model, Annotation annotation, ScanStatisticScore previousScore, double newScore) {
		super(previousScore, annotation, newScore); //Set the new score without computing
		coordSpace = model.getCoordinateSpace();
		setRegionLength(previousScore.getRegionLength());
		setRegionTotal(previousScore.getRegionTotal());
		setGlobalLength(model.getGlobalLength());
		try {
			setScanPvalue(AlignmentDataModelStats.calculatePVal(new Double(getCount()).intValue(), model.getGlobalLambda(), model.getCoordinateSpace().getSize(annotation), model.getGlobalLength()));
		} catch(Exception e) {
			logger.info("Could not set scan P value for annotation " + annotation.getName());
		}
		getAnnotation().setScore(getScanPvalue());
	}
	
	public ScanStatisticScore(AlignmentModel model, Annotation annotation, double regionTotal, double regionLength) {
		this(model, annotation);
		setRegionTotal(regionTotal);
		setRegionLength(regionLength);
	}
	
	/**
	 * Control region length and total, global length and total
	 * Scan P value depends on global stats
	 * @param model Alignment data
	 * @param annotation Annotation
	 * @param regionTotal Total number of fragments mapping to local region
	 * @param regionLength Total length of local region
	 * @param globalTotal Global number of fragments
	 * @param globalLength Global length
	 */
	public ScanStatisticScore(AlignmentModel model, Annotation annotation, double regionTotal, double regionLength, double globalTotal, double globalLength) {
		super(model, annotation);
		coordSpace = model.getCoordinateSpace();
		setRegionTotal(regionTotal);
		setTotal(globalTotal);
		setRegionLength(regionLength);
		setGlobalLength(globalLength);
		setScanPvalue(AlignmentDataModelStats.calculatePVal(new Double(getCount()).intValue(), getGlobalLambda(), model.getCoordinateSpace().getSize(annotation), getGlobalLength()));
	}
	
	public CoordinateSpace getCoordinateSpace() {
		return coordSpace;
	}
	
	public void refreshScanPvalue(AlignmentModel model) {
		setScanPvalue(AlignmentDataModelStats.calculatePVal(new Double(getCount()).intValue(), getTotal() / globalLength, model.getCoordinateSpace().getSize(annotation), globalLength));
	}
	
	public double getAverageCoverage(AlignmentModel data) { 
		int regionSize = coordSpace.getSize(annotation);
		CloseableIterator<Alignment> readsIter = data.getOverlappingReads(getAnnotation(), false);
		int basesInReads = 0;
		while(readsIter.hasNext()) {
			Alignment read = readsIter.next();
			basesInReads += read.getOverlap(annotation);
		}
		double avgCoverage = (double) basesInReads / (double)regionSize;
		return avgCoverage;
	}

	public double getEnrichment(AlignmentModel data) {
		return getAverageCoverage(data) / getLocalLambda();
	}

	public double getLocalLambda() {
		return getRegionTotal() / getRegionLength();
	}

	public void setGlobalLength(double d) {
		globalLength = d;
	}
	
	public double getGlobalLambda() {
		return getTotal() / getGlobalLength();
	}
	
	public double getGlobalLength() {
		return globalLength;
	}
	
	public void setRegionLength(double regionLength) {
		this.regionLength = regionLength;
	}

	public double getRegionLength() {
		return regionLength;
	}

	public void setScanPvalue(double scanPvalue) {
		this.scanPvalue = scanPvalue;
	}

	public double getScanPvalue() {
		return scanPvalue;
	}

	public void setFullyContainedNumberOfReads(double fullyContainedNumberOfReads) {
		this.fullyContainedNumberOfReads = fullyContainedNumberOfReads;
	}

	public double getFullyContainedNumberOfReads() {
		return fullyContainedNumberOfReads;
	}
	
	
	public String toString() {
		return super.toString() + "\t" + getScanPvalue() + "\t" 
		+ getGlobalLambda() + "\t" + getLocalLambda() + "\t" + getRegionLength();
	}
	
	
	public static class Processor extends WindowProcessor.AbstractProcessor<ScanStatisticScore> {
		protected AlignmentModel model;
		protected double regionTotal = DEFAULT_REGION_TOTAL;
		protected double regionLength = DEFAULT_REGION_TOTAL;
		
		public Processor(AlignmentModel model) {
			this.model = model;
		}
		
		public ScanStatisticScore processWindow(Annotation annotation) {
			return new ScanStatisticScore(model, annotation, regionTotal, regionLength);
		}
		
		public void initRegion(Annotation region) {
			if (region != null) {
				regionTotal = model.getCount(region);
				regionLength = region.length();
			}
		}

		/**
		 * Compute the count using the previous windowScore
		 * @param nextRegion 
		 * @param previousScore The WindowScore before
		 * @return the count of the current window
		 */
		private double computeCount(Annotation nextRegion, CountScore previousScore) {
			//else, get the minus region scores and the plus value scores
			//This is not so simple because we'll need to use the fully contained regions
			double subtractVal=model.getCountExcludingRegion(previousScore.getAnnotation().minus(nextRegion), nextRegion);
			double addVal=model.getCountExcludingRegion(nextRegion.minus(previousScore.getAnnotation()), previousScore.getAnnotation());
			return (previousScore.getCount()-subtractVal)+addVal;
		}
		
		@Override
		public ScanStatisticScore processWindow(Annotation annotation, ScanStatisticScore previousScore) {
			//if the previous score is null or they don't overlap
			if(previousScore==null || !annotation.overlaps(previousScore.getAnnotation())){
				//compute the score directly
				return processWindow(annotation);
			}
			
			double newScore=computeCount(annotation, previousScore);
			return new ScanStatisticScore(model, annotation, previousScore, newScore);
		}
	}

	/**
	 * True iff scan P value and annotation are equal
	 */
	@Override
	public boolean equals(Object o) {
		ScanStatisticScore otherScore = (ScanStatisticScore) o;
		if(scanPvalue != otherScore.getScanPvalue()) return false;
		if(!getAnnotation().equals(otherScore.getAnnotation())) return false;
		return true;
	}
	

	
}
