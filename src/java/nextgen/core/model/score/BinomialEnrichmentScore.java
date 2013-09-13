package nextgen.core.model.score;

import broad.pda.seq.segmentation.AlignmentDataModelStats;
import net.sf.samtools.util.CloseableIterator;
import nextgen.core.alignment.Alignment;
import nextgen.core.annotation.Annotation;
import nextgen.core.coordinatesystem.CoordinateSpace;
import nextgen.core.model.AlignmentModel;
import jsc.distributions.Binomial;

public class BinomialEnrichmentScore extends CountScore {

	private double Pvalue;
	private CoordinateSpace sampleCoordSpace;
	private CoordinateSpace ctrlCoordSpace;
	private double ctrlCount;              
	private double sampleRegionCount;
	private double ctrlRegionCount; //Total bases in gene/region/chrom in control
	
	public BinomialEnrichmentScore(AlignmentModel sample, AlignmentModel ctrl, Annotation a) {
		super(sample, a);
		sampleCoordSpace = sample.getCoordinateSpace();
		ctrlCoordSpace = ctrl.getCoordinateSpace();
		
		if (!sampleCoordSpace.getChromosomeNames().equals(ctrlCoordSpace.getChromosomeNames())) {
			throw new IllegalArgumentException("Sample coordinate space must match control coordinate space");
		}
		
		setCtrlCount(ctrl.getCount(a));
		setCtrlRegionCount(ctrl.getCount(a));
		setSampleRegionCount(sample.getCount(a));
	
		try {
			setPvalue(calculatePVal(getSampleCount(), getCtrlCount(), getSampleRegionCount(), getCtrlRegionCount()));
			} catch(Exception e){
				logger.info("Cound not set P value for annotation " + a.getName());
			}
			getAnnotation().setScore(getPvalue());
	}
	
	public BinomialEnrichmentScore(AlignmentModel sample, Annotation a) {
		this(sample, sample, a);
	}

	public BinomialEnrichmentScore(AlignmentModel sample, AlignmentModel ctrl, Annotation a, double sampleRegCount, double ctrlRegCount) {
		this(sample,ctrl,a);
		setSampleRegionCount(sampleRegCount);
		setCtrlRegionCount(ctrlRegCount);
		setPvalue(calculatePVal(getSampleCount(), getCtrlCount(), getSampleRegionCount(), getCtrlRegionCount()));
		
		getAnnotation().setScore(getScanPvalue());
	}
	
	public BinomialEnrichmentScore(AlignmentModel sample, AlignmentModel sample, Annotation annotation, BinomialEnrichmentScore previousScore, double newSampleCount,double newCtrlCount) {
		super(previousScore, annotation, newScore); //Set the new score without computing
		sampleCoordSpace = sample.getCoordinateSpace();
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
	
	/**
	 * @param a				Sample reads
	 * @param b	        	Control reads
	 * @param sampleRegionCounts	Sample reads in region
	 * @param ctrlRegionCounts	Control reads in region
	 * @return				Binomial P(X>a)
	 */
	public double calculatePVal(double a, double b, double sampleRegionCounts, double ctrlRegionCounts) {
		
		double p = sampleRegionCounts/ctrlRegionCounts;
		long n = (long) (a + b);
		Binomial C = new Binomial(n,p);
		double pval = 1 - C.cdf(a);
		return pval;

	}
	
	public static class Processor extends WindowProcessor.AbstractProcessor<BinomialEnrichmentScore> {
		protected AlignmentModel sample;
		protected AlignmentModel ctrl;
		protected double sampleRegionCount = DEFAULT_REGION_TOTAL;
		protected double ctrlRegionCount = DEFAULT_REGION_TOTAL;
		private boolean fullyContainedReads;
		
		public Processor(AlignmentModel sample, AlignmentModel ctrl) {
			this.sample = sample;
			this.ctrl = ctrl;
		}
		
		public Processor(AlignmentModel sample) {
			this.sample = sample;
			this.ctrl = sample;
		}
		
		public void initRegion(Annotation a){
			if (a != null){
				sampleRegionCount = sample.getCount(a);
				ctrlRegionCount = ctrl.getCount(a);
			}
		}
		
		public BinomialEnrichmentScore processWindow(Annotation a){
			return new BinomialEnrichmentScore(sample, ctrl, a, sampleRegionCount, ctrlRegionCount);
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
			double subtractVal=sample.getCountExcludingRegion(previousScore.getAnnotation().minus(nextRegion), nextRegion);
			double addVal=sample.getCountExcludingRegion(nextRegion.minus(previousScore.getAnnotation()), previousScore.getAnnotation());
			return (previousScore.getCount()-subtractVal)+addVal;
		}
		
		@Override
		public BinomialEnrichmentScore processWindow(Annotation annotation, BinomialEnrichmentScore previousScore) {
			//if the previous score is null or they don't overlap
			if(previousScore==null || !annotation.overlaps(previousScore.getAnnotation())){
				//compute the score directly
				return processWindow(annotation);
			}
				
			double newScore=computeCount(annotation, previousScore);
			return new BinomialEnrichmentScore(sample, annotation, previousScore, newScore);
			}
		
	}
	
	public double getAverageCoverage(AlignmentModel data, CoordinateSpace coordSpace) { 
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
	
	public void setPvalue(double scanPvalue) { this.Pvalue = scanPvalue; }
	public double getPvalue() { return Pvalue; }

	public void setSampleCount(double d) { setCount(d); }
	public void setCtrlCount(double d) { ctrlCount = d; }

	public void setSampleRegionCount(double d) { sampleRegionCount = d; }
	public void setCtrlRegionCount(double d) { ctrlRegionCount = d; }
	
	public double getSampleCount() { return getCount(); }
	public double getCtrlCount() { return ctrlCount; }

	public double getCtrlRegionCount() { return ctrlRegionCount; }
	public double getSampleRegionCount() { return sampleRegionCount; }

}
