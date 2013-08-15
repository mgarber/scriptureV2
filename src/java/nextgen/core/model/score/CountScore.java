package nextgen.core.model.score;

import org.apache.log4j.Logger;

import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.AnnotationCollection;
import nextgen.core.annotation.BasicAnnotation;
import nextgen.core.model.AlignmentModel;

/**
 * @author engreitz
 * This class represents a simple count scoring function for the DataAlignmentModel.
 */
public class CountScore extends WindowScore.AbstractWindowScore implements Comparable<CountScore> {

	static Logger logger = Logger.getLogger(CountScore.class.getName());
	public static double DEFAULT_REGION_TOTAL = -1.0;
	
	private double count;                // use a double to allow for weighted read counting
	private double total;                // whole "genome"
	private double regionTotal;          // whole chromosome / gene / region under consideration / etc.

	public CountScore(Annotation a) {
		super(a);
	}
	
	public CountScore(CountScore other) {
		super(other.getAnnotation());
		this.count = other.count;
		this.total = other.total;
		this.regionTotal = other.regionTotal;
	}
	
	public CountScore(CountScore other, Annotation newAnnotation, double newCount) {
		super(newAnnotation);
		this.count = newCount;
		this.total = other.total;
		this.regionTotal = other.regionTotal;
	}
	
	public CountScore(AnnotationCollection<? extends Annotation> model, Annotation annotation) {
		this(model, annotation, false);
	}
	
	public CountScore(AnnotationCollection<? extends Annotation> model, Annotation annotation, double regionTotal) {
		this(model, annotation, regionTotal, false);
	}

	
	public CountScore(AnnotationCollection<? extends Annotation> model, Annotation annotation, boolean fullyContained) {
		super(annotation);
		setCount(model.getCount(annotation, fullyContained));
		setTotal(model.getGlobalCount());
		getAnnotation().setScore(getCount());
	}
	
	public CountScore(AnnotationCollection<? extends Annotation> model, Annotation annotation, double regionTotal, boolean fullyContained) {
		this(model, annotation, fullyContained);
		setRegionTotal(regionTotal);
	}
	
	public CountScore(Annotation annotation, double count, double regionTotal, double total) {
		super(annotation);
		setCount(count);
		setTotal(total);
		setRegionTotal(regionTotal);
	}
	
	public double getCount() { 
		return count; 
	}
	
	@Override
	public double getScore() { 
		return getCount();
	}
	
	public double getTotal() { return total; }
	public double getRegionTotal() { return regionTotal; }
	public double getRPKM() { 
		return asRPKM(count, total, getAnnotation().length()); 
	}
	
	public double getRPK() { 
		return count/(getAnnotation().length() / 1000.0); 
	}
	
	public void setCount(double count) {
		this.count = count;
	}
	
	public void setTotal(double total) {
		this.total = total;
	}
	
	public void setRegionTotal(double regionTotal) {
		this.regionTotal = regionTotal;
	}
	
	public static double asRPKM(double count, double total, int windowSize) {
		return count / total * 1000000.0 / (windowSize / 1000.0);
	}
	
	public String toString() {
		annotation.setScore(getScore());
		return annotation.toBED() + "\t" + 
				getCount() + "\t" + 
				getRPKM() + "\t" + 
				getRegionTotal() + "\t" + 
				getTotal();
	}
	
	/**
	 * True iff count and annotation are equal
	 */
	@Override
	public boolean equals(Object o) {
		ScanStatisticScore otherScore = (ScanStatisticScore) o;
		if(count != otherScore.getCount()) return false;
		if(!getAnnotation().equals(otherScore.getAnnotation())) return false;
		return true;
	}
	
	/**
	 * First compare counts
	 * Then compare annotations
	 */
	@Override
	public int compareTo(CountScore o) {
		
		// First compare counts
		double otherCount = o.getCount();
		if(count < otherCount) return -1;
		if(count > otherCount) return 1;
		
		// Then compare annotations
		return getAnnotation().compareTo(o.getAnnotation());
	}

	
	
	public static class Processor extends WindowProcessor.AbstractProcessor<CountScore> {
		protected AnnotationCollection<? extends Annotation> model;
		protected double regionTotal = DEFAULT_REGION_TOTAL;
		protected boolean skipInit = false;
		private boolean fullyContainedReads;
		
		public Processor(AnnotationCollection<? extends Annotation> model) {
			this(model, false);
		}
		
		
		public Processor(AnnotationCollection<? extends Annotation> model, boolean fullyContained) {
			this.model = model;
			this.fullyContainedReads = fullyContained;
		}
		
		public Processor(AnnotationCollection<? extends Annotation> model, boolean skipInit, boolean fullyContained) {
			this(model, fullyContained);
			this.skipInit = skipInit;
		}
		
		
		public CountScore processWindow(Annotation annotation) {
			return new CountScore(model, annotation, regionTotal, fullyContainedReads);
		}
		
		public CountScore processWindow(Annotation annotation, CountScore previousScore) {
			// This previousScore method is broken because getCountExcludingRegion filters out reads
			// that overlap both the requested interval as well as the excluded interval (which should be counted)
			// Also, it ends up using twice as many cache updates.  
			// Fix these before re-applying. -JE
			
			return processWindow(annotation);
			
			/*
			//if the previous score is null or they don't overlap by a significant amount (or at all)
			if (previousScore==null || annotation.intersect(previousScore.getAnnotation()).length() <= annotation.size()/2){
				//compute the score directly
				return processWindow(annotation);
			} 
	
			double newScore=computeCount(annotation, previousScore);
			return new CountScore(previousScore, annotation, newScore);
			*/
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
			
			if (previousScore.getAnnotation().equals(nextRegion)) {
				return previousScore.getCount();
			}
			
			double subtractVal = 0, addVal = 0;
			if (!nextRegion.contains(previousScore.getAnnotation())) {
				subtractVal = model.getCountExcludingRegion(previousScore.getAnnotation().minus(nextRegion), nextRegion);	
			}
			if (!previousScore.getAnnotation().contains(nextRegion)) {
				addVal=model.getCountExcludingRegion(nextRegion.minus(previousScore.getAnnotation()), previousScore.getAnnotation());
			}
			return (previousScore.getCount()-subtractVal)+addVal;
		}
		
		public void initRegion(Annotation region) {
			if (region != null && !skipInit) {
				regionTotal = model.getCount(region);
			}
		}
		
	}
	
	
	public static class Factory implements nextgen.core.general.TabbedReader.Factory<CountScore> {
		public CountScore create(String[] rawFields) {
			BasicAnnotation a = new BasicAnnotation.Factory().create(rawFields);
			return new CountScore(a, Double.parseDouble(rawFields[12]), Double.parseDouble(rawFields[14]), Double.parseDouble(rawFields[15]));
		}
	}


	

}
