/**
 * 
 */
package broad.pda.seq.protection;

/**
 * @author prussell
 */
public class PairedSampleWindowScore {

	private double backgroundLambda;
	private double signalLambda;
	private double backgroundCount;
	private double signalCount;
	private double countDifference;
	private double pValue;
	private double fdrCorrectedPvalue;
	
	/**
	 * Constructor initializes all fields to -99 except P value and corrected P value which are initialized to 99
	 */
	public PairedSampleWindowScore() {
		this.backgroundLambda = -99;
		this.signalLambda = -99;
		this.backgroundCount = -99;
		this.signalCount = -99;
		this.countDifference = -99;
		this.pValue = 99;
		this.fdrCorrectedPvalue = 99;
	}
	
	/**
	 * The names of the fields in the order they are rendered by toString()
	 */
	public static String SCORE_FIELDS = "backgroundLambda\tsignalLambda\tbackgroundCount\tsignalCount\tdifference\tpValue\tfdrCorrectedPvalue";
	
	@Override
	public String toString() {
		String rtrn = Double.valueOf(this.backgroundLambda).toString() + "\t";
		rtrn += Double.valueOf(this.signalLambda).toString() + "\t";
		rtrn += Double.valueOf(this.backgroundCount).toString() + "\t";
		rtrn += Double.valueOf(this.signalCount).toString() + "\t";
		rtrn += Double.valueOf(this.countDifference).toString() + "\t";
		rtrn += Double.valueOf(this.pValue) + "\t";
		rtrn += Double.valueOf(this.fdrCorrectedPvalue);
		return rtrn;
	}	
	
	/**
	 * Get lambda for background library
	 * @return Lambda for background library
	 */
	public double getBackgroundLambda() {return this.backgroundLambda;}
	
	/**
	 * Get lambda for signal library
	 * @return Lambda for signal library
	 */
	public double getSignalLambda() {return this.signalLambda;}
	
	/**
	 * Get read count for background library
	 * @return Read count for background library
	 */
	public double getBackgroundCount() {return this.backgroundCount;}
	
	/**
	 * Get read count for signal library
	 * @return Read count for signal library
	 */
	public double getSignalCount() {return this.signalCount;}
	
	/**
	 * Get difference in read counts background minus signal
	 * @return Background read count minus signal read count
	 */
	public double getCountDifference() {return this.countDifference;}
	
	/**
	 * Get window P value
	 * @return The P value
	 */
	public double getPvalue() {return this.pValue;}
	
	/**
	 * Get window FDR corrected P value
	 * @return The corrected P value
	 */
	public double getFdrCorrectedPvalue() {return this.fdrCorrectedPvalue;}
	
	/**
	 * Set lambda for background library
	 * @param lambda Lambda for background library
	 */
	public void setBackgroundLambda(double lambda) {this.backgroundLambda = lambda;}
	
	/**
	 * Set lambda for signal library
	 * @param lambda Lambda for signal library
	 */
	public void setSignalLambda(double lambda) {this.signalLambda = lambda;}
	
	/**
	 * Set read count for background libary and recalculate count difference
	 * @param count Read count for background library
	 */
	public void setBackgroundCount(double count) {
		this.backgroundCount = count;
		if(this.signalCount > -99) this.countDifference = this.backgroundCount - this.signalCount;
	}
	
	/**
	 * Set read count for signal library and recalculate count difference
	 * @param count Read count for signal library
	 */
	public void setSignalCount(double count) {
		this.signalCount = count;
		if(this.backgroundCount > -99) this.countDifference = this.backgroundCount - this.signalCount;
	}
	
	/**
	 * Set window P value
	 * @param p The P value
	 */
	public void setPvalue(double p) {this.pValue = p;}
	
	/**
	 * Set window FDR corrected P value
	 * @param p The corrected P value
	 */
	public void setFdrCorrectedPvalue(double p) {this.fdrCorrectedPvalue = p;}
	
	
}
