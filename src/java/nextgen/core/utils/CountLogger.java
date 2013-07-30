/**
 * 
 */
package nextgen.core.utils;

import java.text.DecimalFormat;

import org.apache.log4j.Logger;

/**
 * @author prussell
 * Log progress of an incrementing count
 */
public class CountLogger {
	
	private int totalCount;
	private int totalMessages;
	private int numDone;
	private int step;
	private static Logger logger = Logger.getLogger(CountLogger.class.getName());
	private DecimalFormat decimalFormat;
	
	/**
	 * @param overallTotal The eventual total for the incrementing count
	 */
	public CountLogger(int overallTotal) {
		this(overallTotal, 100);
	}
	
	/**
	 * @param overallTotal The eventual total for the incrementing count
	 * @param numMessages The total number of messages to print
	 */
	public CountLogger(int overallTotal, int numMessages) {
		totalCount = overallTotal;
		totalMessages = numMessages;
		numDone = 0;
		step = Math.max(1, totalCount / totalMessages);
		decimalFormat = new DecimalFormat("#.##");
	}
	
	/**
	 * Increment the count and print progress message if appropriate
	 */
	public void advance() {
		numDone++;
		if(numDone % step == 0) {
			String pct = decimalFormat.format(100 * (double)numDone / totalCount);
			logger.info("Finished " + pct + "%");
		}
	}
	
}
