/**
 * 
 */
package broad.pda.methylation;

import java.util.*;
import broad.core.annotation.*;

/**
 * @author engreitz
 *
 */
public final class PadlockProbe {
	private final LightweightGenomicAnnotation target, probe;
	private final String leftCapture, rightCapture, sequence, name;
	private final int leftTm, rightTm;
	private final boolean strand; // TRUE = (W)atson, FALSE = (C)rick.
	
	PadlockProbe(String[] fields) {
		target = BasicLightweightAnnotation.createFromUCSC(fields[0]);
		int offset = Integer.parseInt(fields[2]);
		probe = new BasicLightweightAnnotation(target.getChromosome(), 
							target.getStart() + offset + Integer.parseInt(fields[3]),
							target.getStart() + offset + Integer.parseInt(fields[4]));
		leftCapture = fields[5];
		leftTm = Integer.parseInt(fields[6]);
		rightCapture = fields[7];
		rightTm = Integer.parseInt(fields[8]);
		strand = fields[10].equals("W");
		sequence = fields[11];
		name = fields[0] + "_" + fields[3];
	}
	
	public LightweightGenomicAnnotation getTarget() { return target; }
	public LightweightGenomicAnnotation getProbe() { return probe; }
	public String getLeftCapture() { return leftCapture; }
	public String getRightCapture() { return rightCapture; }
	public String getCapturedSequence() { return sequence; }
	public int getLeftTm() { return leftTm; }
	public int getRightTm() { return rightTm; }
	public boolean getStrand() { return strand; }
	public String toString() { return name; }
	
	public boolean equals(PadlockProbe other) {
		return leftCapture.equals(other.leftCapture) & rightCapture.equals(other.rightCapture);
	}
}
