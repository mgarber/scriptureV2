package broad.core.primer3;

import java.io.BufferedWriter;
import java.io.IOException;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.List;

import broad.core.sequence.Sequence;

public class PrimerPair implements Comparable{
	private String leftPrimer;
	private String rightPrimer;
	private float leftPrimerTM;
	private float rightPrimerTM;
	private int leftPrimerPosition;
	private int rightPrimerPosition;
	private String productSequence;
	private String primerPairId;
	private int productSize;
	private String comment;
	private Primer3SequenceInputTags usedInput;
	private Primer3Configuration configUsed;
	private float primerPairPenalty;
	private float leftPrimerPenalty;
	private float rightPrimerPenalty;

	static private NumberFormat nf = NumberFormat.getPercentInstance();


	public PrimerPair(String [] rawData) {
		int i = 0;
		primerPairId = rawData[i++];
		leftPrimerPosition = Integer.parseInt(rawData[i++]);
		leftPrimer         = rawData[i++];
		rightPrimerPosition= Integer.parseInt(rawData[i++]);
		rightPrimer        = rawData[i++];
		leftPrimerTM       = Float.parseFloat(rawData[i++]);
		rightPrimerTM      = Float.parseFloat(rawData[i++]);
		productSize        = Integer.parseInt(rawData[i++]);
		comment            = rawData[i++];
	}
	
	public String toString(boolean printExpectedProduct) {
		StringBuffer buf = new StringBuffer(primerPairId);
		buf.append("\t").append(leftPrimerPosition)
			.append("\t").append(leftPrimer)
			.append("\t").append(rightPrimerPosition)
			.append("\t").append(rightPrimer)
			.append("\t").append(leftPrimerTM)
			.append("\t").append(rightPrimerTM)
			.append("\t").append(productSize)
			.append("\t").append(comment == null ? "" : comment)
			.append("\t").append(productSequence == null ? "-" : nf.format(Sequence.computeGCContent(productSequence)));
			if(printExpectedProduct) {
				buf.append("\t").append(productSequence == null ? "-" : productSequence);
			}
		return buf.toString();
	}
	
	public String toString() { return toString(false); }
	
	public static String getPrimerFieldsAsString() {
		StringBuffer buf = new StringBuffer("Primer Pair Id");
		buf.append("\t").append("Left Primer Position")
			.append("\t").append("Left Primer")
			.append("\t").append("Right Primer Position")
			.append("\t").append("Right Primer")
			.append("\t").append("Left PrimerTM")
			.append("\t").append("Right PrimerTM")
			.append("\t").append("Product Size")
			.append("\t").append("Comment")
			.append("\t").append("Product GC%");
		return buf.toString();
	}
	

	/**
	 * The string hash code from the string representation
	 * @return
	 */
	@Override
	public int hashCode() {
		return (this.getLeftPrimer() + "_" + this.getRightPrimer()).hashCode();
	}
	
	/**
	 * Equal if left and right primers sequences are the same
	 */
	@Override
	public boolean equals(Object o) {
		if(o.getClass() != this.getClass()) return false;
		PrimerPair p = (PrimerPair)o;
		if(this.getLeftPrimer().equals(p.getLeftPrimer()) && this.getRightPrimer().equals(p.getRightPrimer())) return true;
		return false;
	}
	
	/**
	 * Get the fields in the order accepted by the constructor
	 * @return
	 */
	public String getPrimerFieldsAsStringForConstructor() {
		String data = this.getPrimerPairId();
		data += " ";
		data = Integer.valueOf(this.getLeftPrimerPosition()).toString();
		data += " ";
		data +=  this.getLeftPrimer();
		data += " ";
		data +=  Integer.valueOf(this.getRightPrimerPosition()).toString();
		data += " ";
		data +=  this.getRightPrimer();
		data += " ";
		data +=  Float.valueOf(this.getLeftPrimerTM()).toString();
		data += " ";
		data += Float.valueOf(this.getRightPrimerTM()).toString();
		data += " ";
		data +=  Integer.valueOf(this.getProductSize()).toString();
		data += " ";
		data += this.getComment();
		return data;
	}
	
	/**
	 * Get basic information about the primer pair in a tab-delimited string
	 * @return
	 */
	public String getPrimerPairInformation() {
		String data = Integer.valueOf(this.getLeftPrimerPosition()).toString();
		data += "\t";
		data +=  this.getLeftPrimer();
		data += "\t";
		data +=  Float.valueOf(this.getLeftPrimerTM()).toString();
		data += "\t";
		data +=  Float.valueOf(this.getLeftPrimerPenalty()).toString();
		data += "\t";
		data +=  Integer.valueOf(this.getRightPrimerPosition()).toString();
		data += "\t";
		data +=  this.getRightPrimer();
		data += "\t";
		data += Float.valueOf(this.getRightPrimerTM()).toString();
		data += "\t";
		data +=  Float.valueOf(this.getRightPrimerPenalty()).toString();
		data += "\t";
		data +=  Integer.valueOf(this.getProductSize()).toString();
		data += "\t";
		data +=  Float.valueOf(this.getPrimerPairPenalty()).toString();

		return data;
	}

	/**
	 * The names of the fields reported by getPrimerPairInformation()
	 * @return
	 */
	public static String getPrimerPairInformationFieldNames() {
		return "LeftPrimerPos\tLeftPrimer\tLeftPrimerTM\tLeftPrimerPenalty\tRightPrimerPos\tRightPrimer\tRightPrimerTM\tRightPrimerPenalty\tProductSize\tPrimerPairPenalty";
	}
	
	public String getComment() {
		return comment;
	}

	protected void setComment(String comment) {
		this.comment = comment;
	}

	public PrimerPair() {
		super();
		// TODO Auto-generated constructor stub
	}

	public String getLeftPrimer() {
		return leftPrimer;
	}

	public int getLeftPrimerPosition() {
		return leftPrimerPosition;
	}

	public float getLeftPrimerTM() {
		return leftPrimerTM;
	}

	public String getPrimerPairId() {
		String id = null;
		if(primerPairId != null) {
			id =  primerPairId;
		} else if (usedInput != null && usedInput.getSequence() != null) {
			id = usedInput.getSequence().getId();
		}
		return id;
	}

	public String getProductSequence() {
		return productSequence;
	}
	
	public String getRightPrimer() {
		return rightPrimer;
	}

	public int getRightPrimerPosition() {
		return rightPrimerPosition;
	}

	public float getRightPrimerTM() {
		return rightPrimerTM;
	}

	public Primer3SequenceInputTags getUsedInput() {
		return usedInput;
	}

	public void setUsedInput(Primer3SequenceInputTags usedInput) {
		this.usedInput = usedInput;
	}
	
	public Primer3Configuration getConfigurationUsed() {
		return configUsed;
	}
	
	public void setConfigurationUsed(Primer3Configuration usedConfig) {
		this.configUsed = usedConfig;
		
	}
	
	protected void setLeftPrimer(String leftPrimer) {
		this.leftPrimer = leftPrimer;
	}

	public void setLeftPrimerPosition(int leftPrimerPosition) {
		this.leftPrimerPosition = leftPrimerPosition;
	}

	public void setPrimerPairId(String primerPairId) {
		this.primerPairId = primerPairId;
	}

	protected void setProductSequence(String productSequence) {
		this.productSequence = productSequence;
	}

	protected void setRightPrimer(String rightPrimer) {
		this.rightPrimer = rightPrimer;
	}

	public void setRightPrimerPosition(int rightPrimerPosition) {
		this.rightPrimerPosition = rightPrimerPosition;
	}

	protected void setRightPrimerTM(float rightPrimerTM) {
		this.rightPrimerTM = rightPrimerTM;
	}

	public int getProductSize() {
		return productSize;
	}

	protected void setLeftPrimerTM(float leftPrimerTM) {
		this.leftPrimerTM = leftPrimerTM;
	}

	protected void setProductSize(int productSize) {
		this.productSize = productSize;
	}

	protected void setRightPrimerTM(Float rightPrimerTM) {
		this.rightPrimerTM = rightPrimerTM;
	}

	public boolean isProductRestricted(String recognizedSequence) {
		if(productSequence == null) {
			System.out.println("Returning False because product is null");
			return false;  
		}
		//System.out.println("recognized sequence: "+recognizedSequence);
		return productSequence.contains(recognizedSequence);
	}
	
	public boolean hasPrimers() {
		return leftPrimer != null && leftPrimer.length() > 0;
	}

	public void writePrimer3Record(BufferedWriter bw, Primer3IO p3io) throws IOException {
		p3io.writeRecord(bw, getUsedInput(), getConfigurationUsed());
		
	}
	
	/**
	 * Write the whitespace separated data line in the order that the constructor uses
	 */
	public void print() {
		String data = Integer.valueOf(this.getLeftPrimerPosition()).toString();
		data += " ";
		data +=  this.getLeftPrimer();
		data += " ";
		data +=  Integer.valueOf(this.getRightPrimerPosition()).toString();
		data += " ";
		data +=  this.getRightPrimer();
		data += " ";
		data +=  Float.valueOf(this.getLeftPrimerTM()).toString();
		data += " ";
		data += Float.valueOf(this.getRightPrimerTM()).toString();
		data += " ";
		data +=  Integer.valueOf(this.getProductSize()).toString();
		data += " ";
		data += this.getComment();
		System.out.println(data);
	}
	
	
	public List<Sequence> getPrimerSequences() {
		ArrayList<Sequence> result = new ArrayList<Sequence>();
		if(hasPrimers()) {
			Sequence leftPrimerSeq = new Sequence(getPrimerPairId() + "_left");
			leftPrimerSeq.setSequenceBases(leftPrimer);
			result.add(leftPrimerSeq);
			
			Sequence rightPrimerSeq = new Sequence(getPrimerPairId() + "_right");
			rightPrimerSeq.setSequenceBases(rightPrimer);
			result.add(rightPrimerSeq);
		}
		
		return result;
	}

	public void setPrimerPairPenalty(float f) {
		this.primerPairPenalty = f;	
	}
	
	public float getPrimerPairPenalty(){
		return this.primerPairPenalty;
	}

	public void setFields(String primer, String tag, String value) {
		if(primer.equalsIgnoreCase("LEFT")){
			setLeftValue(tag, value);
		}
		else if(primer.equalsIgnoreCase("RIGHT")){
			setRightValue(tag, value);
		}
		else if(primer.equalsIgnoreCase("PAIR")){
			setPairValue(tag, value);
		}
		else if(primer.equalsIgnoreCase("INTERNAL")){
			setInternalValue(tag, value);
		}
		else if(primer.equalsIgnoreCase("MIN") && tag.equalsIgnoreCase("PRIME_OVERLAP_OF_JUNCTION_")){}
		else{System.err.println(primer+" "+tag+" "+value); throw new IllegalArgumentException("Primer oligo has to be Left, Right, Pair, or Internal");}
					
	}

	private void setLeftValue(String tag, String value) {
		if(tag.equalsIgnoreCase("")){
			//then set primer position
			setLeftPrimerPosition(Integer.parseInt(value.split(",")[0]));
		}
		else if(tag.equalsIgnoreCase("PENALTY_")){
			setLeftPrimerPenalty(Float.parseFloat(value));
		}
		else if(tag.equalsIgnoreCase("SEQUENCE_")){
			setLeftPrimer(value);
		}
		else if(tag.equalsIgnoreCase("TM_")){
			setLeftPrimerTM(Float.parseFloat(value));
		}
		else if(tag.equalsIgnoreCase("GC_")){}
		else if(tag.equalsIgnoreCase("SELF_ANY_")){}
		else if(tag.equalsIgnoreCase("SELF_END_")){}
		else if(tag.equalsIgnoreCase("LIBRARY_MISPRIMING_")){}
		else if(tag.equalsIgnoreCase("END_STABILITY_")){}
		else if(tag.equalsIgnoreCase("COMPL_ANY_")){}
		else if(tag.equalsIgnoreCase("COMPL_END_")){}
		else{}
	}
	
	private void setRightValue(String tag, String value) {
		if(tag.equalsIgnoreCase("")){
			//then set primer position
			setRightPrimerPosition(Integer.parseInt(value.split(",")[0]));
		}
		else if(tag.equalsIgnoreCase("PENALTY_")){
			setRightPrimerPenalty(Float.parseFloat(value));
		}
		else if(tag.equalsIgnoreCase("SEQUENCE_")){
			setRightPrimer(value);
		}
		else if(tag.equalsIgnoreCase("TM_")){
			setRightPrimerTM(Float.parseFloat(value));
		}
		else if(tag.equalsIgnoreCase("GC_")){}
		else if(tag.equalsIgnoreCase("SELF_ANY_")){}
		else if(tag.equalsIgnoreCase("SELF_END_")){}
		else if(tag.equalsIgnoreCase("LIBRARY_MISPRIMING_")){}
		else if(tag.equalsIgnoreCase("END_STABILITY_")){}
		else if(tag.equalsIgnoreCase("COMPL_ANY_")){}
		else if(tag.equalsIgnoreCase("COMPL_END_")){}
		else{}
	}
	
	private void setPairValue(String tag, String value) {
		if(tag.equalsIgnoreCase("PENALTY_")){
			setPrimerPairPenalty(Float.parseFloat(value));
		}
		else if(tag.equalsIgnoreCase("PRODUCT_SIZE_")){
			setProductSize(Integer.parseInt(value));
		}
		else{}
	}

	private void setLeftPrimerPenalty(float penalty) {
		leftPrimerPenalty=penalty;	
	}
	private void setRightPrimerPenalty(float penalty) {
		rightPrimerPenalty=penalty;	
	}
	
	
	public float getLeftPrimerPenalty(){
		return this.leftPrimerPenalty;
	}
	
	public float getRightPrimerPenalty(){
		return this.rightPrimerPenalty;
	}
	
	private void setInternalValue(String tag, String value) {
		// TODO Add support for internal oligos later
		
	}

	
	public int compareTo(Object o) {
		PrimerPair other=(PrimerPair)o;
		//First sort by penalty
		if(other.getPrimerPairPenalty()!=this.getPrimerPairPenalty()){return new Double(this.getPrimerPairPenalty()).compareTo(new Double(other.getPrimerPairPenalty()));}
		//then by left position
		if(other.getLeftPrimerPosition()!=this.getLeftPrimerPosition()){return this.getLeftPrimerPosition()-other.getLeftPrimerPosition();}
		//then by right position
		if(other.getRightPrimerPosition()!=this.getRightPrimerPosition()){return this.getRightPrimerPosition()-other.getRightPrimerPosition();}
		
		//then by left 
		//then by right
		
		return 0;
	}


}
