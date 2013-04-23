package broad.pda.feature.genome;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import nextgen.core.annotation.Annotation;

import broad.core.sequence.SequenceRegion;


/**
 * Represents a chromosome in a genome
 * @author engreitz
 *
 */
public class Chromosome {
	final private String name;
	// private ShortBEDReader maskedRegions; // TODO
	private long length; //TODO
	Annotation centromere;
	List<AgpEntry> gaps;
	List<AgpEntry> clones;
	private AgpEntry shortArm;
	
	/**
	 * @param The name of the chromosome    e.g. "chr12"
	 */
	public Chromosome(final String name) {
		this.name = name;
		gaps = new ArrayList<AgpEntry>();
		clones = new ArrayList<AgpEntry>();
	}
	
	/**
	 * @param The name of the chromosome    e.g. "chr12"
	 * @param An AGP file describing the chromosome assembly
	 */
	public Chromosome(String name, String agpFile) throws Exception {
		this(name);
		length = 0;
		File source = new File(agpFile);
		BufferedReader br = new BufferedReader(new FileReader(source));
		String line;

		try {
			while((line = br.readLine()) != null) {
				if(line.startsWith("#") || line.trim().length() ==0){
					continue;
				}
				String[] lineSplit = line.split("\t");

				AgpEntry entry = AgpEntryFactory.createEntry(lineSplit);
				switch (entry.getType()) {
				case AgpEntry.CLONE_TYPE :
					clones.add(entry);
					break;
				case AgpEntry.GAP_TYPE :
					gaps.add(entry);
					break;
				case AgpEntry.CENTROMERE_TYPE :
					centromere = entry;
					break;
				case AgpEntry.SHORT_ARM_TYPE : 
					shortArm = entry;
					break;
				}
				length = Math.max(length, entry.getEnd());
			}
		}  finally {
			try {
				br.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}
	
	public boolean isSexChromosome() {
		return isSexChromosome(name);
	}
	
	public Annotation getShortArm() { return this.shortArm;}
	
	public static boolean isSexChromosome(String name) {
		return name.equals("chrX") || name.equals("chrY");
	}
	
	public boolean isAutosome() {
		return isAutosome(name);
	}
	
	public static boolean isAutosome(String name) {
		boolean result = false;
		try {
			Integer.parseInt(getChromosomeNumber(name));
			result = true;
		} catch (NumberFormatException n) {
			// do nothing
		}
		return result;
	}
	
	public String getName() {
		return name;
	}
	
	public String getChromosomeNumber() {
		return getChromosomeNumber(name);
	}
	
	public static String getChromosomeNumber(String name) {
		String num = name.substring(3);
		return num;
	}
	
	public static int compareNames(String name1, String name2) {
		return new Chromosome(name1).compareTo(new Chromosome(name2));
	}
	
	public int compareTo(Object o) {
		if (!this.getClass().equals(o.getClass())) {
			throw new IllegalArgumentException("Cannot compare across classes.");
		}
		
		Chromosome other = (Chromosome) o;
		if (this.getName().equals(other.getName())) return 0;
	
		boolean isAutosome = this.isAutosome();
		boolean otherAutosome = other.isAutosome();
		String thisNum = this.getChromosomeNumber();
		String otherNum = other.getChromosomeNumber();
		
		int result;
		if (isAutosome && otherAutosome) {
			result = -1;
		} else if (!isAutosome() && otherAutosome) {
			result = 1;
		} else if (isAutosome && otherAutosome) {
			result = new Integer(thisNum).compareTo(new Integer(otherNum));
		} else {
			result = thisNum.compareTo(otherNum);
		}
		return result;
	}
	
	public long getUngappedSize() {
		long totalGaps = shortArm != null ? shortArm.getLength() : 0;
		
		Iterator<AgpEntry> gapIt = gaps.iterator();
		while(gapIt.hasNext()) {
			totalGaps += gapIt.next().getLength();
		}
		
		return length - totalGaps;
	}
	
	
	public static class AgpEntry extends SequenceRegion{

		static public final int CLONE_TYPE       = 0;
		static public final int CENTROMERE_TYPE  = 1;
		static public final int GAP_TYPE         = 2;
		static public final int CONTIG_TYPE      = 4;
		static public final int OTHER_TYPE		 = 3;
		public static final int SHORT_ARM_TYPE   = 5;
		private boolean reversedOrientation;
		private int type;
		private int number;
		
		public AgpEntry(String parentSequence) {
			super(parentSequence);
		}
		

		protected void setInReverseOrientation(boolean isInreverseOrientation) {
			this.reversedOrientation = isInreverseOrientation;
		}
		public boolean inReversedOrientation() {
			return reversedOrientation;
		}
			
		public void setType(int type) {
			this.type = type;
		}
		
		public int getType() {
			return type;
		}
		
		public int getLength() { return super.getLength() + 1;}
		
		public String toString() {
			StringBuffer buf = new StringBuffer("chr"+getContainingSequenceId());
			
			buf.append("\t")
				.append(getStart())
				.append("\t")
				.append(getEnd())
				.append("\t")
				.append(getNumber())
				.append("\t");
			
			switch (getType()) {
			case AgpEntry.CLONE_TYPE:
				buf.append("F\t");
				break;
			case AgpEntry.CENTROMERE_TYPE:
			case AgpEntry.SHORT_ARM_TYPE:
			case AgpEntry.GAP_TYPE:
				buf.append("N\t");
				break;
			case AgpEntry.CONTIG_TYPE:
				buf.append("W\t");
				break;
			default:
				buf.append("O\t");
				break;
			}
			
			buf.append(getName())
				.append("\t")
				.append("0\t")
				.append(getLength() + 1)
				.append("\t+");
			
			return buf.toString();
		}


		public int getNumber() {
			return number;
		}


		public void setNumber(int number) {
			this.number = number;
		}

	}
	
	public static class AgpEntryFactory {

		public static AgpEntry createEntry(String [] rawInfo) {
			AgpEntry entry = null;
			String parentName = rawInfo[0];
			
			entry = new AgpEntry(parentName);
			entry.setName(rawInfo[5]);
			entry.setInReverseOrientation(rawInfo.length == 9 && "-".equals(rawInfo[8]));	
			entry.setStart(Integer.parseInt(rawInfo[1]));
			entry.setEnd(Integer.parseInt(rawInfo[2]));
			entry.setNumber(Integer.parseInt(rawInfo[3]));
			if(rawInfo[0].startsWith("Un") || rawInfo[0].startsWith("un")) {
				entry.setChromosome("Un");
			} else {
				entry.setChromosome(rawInfo[0].length() < 4 ? rawInfo[0] : rawInfo[0].substring(3)); //try to handle both chrNN and NN notations for chromosomes
				//System.out.println("Entry's chromosome " + entry.getChromosome());
				if(entry.getChromosome().startsWith("0")) {
					entry.setChromosome(entry.getChromosome().substring(1));
				}
			}
			if ("F".equals(rawInfo[4])) {			
				entry.setName(rawInfo[5]);
				entry.setType(AgpEntry.CLONE_TYPE);
			} else if ("centromere".equalsIgnoreCase(rawInfo[6])) {
				entry.setName("centromere");
				entry.setType(AgpEntry.CENTROMERE_TYPE);	
			} else if ("short_arm".equalsIgnoreCase(rawInfo[6])) {
				entry.setName("short_arm");
				entry.setType(AgpEntry.SHORT_ARM_TYPE);	
			} else if ("N".equals(rawInfo[4])) {
				entry.setType(AgpEntry.GAP_TYPE);
			} else if ("W".equalsIgnoreCase(rawInfo[4])) {
				entry.setType(AgpEntry.CONTIG_TYPE);
			} else {
				entry.setName("other");
				entry.setType(AgpEntry.OTHER_TYPE);
				//System.out.print("Can't determine what this AGP entry is type<" + rawInfo[4] +"> name <  rawInfo<"+rawInfo[5]);
			}
			
			//System.out.println("created entry " + entry);
			return entry;
		}


	}

}
