package broad.pda.feature.genome;


/**
 * @author engreitz
 *
 */
public class Chromosome {
	final private String name;
	// private ShortBEDReader maskedRegions; // TODO
	private long length; //TODO
	
	/**
	 * @param name   e.g. "chr12"
	 */
	public Chromosome(final String name) {
		this.name = name;
	}
	
	public boolean isSexChromosome() {
		return isSexChromosome(name);
	}
	
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
}
