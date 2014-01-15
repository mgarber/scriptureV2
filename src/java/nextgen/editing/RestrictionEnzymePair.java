package nextgen.editing;

/**
 * A pair of restriction enzymes with a left and right, for double cuts
 * @author prussell
 */
public class RestrictionEnzymePair {
	
	private RestrictionEnzyme left;
	private RestrictionEnzyme right;
	
	/**
	 * @param leftEnzyme Enzyme for left cut
	 * @param rightEnzyme Enzyme for right cut
	 */
	public RestrictionEnzymePair(RestrictionEnzyme leftEnzyme, RestrictionEnzyme rightEnzyme) {
		left = leftEnzyme;
		right = rightEnzyme;
	}
	
	public RestrictionEnzyme getLeft() {
		return left;
	}
	
	public RestrictionEnzyme getRight() {
		return right;
	}
	
}
