package nextgen.synbio;

import broad.core.sequence.Sequence;

/**
 * @author prussell
 * A circular plasmid incorporating a deletion to a transcript
 * Plasmid sequence starts immediately 3' of a deletion and ends immediately 5' of the deletion
 * Includes only the actual transcript sequence, no special "plasmid stuff"
 */
public class DeletionPlasmid {
	
	private Sequence originalSequence;
	private int firstDeletedPos;
	private int lastDeletedPos;
	private String id;
	
	/**
	 * @param fullTranscriptSequence Original transcript
	 * @param firstDeletedPosition First position contained in deletion, in transcript coordinates
	 * @param lastDeletedPosition Last position contained in deletion (inclusive), in transcript coordinates
	 */
	public DeletionPlasmid(Sequence fullTranscriptSequence, int firstDeletedPosition, int lastDeletedPosition) {
		this(fullTranscriptSequence, firstDeletedPosition, lastDeletedPosition, generateID(fullTranscriptSequence.getId(), firstDeletedPosition, lastDeletedPosition));
	}
	
	/**
	 * @param fullTranscriptSequence Original transcript
	 * @param firstDeletedPosition First position contained in deletion, in transcript coordinates
	 * @param lastDeletedPosition Last position contained in deletion (inclusive), in transcript coordinates
	 * @param plasmidID Plasmid ID to assign
	 */
	public DeletionPlasmid(Sequence fullTranscriptSequence, int firstDeletedPosition, int lastDeletedPosition, String plasmidID) {
		originalSequence = fullTranscriptSequence;
		firstDeletedPos = firstDeletedPosition;
		lastDeletedPos = lastDeletedPosition;
		id = plasmidID;
	}
	
	/**
	 * @return ID
	 */
	public String getId() {
		return id;
	}
	
	protected static String generateID(String transcriptName, int firstDeletedPos, int lastDeletedPos) {
		return transcriptName + "_deletion_" + firstDeletedPos + "_" + lastDeletedPos;
	}
	
	/**
	 * @return The original sequence not incorporating the deletion
	 */
	public Sequence getOriginalSequence() {
		return originalSequence;
	}
	
	/**
	 * @return First position contained in deletion, in transcript coordinates
	 */
	public int getFirstDeletedPosition() {
		return firstDeletedPos;
	}
	
	/**
	 * @return Last position contained in deletion (inclusive), in transcript coordinates
	 */
	public int getLastDeletedPosition() {
		return lastDeletedPos;
	}
	
	/**
	 * @return The sequence of the plasmid incorporating the deletion
	 * Sequence starts immediately 3' of the deletion and ends immediately 5' of the deletion
	 */
	public Sequence getPlasmidSequence() {
		if(firstDeletedPos <= lastDeletedPos) {
			Sequence beforeDeletion = originalSequence.getSubSequence(null, 0, firstDeletedPos);
			Sequence afterDeletion = originalSequence.getSubSequence(id, lastDeletedPos + 1, originalSequence.getLength());
			afterDeletion.append(beforeDeletion.getSequenceBases());
			afterDeletion.uppercase();
			return afterDeletion;
		}
		Sequence afterDeletion = originalSequence.getSubSequence(id, lastDeletedPos + 1, firstDeletedPos);
		afterDeletion.uppercase();
		return afterDeletion;
	}
	
	/**
	 * @return Plasmid size
	 */
	public int getPlasmidSize() {
		return getPlasmidSequence().getLength();
	}
	
	/**
	 * @return The size of the deletion
	 */
	public int getDeletionSize() {
		return originalSequence.getLength() - getPlasmidSequence().getLength();
	}
	
	/**
	 * Get the sequence surrounding the deletion
	 * @param numBasesUpstream Number of bases upstream of deletion to include
	 * @param numBasesDownstream Number of bases downstream of deletion to include
	 * @param includeDeletionAsN Whether to include the deletion region as Ns
	 * @return Flanking sequence surrounding the deletion
	 */
	public Sequence getSequenceFlankingDeletion(int numBasesUpstream, int numBasesDownstream, boolean includeDeletionAsN) {
		Sequence plasmidSequence = getPlasmidSequence();
		String plasmidBases = plasmidSequence.getSequenceBases();
		int plasmidSize = plasmidBases.length();
		String upstream = plasmidBases.substring(plasmidSize - numBasesUpstream);
		String downstream = plasmidBases.substring(0, numBasesDownstream);
		String flankingSequence = upstream;
		if(includeDeletionAsN) {
			String deletion = "";
			for(int i = 0; i < getDeletionSize(); i++) {
				deletion += "N";
			}
			flankingSequence += deletion;
		}
		flankingSequence += downstream;
		String rtrnID = getId() + "_flanking_sequence_" + numBasesUpstream + "_upstream_" + numBasesDownstream + "_downstream";
		Sequence rtrn = new Sequence(rtrnID);
		rtrn.setSequenceBases(flankingSequence);
		return rtrn;
	}
	
	@Override
	public String toString() {
		return id + "_" + originalSequence.getId() + "_" + firstDeletedPos + "_" + lastDeletedPos;
	}
	
	
}
