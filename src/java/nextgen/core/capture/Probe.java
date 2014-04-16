package nextgen.core.capture;

import nextgen.core.capture.arrayscheme.ProbeLayout;
import broad.core.sequence.Sequence;

/**
 * @author prussell
 * The probe itself without primers
 */
public class Probe implements Comparable<Probe> {
	
	private Sequence parent = null;
	private String parentSequenceId;
	private String id;
	private int start;
	private int size;
	private boolean antisense;
	private ProbeLayout probeLayout = null;
	private String probeLayoutString = null;
	
	
	/**
	 * @param parentTranscript Parent transcript
	 * @param layout Probe layout that produced the probe
	 * @param probeID Unique probe identifier
	 * @param startPosOnTranscript First transcript position included in probe
	 * @param probeSize Probe length
	 * @param antisenseToTranscript Whether the probe is antisense to transcript
	 */
	public Probe(Sequence parentTranscript, ProbeLayout layout, String probeID, int startPosOnTranscript, int probeSize, boolean antisenseToTranscript) {
		parent = parentTranscript;
		id = probeID;
		start = startPosOnTranscript;
		size = probeSize;
		antisense = antisenseToTranscript;
		probeLayout = layout;
	}
	
	public Probe(String parentSequenceId, String layoutString, String probeID, int startPosOnTranscript, int probeSize, boolean antisenseToTranscript) {
		this.parentSequenceId = parentSequenceId;
		id = probeID;
		start = startPosOnTranscript;
		size = probeSize;
		antisense = antisenseToTranscript;
		probeLayoutString = layoutString;
	}
	
	
	/**
	 * @return Probe ID
	 */
	public String getID() {
		return id;
	}
	
	/**
	 * @return Parent transcript
	 */
	public Sequence getParentTranscript() {
		return parent;
	}
	
	public String getParentTranscriptId() {
		return (parent == null) ? parentSequenceId : parent.getId();
	}
	
	/**
	 * @return Probe start position on transcript
	 */
	public int getStartPosOnTranscript() {
		return start;
	}
	
	/**
	 * @return First transcript position after probe
	 */
	public int getEndPosOnTranscript() {
		return start + size;
	}
	
	/**
	 * @return Probe size
	 */
	public int getSize() {
		return size;
	}
	
	/**
	 * @return The probe layout that created the probe
	 */
	public ProbeLayout getProbeLayout() {
		return probeLayout;
	}
	
	public String getProbeLayoutString() {
		if (probeLayout == null) return probeLayoutString;
		else return probeLayout.toString();
	}
	
	/**
	 * @return Probe tiling path
	 */
	public int getTilingPath() {
		return start % size;
	}
	
	/**
	 * @return Whether the probe sequence is antisense to the transcript
	 */
	public boolean isAntisenseToTranscript() {
		return antisense;
	}
	
	/**
	 * @return The sequence of the probe
	 */
	public String getProbeSequence() {
		Sequence probeSeq = parent.getSubSequence("", start, start + size);
		if(antisense) {
			probeSeq.reverse();
		}
		return probeSeq.getSequenceBases();
	}
	
	@Override
	public boolean equals(Object o) {
		if(!o.getClass().equals(getClass())) {
			return false;
		}
		Probe p = (Probe)o;
		return p.hashCode() == hashCode();
	}
	
	@Override
	public int hashCode() {
		String str = parent.getId() + "_" + id + "_" + start + "_" + size + "_" + Boolean.valueOf(antisense).toString() + "_" + probeLayout.toString();
		return str.hashCode();
	}
	
	@Override
	public int compareTo(Probe other) {
		if(!parent.equals(other.getParentTranscript())) {
			return parent.getId().compareTo(other.getID());
		}
		if(start == other.getStartPosOnTranscript()) {
			if(size == other.getSize()) {
				if(antisense && !other.isAntisenseToTranscript()) return -1;
				if(!antisense && other.isAntisenseToTranscript()) return 1;
				return probeLayout.toString().compareTo(other.getProbeLayout().toString());
			}
			return size - other.getSize();
		}
		return start - other.getStartPosOnTranscript();
	}
	
}
