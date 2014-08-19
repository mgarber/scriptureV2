package broad.core.primer3;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;

import broad.core.sequence.Sequence;

public class Primer3SequenceInputTags {
	private Sequence sequence;
	private ArrayList<SequenceRegionCoordinates> excludedRegions = new ArrayList<SequenceRegionCoordinates>();
	private ArrayList<SequenceRegionCoordinates> includedRegions = new ArrayList<SequenceRegionCoordinates>();
	private ArrayList<SequenceRegionCoordinates> targets = new ArrayList<SequenceRegionCoordinates>();
	private ArrayList<Integer> junctions=new ArrayList<Integer>();
	private String primerSequenceId;
	private String primerLeftInput;
	private String primerRightInput;
	private int    primerStartCodonPosition = -1000000;


	public Primer3SequenceInputTags(Sequence seq) {
		super();
		sequence = seq;
	}
	public static class SequenceRegionCoordinates {
		int start;
		int end;
		
		public SequenceRegionCoordinates(int start, int end) {
			this.start = start;
			this.end = end;			
		}
		
		public String toString() {
			return start + "," + length();
		}

		public int length() {
			return end - start;
		}
	}
	
	public int getPrimerStartCodonPosition() {
		return primerStartCodonPosition;
	}
	public void setPrimerStartCodonPosition(int primerStartCodonPosition) {
		this.primerStartCodonPosition = primerStartCodonPosition;
	}
	public List<SequenceRegionCoordinates> getExcludedRegions() {
		return excludedRegions;
	}
	public void addExcludedRegions(List<SequenceRegionCoordinates> excludedRegions) {
		this.excludedRegions.addAll(excludedRegions);
	}
	public void addExcludedRegion(SequenceRegionCoordinates excludedRegion) {
		this.excludedRegions.add(excludedRegion);
	}
	public List<SequenceRegionCoordinates> getIncludedRegions() {
		return includedRegions;
	}
	public void addIncludedRegions(List<SequenceRegionCoordinates> includedRegions) {
		this.includedRegions.addAll(includedRegions);
	}
	public void addIncludedRegion(SequenceRegionCoordinates includedRegion) {
		this.includedRegions.add(includedRegion);
	}
	public String getPrimerLeftInput() {
		return primerLeftInput;
	}
	public void setPrimerLeftInput(String primerLeftInput) {
		this.primerLeftInput = primerLeftInput;
	}
	public String getPrimerRightInput() {
		return primerRightInput;
	}
	public void setPrimerRightInput(String primerRightInput) {
		this.primerRightInput = primerRightInput;
	}
	public String getPrimerSequenceId() {
		return primerSequenceId == null ? sequence.getId() : primerSequenceId;
	}
	public void setPrimerSequenceId(String primerSequenceId) {
		this.primerSequenceId = primerSequenceId;
	}
	public List<SequenceRegionCoordinates> getTargets() {
		return targets;
	}
	public void addTargets(List<SequenceRegionCoordinates> targets) {
		this.targets.addAll(targets);
	}
	
	public void addTarget(SequenceRegionCoordinates coordinates) {
		targets.add(coordinates);
	}
	
	public Sequence getSequence() {
		return sequence;
	}
	public String regionListToString(List<SequenceRegionCoordinates> regions) {
		Iterator<SequenceRegionCoordinates> regionIt = regions.iterator();
		StringBuffer buf = new StringBuffer();
		while(regionIt.hasNext()) {
			buf.append(regionIt.next().toString());
			if(regionIt.hasNext()) {
				buf.append(" ");
			}
		}
		return buf.toString();
	}
	
	public String regionListToString(List<SequenceRegionCoordinates> regions, int maxNum) {
		Iterator<SequenceRegionCoordinates> regionIt = regions.iterator();
		StringBuffer buf = new StringBuffer();
		int counter=0;
		while(regionIt.hasNext() && counter<maxNum) {
			buf.append(regionIt.next().toString());
			if(regionIt.hasNext()) {
				buf.append(" ");
			}
			counter++;
		}
		return buf.toString();
	}
	public String getJunctions() {
		String str="";
		for(Integer junction: this.junctions){
			str+=""+(junction+1)+" ";
		}
		return str;
	}

	public void addJunctions(Collection<Integer> junctions){
		this.junctions.addAll(junctions);
	}

	
}
