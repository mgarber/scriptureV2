package nextgen.core.annotation;

import nextgen.core.general.TabbedReader;
import broad.core.error.ParseException;



/**
 * @author engreitz
 * Class that I started writing but didn't finish to represent a PSL annotation file (output by BLAT)
 * NOTE:  Checking in because it might be helpful down the line, but it is not currently used anywhere
 */
public class PSL extends BasicAnnotation {

	protected int match, mismatch, repMatch, Ncount, queryGap, queryGapBases, targetGap, targetGapBases, querySize, queryStart, queryEnd, targetSize;
	protected int[] queryStarts;
	
	public PSL(BasicAnnotation annotation) {
		super(annotation);
	}
	
	public void setMatch(int m) { match = m; }
	public void setMismatch(int m) { mismatch = m; }
	public void setRepMatch(int m) { repMatch = m; }
	public void setNcount(int c) { Ncount = c; }
	public void setQueryGap(int x) { queryGap = x; }
	public void setQueryGapBases(int x) { queryGapBases = x; }
	public void setTargetGap(int x) { targetGap = x; }
	public void setTargetGapBases(int x) { targetGapBases = x; }
	public void setQuerySize(int x) { querySize = x; }
	public void setQueryStart(int x) { queryStart = x; }
	public void setQueryEnd(int x) { queryEnd = x; }
	public void setTargetSize(int x) { targetSize = x; }
	
	public int getMatch() { return match; }

	
	public static class Factory implements TabbedReader.Factory<PSL> {
		@Override
		public PSL create(String[] rawFields) throws ParseException {
			if (rawFields.length < 21) {
				throw new IllegalArgumentException("Cannot create BasicAnnotation from less than 21 fields");
			}
			
			// Assemble the annotation from the block information
			int nBlocks = Integer.parseInt(rawFields[17]);
			
			String [] sizes = rawFields[18].split(",");
			String [] queryStarts = rawFields[19].split(",");
			String [] targetStarts = rawFields[20].split(",");
			if (sizes.length != nBlocks || queryStarts.length != nBlocks || targetStarts.length != nBlocks) {
				throw new ParseException("BAD BED FORMAT apparently the number of start ("+rawFields[19] +") and size ("+rawFields[18] +
						") items does not agree with the blockCount " + nBlocks);
			}
			
			PSL p = null;
			int[] qstarts = new int[nBlocks];
			String chr = rawFields[13];
			int start = Integer.parseInt(rawFields[15]);

			for (int i = 0; i < nBlocks; i++) {
				int blockSize = Integer.parseInt(sizes[i]);
				int queryStart = Integer.parseInt(queryStarts[i]);
				qstarts[i] = queryStart;
				int targetStart = Integer.parseInt(targetStarts[i]);
				
				if (i == 0) {
					p = new PSL(new BasicAnnotation(chr, start, start + blockSize));
				} else {
					p.addBlocks(new BasicAnnotation(chr, start + targetStart, start + targetStart + blockSize));
				}
			}
			if (p.getEnd() != Integer.parseInt(rawFields[16])) throw new IllegalArgumentException("End specified by blocks does not match BED end");

			p.setMatch(Integer.parseInt(rawFields[0]));
			p.setMismatch(Integer.parseInt(rawFields[1]));
			p.setRepMatch(Integer.parseInt(rawFields[2]));
			p.setNcount(Integer.parseInt(rawFields[3]));
			p.setQueryGap(Integer.parseInt(rawFields[4]));
			p.setQueryGapBases(Integer.parseInt(rawFields[5]));
			p.setTargetGap(Integer.parseInt(rawFields[6]));
			p.setTargetGapBases(Integer.parseInt(rawFields[7]));
			p.setOrientation(Annotation.Strand.fromString(rawFields[8]));
			p.setName(rawFields[9]);
			p.setQuerySize(Integer.parseInt(rawFields[10]));
			p.setQueryStart(Integer.parseInt(rawFields[11]));
			p.setQueryEnd(Integer.parseInt(rawFields[12]));
			p.setTargetSize(Integer.parseInt(rawFields[14]));
			
			return p;
		}
	}
}
