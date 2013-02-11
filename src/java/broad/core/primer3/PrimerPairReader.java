package broad.core.primer3;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.SequenceRegion;

public class PrimerPairReader {
	ArrayList<PrimerPair> primerPairs;

	public PrimerPairReader(String fileName) throws IOException {
		super();
		primerPairs = new ArrayList<PrimerPair>();
		File source = new File(fileName);
		BufferedReader br = new BufferedReader(new FileReader(source));
		String line;
		while((line = br.readLine()) != null) {
			if(line.startsWith("#")){
				continue;
			}
			String[] splitLine = line.split("\t");
			primerPairs.add(new PrimerPair(splitLine));
		}
		System.out.print("Closing "+fileName);
		br.close();
		System.out.print(" ..... Closed\n");	
	}
	
	public List<PrimerPair> getPrimerPairs() { return primerPairs; }
	
	public List<SequenceRegion> extractProductSequences(String sequenceFile) throws IOException {
		FastaSequenceIO fsio = new FastaSequenceIO(sequenceFile);		
		
		ArrayList<SequenceRegion> seqs = new ArrayList<SequenceRegion>();
		Iterator<PrimerPair> ppIt = primerPairs.iterator();
		while(ppIt.hasNext()) {
			PrimerPair pp = ppIt.next();
			SequenceRegion sr = new SequenceRegion(pp.getPrimerPairId());
			sr.setRegionStart(pp.getLeftPrimerPosition());
			sr.setRegionEnd(pp.getRightPrimerPosition());
			seqs.add(sr);
			System.out.println("region start " + sr.getRegionStart() + " region end " + sr.getRegionEnd());
		}
		
		fsio.extractRegions(seqs);
		
		return seqs;
	}
	
	public static void main(String [] args) throws IOException {
		if(args.length != 3) {
			System.err.println("Ussage: PrimerPairReader <primerPairFile> <sequence file> <output file>");
			return;
		}
		
		PrimerPairReader ppr = new PrimerPairReader(args[0]);
		List<SequenceRegion> extractedSeqs = ppr.extractProductSequences(args[1]);
		FastaSequenceIO fsio = new FastaSequenceIO(args[2]);
		fsio.write(extractedSeqs);
	}

}
