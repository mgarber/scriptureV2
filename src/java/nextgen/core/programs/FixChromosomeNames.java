package nextgen.core.programs;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import net.sf.picard.util.Log;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMProgramRecord;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;

/**
 * Some pipelines align to genome builds where chromosomes names do not 
 * have the "chr" prefix 
 * @author mgarber
 *
 */
public class FixChromosomeNames {
	private static final Log log = Log.getInstance(FixChromosomeNames.class);
	static final String USAGE = "Convert the sequence names of a bam file by adding the chr prefix: " + 
			"\n\t-in <Path to the  BAM or SAM alignment file (NO standard input is supported must input a file path)>" +
			"\n\t-seqsToFix <Comma sparated sequences to fix, all others will be left untouched>" +
			"\n\t-out <Ouput BAM/SAM file or standard out if none is provided>" +
			"\n";
	
	

	public static void main (String [] args) throws IOException {
		ArgumentMap argMap = CLUtil.getParameters(args,USAGE , "default");
		String alnFile = argMap.getInput();
		String out = argMap.getOutput();
		List<String> seqsToChange = new ArrayList<String>();
		if(argMap.containsKey("seqsToFix")) {
			seqsToChange = CLUtil.listFromArray(argMap.getMandatory("seqsToFix").split(","));
		}

		SAMFileReader reader = new SAMFileReader(new File(alnFile));
		SAMFileWriterFactory samfwf = new SAMFileWriterFactory();
		SAMFileHeader oldHeader = reader.getFileHeader();
		SAMFileHeader header = oldHeader.clone();
		header.setSortOrder(oldHeader.getSortOrder());
		oldHeader.getAttributes();
		for (Map.Entry<String, String> he : oldHeader.getAttributes()) { 
			header.setAttribute(he.getKey(), he.getValue());
		}
		
		
		
		header.addProgramRecord(new SAMProgramRecord("nextgen.core.programs.FixChromosomeNames -in " + alnFile + " -out " + out));
		SAMSequenceDictionary ssd = header.getSequenceDictionary();
		SAMSequenceDictionary newSSD = new SAMSequenceDictionary();
		List<SAMSequenceRecord > sequences = ssd.getSequences();
		for (SAMSequenceRecord s : sequences) {
			if(seqsToChange.isEmpty() || seqsToChange.contains(s.getSequenceName())) {
				SAMSequenceRecord ns = new SAMSequenceRecord("chr" + s.getSequenceName(), s.getSequenceLength());
				ns.setAssembly(s.getAssembly());
				ns.setSpecies(ns.getSpecies());
				for (Map.Entry<String, String> e : s.getAttributes()) { 
					ns.setAttribute(e.getKey(), e.getValue());
				}
				newSSD.addSequence(ns);
			} else {
				newSSD.addSequence(s);
			}
		}
		header.setSequenceDictionary(newSSD);
		
		samfwf.setCreateIndex(true);
		SAMFileWriter writer = samfwf.makeBAMWriter(header, true, new File(out))  ;
		
		Iterator<SAMRecord> alnIt = reader.iterator();
		while (alnIt.hasNext()) {
			SAMRecord a = alnIt.next();
			if(seqsToChange.isEmpty() || seqsToChange.contains(a.getReferenceName())) {
				a.setReferenceName("chr"+a.getReferenceName());
				a.setHeader(header);
			}
			writer.addAlignment(a);
		}
		
		reader.close();
		writer.close();

	}
}

	



