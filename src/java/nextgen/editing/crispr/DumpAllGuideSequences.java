package nextgen.editing.crispr;

import java.io.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.io.IoUtil;
import net.sf.picard.util.Log;
import broad.core.motif.SequenceMotif;
import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;
import broad.core.sequence.SequenceRegion;
import broad.pda.seq.rap.CRISPRDesigner;

/**
 * @author engreitz
 * This class dumps all guide sequences in a given FASTA file to a file.  Build
 * to be completely lightweight, since I'm using this to dump all guides in the
 * entire genome.
 */
public class DumpAllGuideSequences extends CommandLineProgram {
    private static final Log log = Log.getInstance(DumpAllGuideSequences.class);

	@Option(doc="Genome fasta file")
	public File GENOME_FASTA = new File("/seq/lincRNA/data/mm9.nonrandom.fa");
	
	@Option(doc="Bed file containing subsets of genome to output", optional=true)
	public File REGIONS = null;
	
	@Option(doc="Output file with 20mer guide sequences, one per line")
	public File OUTPUT;
	
	/**
	 * Stock main method.
	 *
	 * @param args main arguments
	 */
	public static void main(final String[] args) {
		System.exit(new DumpAllGuideSequences().instanceMain(args));
	}
	
	
    @Override
    protected int doWork() {
        IoUtil.assertFileIsReadable(GENOME_FASTA);
        IoUtil.assertFileIsWritable(OUTPUT);
        
        try {
        	BufferedWriter bw = new BufferedWriter(new FileWriter(OUTPUT));
			Map<String,Sequence> chrs = FastaSequenceIO.getChrSequencesFromFasta(GENOME_FASTA.getAbsolutePath());
			
			if (REGIONS != null) {
				// Not implemented yet
			} else {
				for (String chr : chrs.keySet()) {
					log.info("Starting " + chr);
					Sequence seq = chrs.get(chr);
					outputAll(seq, bw);
				}
			}
			bw.close();
        } catch (Exception e) {
        	log.error(e);
        }

        return 0;
    }
    
    
    private void outputAll(Sequence seq, BufferedWriter bw) throws IOException {
    	seq.setSequenceBases(seq.getSequenceBases().toUpperCase());
		Matcher m = Pattern.compile(".{21}GG").matcher(seq.getSequenceBases());
		for (int i = 0; i < seq.getLength() - 23; i++) {
			m.region(i,i+23);
			if (m.find()) {
				bw.write(seq.getId() + "\t" + i + "\t" + (i+23) + "\t" + m.group() + "\t0\t+");
				bw.newLine();
			}
		}

		m = Pattern.compile("CC.{21}").matcher(seq.getSequenceBases());
		for (int i = 0; i < seq.getLength() - 23; i++) {
			m.region(i,i+23);
			if (m.find()) {
				bw.write(seq.getId() + "\t" + i + "\t" + (i+23) + "\t" + m.group() + "\t0\t-");
				bw.newLine();
			}
		}
    }
}
