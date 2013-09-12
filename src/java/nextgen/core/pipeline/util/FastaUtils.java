package nextgen.core.pipeline.util;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;


import nextgen.core.job.LSFJob;

import org.apache.log4j.Logger;

import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;

/**
 * @author prussell
 *
 */
public class FastaUtils {

	private static Logger logger = Logger.getLogger(FastaUtils.class.getName());
	
	/**
	 * Index a fasta file using samtools faidx and write bsub output file to working directory
	 * @param fastaFileName The fasta file to index
	 * @param samtoolsExecutable Path to samtools executable
	 * @param scheduler Name of scheduler e.g. "LSF" or "SGE"
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public static void indexFastaFile(String fastaFileName, String samtoolsExecutable, String scheduler) throws IOException, InterruptedException {
		indexFastaFile(fastaFileName, samtoolsExecutable, ".", scheduler);
	}
	
	/**
	 * Index a fasta file using samtools faidx
	 * @param fastaFileName The fasta file to index
	 * @param samtoolsExecutable Path to samtools executable
	 * @param bsubOutDir Directory to write bsub output to
	 * @param scheduler Name of scheduler e.g. "LSF" or "SGE"
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public static void indexFastaFile(String fastaFileName, String samtoolsExecutable, String bsubOutDir, String scheduler) throws IOException, InterruptedException {
		String cmmd = samtoolsExecutable + " faidx " + fastaFileName;
		logger.info("Running samtools command: " + cmmd);
		if(scheduler.equals("LSF")) {
		String jobID = Long.valueOf(System.currentTimeMillis()).toString();
		logger.info("LSF job ID is " + jobID + ".");
		LSFJob job = new LSFJob(Runtime.getRuntime(), jobID, cmmd, bsubOutDir + "/index_fasta_" + jobID + ".bsub", "hour", 4);
		job.submit();
		logger.info("Waiting for samtools faidx job to finish...");
		job.waitFor();
		} else {
			throw new IllegalArgumentException("Scheduler " + scheduler + " not supported.");
		}
	}
	
	/**
	 * Write size file for a fasta file
	 * @param refFasta The fasta file
	 * @return The name of the size file that has been written
	 * @throws IOException
	 */
	public static String writeSizeFile(String refFasta) throws IOException {
		String sizeFileName = refFasta + ".sizes";
		File chrSizeFile = new File(sizeFileName);
		if(!chrSizeFile.exists()) {
			FileWriter w = new FileWriter(sizeFileName);
			logger.info("Writing chromosome sizes to file " + sizeFileName);
			FastaSequenceIO fsio = new FastaSequenceIO(refFasta);
			Collection<Sequence> seqs = fsio.loadAll();
			for(Sequence seq : seqs) {
				w.write(seq.getId() + "\t" + seq.getLength() + "\n");
			}
			w.close();
		}
		return sizeFileName;
	}


}
