package nextgen.core.pipeline.util;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;


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
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public static void indexFastaFile(String fastaFileName, String samtoolsExecutable) throws IOException, InterruptedException {
		indexFastaFile(fastaFileName, samtoolsExecutable, ".");
	}
	
	/**
	 * Index a fasta file using samtools faidx
	 * @param fastaFileName The fasta file to index
	 * @param samtoolsExecutable Path to samtools executable
	 * @param bsubOutDir Directory to write bsub output to
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public static void indexFastaFile(String fastaFileName, String samtoolsExecutable, String bsubOutDir) throws IOException, InterruptedException {
		String cmmd = samtoolsExecutable + " faidx " + fastaFileName;
		logger.info("Running samtools command: " + cmmd);
		String jobID = Long.valueOf(System.currentTimeMillis()).toString();
		logger.info("LSF job ID is " + jobID + ".");
		LSFUtils.bsubProcess(Runtime.getRuntime(), jobID, cmmd, bsubOutDir + "/index_fasta_" + jobID + ".bsub", "hour", 4);
		logger.info("Waiting for samtools faidx job to finish...");
		LSFUtils.waitForJobs(jobID, Runtime.getRuntime());
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
