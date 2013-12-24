/**
 * 
 */
package broad.pda.seq.clip;

import java.io.File;
import java.io.IOException;
import java.util.Map;

import org.apache.log4j.Level;
import org.ggf.drmaa.DrmaaException;

/**
 * @author shari
 *
 */
public class BatchedMultiSamplePeakCaller extends MultiSamplePeakCaller {

	private String sampleName;
	private String chr;
	private SampleData sampleData;
	
	private BatchedMultiSamplePeakCaller(MultiSamplePeakCaller m, String sample, String chrName) throws IOException, DrmaaException {
		super(m);
		sampleName = sample;
		chr = chrName;
				
		boolean foundSample = false;
		for(SampleData s : allSamples) {
			if(s.getSampleName().equals(sampleName)) {
				if(foundSample) {
					throw new IllegalArgumentException("Found sample name " + sampleName + " twice.");
				}
				sampleData = s;
				foundSample = true;
			}
		}
		if(foundSample) {
			logger.info("Found sample " + sampleName + ".");
		}
		
	}

	public static String[] extendSuperArgsForSampleAndChr(String[] commandArgs, String sampleName, String chrName) {
		String[] rtrn = new String[commandArgs.length + 2];
		for(int i=0; i < commandArgs.length; i++) {
			rtrn[i] = commandArgs[i];
		}
		rtrn[commandArgs.length] = sampleName;
		rtrn[commandArgs.length + 1] = chrName;
		return rtrn;
	}
	
	private static String[] getSuperCommandArgs(String[] extendedCommandArgs) {
		String[] rtrn = new String[extendedCommandArgs.length - 2];
		for(int i=0; i < rtrn.length; i++) {
			rtrn[i] = extendedCommandArgs[i];
		}
		return rtrn;
	}
	
	private static String getSampleName(String[] extendedCommandArgs) {
		return extendedCommandArgs[extendedCommandArgs.length - 2];
	}
	
	private static String getChrName(String[] extendedCommandArgs) {
		return extendedCommandArgs[extendedCommandArgs.length - 1];
	}
	
	private void writePeaks(String outDir) throws IOException {
		Map<String,String> outFiles = getPeakBedFileName(sampleData, outDir, chr);
		writeSingleSamplePeaks(sampleData, outFiles, chr);
	}
	
	/**
	 * @param args
	 * @throws IOException 
	 * @throws DrmaaException 
	 */
	public static void main(String[] args) throws IOException, DrmaaException {
		
		
		String[] superArgs = getSuperCommandArgs(args);
		
		String sampleName = getSampleName(args);
		String chrName = getChrName(args);
		MultiSamplePeakCaller m = createFromCommandArgs(superArgs);
		BatchedMultiSamplePeakCaller b = new BatchedMultiSamplePeakCaller(m, sampleName, chrName);
		String outDir = commandLineOutDir(superArgs);
		File o = new File(outDir);
		@SuppressWarnings("unused")
		boolean madeDir = o.mkdir();
		if(!o.exists()) {
			throw new IOException("Could not create directory " + outDir);
		}
		
		b.initializeFilterRejectWriters(chrName, outDir + "/" + FILTER_REJECT_DIR);
		if(commandLineHasDebugFlag(superArgs)) {
			b.setLoggerLevel(Level.DEBUG);
		}
		b.writePeaks(outDir);
		b.closeFilterRejectWriters();
	}

}
