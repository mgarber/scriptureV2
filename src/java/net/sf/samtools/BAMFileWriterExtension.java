package net.sf.samtools;

import java.io.File;


/**
 * This wraper is to enable setting the temporary directory.
 * The visibility is default so it can't be accessed by client classes
 * @author mgarber
 *
 */
public class BAMFileWriterExtension extends BAMFileWriter {

	public BAMFileWriterExtension(File path) {
		super(path);
		// TODO Auto-generated constructor stub
	}
	
	public void setTemporaryDirectory(File tmpDirPath) {
		super.setTempDirectory(tmpDirPath);
	}
	
}
