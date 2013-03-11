package xp.test.command;

import java.io.*;
import java.util.Iterator;

import org.apache.log4j.Logger;

import xp.test.DBI.ShortBEDTabixDBI;

import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.Usage;
//import net.sf.samtools.BAMFileWriter;
//import net.sf.samtools.SAMFileHeader;
import nextgen.core.coordinatesystem.GenomicSpace;
import nextgen.core.feature.Window;
//import nextgen.core.model.AlignmentModel;

public class  WindowTabixCount extends CommandLineProgram {
	static Logger logger = Logger.getLogger(WindowTabixCount.class.getName());	
	private static final String PROGRAM_VERSION = "0.01";
	@Usage

    @Option(doc="This is signal bam", shortName="s") public File SIGNAL;
	@Option(doc="Bin Size",shortName="w") public int WINDOWSIZE;
	@Option(doc="output file",shortName="o") public String FOUT;
	
    public static void main(String[] argv){
    	System.exit(new WindowTabixCount().instanceMain(argv));
        
    }
    
    protected int doWork() 
    {
		
    	
    	/*
    	 * Construct chromosome space
    	 */
    	 logger.info("initialize genomic space");
    	 GenomicSpace GENOMESPACE = new GenomicSpace("/Users/zhuxp/Data/genome/mm9.chromSizes");
    	 logger.info("initialize genomic space done");
    	
    	 logger.info(GENOMESPACE.getReferenceNames());
    	
    	/* init ChIPSeqPeakCaller
    	 * 
    	 */
    	try {
    	 FileWriter fout;
		
			
    	
    	fout = new FileWriter(FOUT);
		
    	long GENOMESIZE=GENOMESPACE.getLength();
    	Iterator<Window> iter=GENOMESPACE.getWindowIterator(WINDOWSIZE, 0);
    	ShortBEDTabixDBI signalData=new ShortBEDTabixDBI(SIGNAL.getAbsolutePath(),GENOMESPACE);
       	logger.info("counting global stat.");
        double globalCount=signalData.getGlobalCount();
        
       	
	    fout.write("# INPUT TABIX FILE:  "+ SIGNAL.toString()+"\n");
		
       	fout.write("# GLOBAL COUNT  :" + String.format("%f", globalCount)+"\n");
       	fout.write("# WINDOW SIZE   :" +String.format("%d", WINDOWSIZE)+"\n");
       	fout.write("# GENOME SIZE   :" + String.format("%d",GENOMESIZE)+"\n");
       	fout.write("# AVERAGE       :" + String.format("%f", (float)globalCount*WINDOWSIZE/GENOMESIZE)+"\n");
       	
       	int i=0;
        while(iter.hasNext())
        {
        	i++;
        	if (i%100000==0) {logger.info("processing " + i +" windows"); fout.flush();}
        	Window window=iter.next();
        	double count=signalData.getCount(window.getChr(),window.getStart(),window.getEnd());
			fout.write(String.format("%f",count)+"\n");
        }
        fout.close();
    	}
    	catch (IOException e)
    	{
    	//To do	
    	}
    	return 0;
    }
}