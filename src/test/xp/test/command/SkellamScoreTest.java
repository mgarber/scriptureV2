package xp.test.command;

import java.io.*;
import java.util.Iterator;

import org.apache.log4j.Logger;

import xp.test.DBI.ShortBEDTabixDBI;
import xp.test.ScoreModel.SkellamScoreSystem;

import broad.core.annotation.ShortBED;

import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.Usage;
//import net.sf.samtools.BAMFileWriter;
//import net.sf.samtools.SAMFileHeader;
import nextgen.core.coordinatesystem.GenomicSpace;
import nextgen.core.feature.Window;
//import nextgen.core.model.AlignmentModel;

public class  SkellamScoreTest extends CommandLineProgram {
	static Logger logger = Logger.getLogger(WindowTabixCount.class.getName());	
	private static final String PROGRAM_VERSION = "0.01";
	@Usage

    @Option(doc="This is caseA tabix gzip", shortName="A") public File CASEA;
    @Option(doc="This is caseB tabix gzip", shortName="B") public File CASEB;
    @Option(doc="This is control A tabix gzip", shortName="X") public File CONTROLA;
    @Option(doc="This is control B tabix gzip", shortName="Y") public File CONTROLB;
	@Option(doc="Bin Size",shortName="w") public int WINDOWSIZE;
	@Option(doc="output file",shortName="o") public String FOUT;
	
    public static void main(String[] argv){
    	System.exit(new SkellamScoreTest().instanceMain(argv));
        
    }
    
    protected int doWork() 
    {
		
    	
    	/*
    	 * Construct chromosome space
    	 */
    	 int NEARBY=50000;
    	 logger.info("initialize genomic space");
    	 GenomicSpace GENOMESPACE = new GenomicSpace("/Users/zhuxp/Data/genome/mm9.chromSizes");
    	 logger.info("initialize genomic space done");
    	
    	 logger.info(GENOMESPACE.getReferenceNames());
    	 try {
    	 ShortBEDTabixDBI caseADBI= new ShortBEDTabixDBI(CASEA.getAbsolutePath(),GENOMESPACE);
    	 ShortBEDTabixDBI caseBDBI= new ShortBEDTabixDBI(CASEB.getAbsolutePath(),GENOMESPACE);
    	 ShortBEDTabixDBI controlADBI= new ShortBEDTabixDBI(CONTROLA.getAbsolutePath(),GENOMESPACE);
    	 ShortBEDTabixDBI controlBDBI= new ShortBEDTabixDBI(CONTROLB.getAbsolutePath(),GENOMESPACE);
    	SkellamScoreSystem scoreSystem = new SkellamScoreSystem(caseADBI,caseBDBI,controlADBI,controlBDBI,NEARBY);
    	/* init ChIPSeqPeakCaller
    	 * 
    	 */
    	 FileWriter fout;
		/*
    	Iterator<ShortBED> tmpout=controlBDBI.query("chr1",1,10000000);
    	
    	while(tmpout.hasNext())
    	{
    		ShortBED a=tmpout.next();
    		System.err.println(a);
    	}
    	*/
    	fout = new FileWriter(FOUT);
		
    	long GENOMESIZE=GENOMESPACE.getLength();
    	Iterator<Window> iter=GENOMESPACE.getWindowIterator(WINDOWSIZE, 0);
		
       	fout.write("# WINDOW SIZE   :" +String.format("%d", WINDOWSIZE)+"\n");
       	fout.write("# GENOME SIZE   :" + String.format("%d",GENOMESIZE)+"\n");
       	
       	int i=0;
       
        while(iter.hasNext())
        {
        	i++;
        	if (i%1000==0) {logger.info("processing " + i +" windows"); fout.flush();}
        	Window window=iter.next();
        	Double count=scoreSystem.getScore(window.getChr(),window.getStart(),window.getEnd());
        	if (count!=0.0)
        	{	
			fout.write(window.getChr()+"\t"+window.getStart()+"\t"+window.getEnd()+"\t"+String.format("%f",count)+"\n");
        	}
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