package xp.test.command;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.io.PrintWriter;

import xp.test.Converter.BamToPairedEndIterator;

import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.samtools.util.CloseableIterator;
import nextgen.core.alignment.Alignment;
import nextgen.core.alignment.PairedEndAlignment;

public class BamToPairedEndBed extends CommandLineProgram {
	
	private static final String PROGRAM_VERSION = "0.01";
	@Usage
    public String USAGE =  "java nextgen.core.programs.BamToPairedEndBed I=file.bam";

    @Option(doc="This is test paired end bam iterater and convert to simple bed format.", shortName=StandardOptionDefinitions.INPUT_SHORT_NAME) public File INPUT;
	@Option(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME) public File OUTPUT;
    public static void main(String[] argv){
    	System.exit(new BamToPairedEndBed().instanceMain(argv));
        
    }
    
    protected int doWork()
    {
    	
    	BamToPairedEndIterator iter = new BamToPairedEndIterator(INPUT);
		//SAMFileHeader header = reader.getFileHeader();
   		PrintStream out;
   		try
   		{
    	out = new PrintStream(OUTPUT);
    	
   		} catch (Exception e)
   		{
   			out=System.out;
   		}
    	int i=0;
		while (iter.hasNext()) {
		    Alignment  record= iter.next();
			//System.out.println("aha");
		    out.println(record);
		    i+=1;
		 
		   if (i%10000==0) {System.err.format("Processing %d reads      BufferSize:%d         \r",i,iter.getBufferSize());}
		}	
		iter.close();
        out.close();
        try{
        bedSort(OUTPUT);
        } catch (Exception e)
        {
        	System.err.println("Exception on bedSort "+OUTPUT.toString()+"\n"+e.getMessage());
        }
		return 0;
    }
    
    private static void bedSort(File OUTPUT) throws IOException, InterruptedException
    {
    	bedSort(OUTPUT.toString());
    }
    private static void bedSort(String OUTPUT)  throws IOException, InterruptedException
    {
    
    	String cmd="/Users/zhuxp/bin/x86_64/bedSort "+OUTPUT+" "+OUTPUT;
		Runtime run=Runtime.getRuntime();
		Process p=run.exec(cmd);
		p.waitFor();// what is this?
    	
    }
    
}

