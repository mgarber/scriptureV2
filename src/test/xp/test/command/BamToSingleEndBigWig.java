package xp.test.command;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.io.PrintWriter;

import xp.test.Converter.BamToSingleEndIterator;

import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.samtools.util.CloseableIterator;
import nextgen.core.alignment.Alignment;
import nextgen.core.alignment.SingleEndAlignment;

public class BamToSingleEndBigWig extends CommandLineProgram {
	
	private static final String PROGRAM_VERSION = "0.01";
	@Usage
    public String USAGE =  "java nextgen.core.programs.BamToSingleEndBed I=file.bam";
	
	private static final String SingleEndBedSuffix = ".SingleEnd.bed";
	private static final String SingleEndWigSuffix = ".SingleEnd.bedGraph";
	private static final String SingleEndBigWigSuffix = ".SingleEnd.bw";

    @Option(doc="This program iterate bam file as single end and convert it into bigwig format.", shortName=StandardOptionDefinitions.INPUT_SHORT_NAME) public File INPUT;
    public static void main(String[] argv){
    	System.exit(new BamToSingleEndBigWig().instanceMain(argv));
        
    }
    
    protected int doWork()
    {
    	
    	BamToSingleEndIterator iter = new BamToSingleEndIterator(INPUT);

    	//SAMFileHeader header = reader.getFileHeader();
    	String filename = INPUT.getAbsolutePath();
    	//String dirname=INPUT.getAbsolutePath();
    	
    	String name = filename.substring(0,filename.lastIndexOf("."));
    	
    	String sBedName = name+SingleEndBedSuffix;
    	System.err.println(sBedName);
    	String sWigName = name+SingleEndWigSuffix;
    	String sBigWigName = name+SingleEndBigWigSuffix;
    	System.err.println(sWigName);
    	System.err.println(sBigWigName);
    	
    	PrintStream fSingleEndBed;
    	try
    	{
    		fSingleEndBed=new PrintStream(sBedName);
    	
    	int i=0;
		while (iter.hasNext()) {
		    Alignment  record= iter.next();
		    fSingleEndBed.println(record.toShortBED());
		    i+=1;
		 
		   if (i%10000==0) {System.err.format("Processing %d reads      \r",i);}
		}	
		iter.close();
        fSingleEndBed.close();
    	}
    	catch (Exception e)
    	{
    		System.err.print("can't open "+ sBedName +"\t"+ e);
    	}
        
        
        try{
          bedSort(sBedName);
        } catch (Exception e)
        {
        	System.err.println("Exception on bedSort "+sBedName+"\n"+e.getMessage());
        }
        try{
        	
        	bedToWig(sBedName,sWigName,sBigWigName);
        }
        catch (Exception e)
        {
        	System.err.println("Exception on Convert bigwig "+sBedName+"\n"+e.getMessage());
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
    private static void bedToWig(String bedFileName,String wigFileName,String bigWigFileName) throws IOException, InterruptedException
    {
    /*
     * genomeCoverageBed -split -bg -i accepted_hits.sorted.bed -g dm3.chrom.sizes > accepted_hits.bedgraph	
     * 
     * wigToBigWig accepted_hits.bedgraph dm3.chrom.sizes myfile.bw
     */
     String home="/Users/zhuxp"	;
     String genomeSize=home+"/Data/genome/mm10.chrom.sizes";
     String genomeCoverageBed = home+"/bin/genomeCoverageBed";
     String cmd=genomeCoverageBed + " -split -bg -i "+bedFileName +" -g "+genomeSize+ " > "+ wigFileName;
    // String cmd=genomeCoverageBed + " -split -bg -i "+bedFileName +" -g "+genomeSize;
     String[] cmds= {"sh","-c",cmd};
     Runtime run1=Runtime.getRuntime();
     System.err.println(cmd);
	 Process p1=run1.exec(cmds);
	 p1.waitFor();
	 
	 String wigToBigWig = home+"/bin/x86_64/wigToBigWig";
	 String cmd2=wigToBigWig + " " + wigFileName + " " +genomeSize+" "+ bigWigFileName;
	 Runtime run2=Runtime.getRuntime();
     System.err.println(cmd2);
	 Process p2=run2.exec(cmd2);
	 p2.waitFor();
	 File wigFile = new File(wigFileName);
	 wigFile.delete();
	 
	 
    		 
    
    }
    
}
