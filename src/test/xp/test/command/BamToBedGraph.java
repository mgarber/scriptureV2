
package xp.test.command;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.PrintStream;
import java.util.Iterator;

import org.apache.log4j.Logger;

import xp.test.Basic.BedGraphMultiScore;
import xp.test.Converter.BamIteratorFactory;
import xp.test.Converter.BedGraphMultiScoreReader;
import xp.test.Utils.JieCodeSortingCollection;

import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.Usage;
import nextgen.core.annotation.Annotation;
import nextgen.core.coordinatesystem.GenomicSpace;

/**
 *  Created on 2013-6-26  
 */



public class BamToBedGraph extends CommandLineProgram {

private static final String PROGRAM_VERSION = "0.01";
@Usage
public static final String USAGE = "Usage: ";
@Option(doc = "input file", shortName = "i")
public File INPUT;
@Option(doc = "output file", shortName = "o")
public String OUT = "stdout";
@Option(doc="Chromosome Size",shortName="G") public String CHROMSIZES;

public static void main(String[] argv) {
	System.exit(new BamToBedGraph().instanceMain(argv));
}

@Override
protected int doWork() {
	
	GenomicSpace GENOMESPACE = new GenomicSpace(CHROMSIZES);
	JieCodeSortingCollection codes=new JieCodeSortingCollection(GENOMESPACE);
	Iterator<? extends Annotation> iter = BamIteratorFactory.makeIterator(INPUT); 
	codes.add(iter,0,false);
	//codes.add(iter,i,ignoreBlocksDetails[i]);
	// TODO Auto-generated method stub
	double globalLambda=(double)codes.getCoverage(0)/GENOMESPACE.getLength();
	PrintStream out = null;
	if(OUT.equalsIgnoreCase("stdout"))
	{
		out=System.out;
	}
	else
	{
		try {
			out=new PrintStream(new File(OUT));
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	BedGraphMultiScoreReader reader=new BedGraphMultiScoreReader(codes);
    
	while(reader.hasNext())
    {
    	BedGraphMultiScore b = reader.next();
    	if(b.getScore()[0]==0) continue;
    	out.println(String.format("%s\t%d\t%d\t%d",codes.getTid2chr().get(b.getTid()),b.getStart(),b.getEnd(),b.getScore()[0]));
    }
    out.close();
	return 0;
}

}
