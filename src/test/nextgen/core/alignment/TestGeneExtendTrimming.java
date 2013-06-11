package nextgen.core.alignment;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Map;


import broad.pda.annotation.BEDFileParser;
import nextgen.core.annotation.BasicAnnotation;
import nextgen.core.annotation.Gene;
import nextgen.core.annotation.Annotation.Strand;
import junit.framework.TestCase;

public class TestGeneExtendTrimming extends TestCase{
	
	public void testGeneExpand() throws IOException{
		
		Map<String,Collection<Gene>> data = BEDFileParser.loadDataByChr(new File("test.bed"));
		
		Collection<Gene> geS50 = new ArrayList<Gene>();
		Collection<Gene> geE50 = new ArrayList<Gene>();
		for(String chr:data.keySet()){
			//For each gene
			for(Gene gene:data.get(chr)){
				//Expand by 50 nts on start
				Gene dummy = gene.copy();
				dummy.expand(50, 0);
				geS50.add(dummy);
				
				//Expand by 50 nts on end
				dummy = gene.copy();
				dummy.expand(0, 50);
				geE50.add(dummy);
			}
		}
		BEDFileParser.writeFullBED("afterFixing.test.bed.expand.start.50.bed", geS50);
		BEDFileParser.writeFullBED("afterFixing.test.bed.expand.end.50.bed", geE50);
		
		//Test Annotations
		BasicAnnotation a = new BasicAnnotation("chr1",1000,2000,Strand.POSITIVE,"pos_candidate");
		BasicAnnotation b = new BasicAnnotation("chr1",1000,2000,Strand.NEGATIVE,"neg_candidate");
		System.out.println("Gene without change:");
		System.out.println(a.toBED());
		System.out.println(b.toBED());
		System.out.println(" Expanding 50nts start");
		a.expand(50, 0);
		b.expand(50, 0);
		System.out.println(a.toBED());
		System.out.println(b.toBED());
		System.out.println(" Expanding 50nts end");
		a.expand(0, 50);
		b.expand(0, 50);
		System.out.println(a.toBED());
		System.out.println(b.toBED());
		System.out.println("Trim 50nts start");
		a.trim(50, 0);
		b.trim(50,0);
		System.out.println(a.toBED());
		System.out.println(b.toBED());
		System.out.println("Trim 50nts end");
		a.trim(0, 50);
		b.trim(0, 50);
		System.out.println(a.toBED());
		System.out.println(b.toBED());
	}

}
