/**
 * 
 */
package nextgen.core.tests;

import java.util.ArrayList;
import java.util.Collection;

import junit.framework.TestCase;

import org.apache.log4j.Logger;

import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.Gene;
import nextgen.core.annotation.Annotation.Strand;

/**
 * @author skadri
 *
 */
public final class TestTrimming extends TestCase{
	
	private static Logger logger = Logger.getLogger(TestTrimming.class.getName());
	
	public void testTrimming() throws Exception {
		Gene exon1 = new Gene("chr1",100,200);
		Gene exon2 = new Gene("chr1",300,400); 
		Collection<Gene> exons = new ArrayList<Gene>();
		exons.add(exon1);
		exons.add(exon2);
		Gene gp = new Gene(exons);
		Gene gn = new Gene(exons);
		gn.setOrientation(Strand.NEGATIVE);
		gp.setOrientation(Strand.POSITIVE);
		logger.info("Testing positive trimming");
		logger.info(gp.toBED());
		logger.info("Trim() function - strand specific");
		Annotation d = gp.trim(10, 50);
		logger.info(d.toBED());
		
		//Test whether cast to Gene works
		Gene d_gene = new Gene(gp.trim(10, 50));
		logger.info(d_gene.toBED());
		logger.info("Testing negative trimming");
		logger.info(gn.toBED());
		logger.info("Trim() function - strand specific");
		d = gn.trim(10, 50);
		logger.info(d.toBED());
	}

}
