
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
public final class TestPositionOnReferenceCoordinate extends TestCase{
	
	private static Logger logger = Logger.getLogger(TestPositionOnReferenceCoordinate.class.getName());
	
	public void testposition() throws Exception {
		Gene exon1 = new Gene("chr1",100,200);
		Gene exon2 = new Gene("chr1",300,400); 
		Collection<Gene> exons = new ArrayList<Gene>();
		exons.add(exon1);
		exons.add(exon2);
		Gene gp = new Gene(exons);
		Gene gn = new Gene(exons);
		gn.setOrientation(Strand.NEGATIVE);
		gp.setOrientation(Strand.POSITIVE);
		logger.info("Testing positive gene");
		logger.info(gp.toBED());
		logger.info("What is position on reference of 350");
		logger.info(gp.getPositionAtReferenceCoordinate(350));
		logger.info("What is position on reference of 100");
		logger.info(gp.getPositionAtReferenceCoordinate(100));
		logger.info("What is position on reference of 399");
		logger.info(gp.getPositionAtReferenceCoordinate(399));
		//Test whether cast to Gene works
		logger.info("Testing negative gene");
		logger.info(gn.toBED());
		logger.info("What is position on reference of 350");
		logger.info(gn.getPositionAtReferenceCoordinate(350));
		logger.info("What is position on reference of 100");
		logger.info(gn.getPositionAtReferenceCoordinate(100));
		logger.info("What is position on reference of 399");
		logger.info(gn.getPositionAtReferenceCoordinate(399));
	}

}
