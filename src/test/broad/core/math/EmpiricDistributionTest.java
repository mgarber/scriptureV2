package broad.core.math;

import java.io.IOException;
import java.util.Random;

import broad.core.math.EmpiricalDistribution;

import junit.framework.TestCase;

public class EmpiricDistributionTest extends TestCase {
	public void testSummary() throws IOException {
		EmpiricalDistribution ed = new EmpiricalDistribution(5, 0, 4);
		Random r = new Random();
		for(int i = 0; i < 50; i++) {
			ed.add(r.nextInt(4)/(double)4);
		}
	}
	
	public void testGetNumberOfValuesLargerThan() throws IOException {
		EmpiricalDistribution ed = new EmpiricalDistribution(500, 0, 10);
		ed.add(1.0);
		
		assertEquals("There is one observation larger than 0.5",1,ed.getNumberOfValuesLargerThan(0.5));
		assertEquals("P(x < 0.5) = 0 since there is only one datapoint" ,0.0, ed.getCummulativeProbability(0.5));
		assertEquals("P(x < 2) = 1 since there is only one datapoint that is less than 2", 1.0, ed.getCummulativeProbability(2));
		
		ed.add(8.9);
		
		assertEquals("There is one observation larger than 0.5",2,ed.getNumberOfValuesLargerThan(0.5));
		assertEquals("There is one observation larger than 8",1,ed.getNumberOfValuesLargerThan(8));
		assertEquals("P(x < 2) = 1/2 since there is only one datapoint" ,0.5, ed.getCummulativeProbability(2));
		assertEquals("P(x < 9) = 1 since both datapoints are less 9", 1.0, ed.getCummulativeProbability(9));
		assertEquals("P(x < 8.87) = 1/2", 0.5, ed.getCummulativeProbability(8.87));
		
		ed.add(8.98);
		
		assertEquals("There is one observation larger than 0.5",3,ed.getNumberOfValuesLargerThan(0.5));
		assertEquals("There is one observation larger than 8",2,ed.getNumberOfValuesLargerThan(8));
		assertEquals("P(x < 2) = 1/2 since there is only one datapoint" ,1.0/3.0, ed.getCummulativeProbability(2));
		assertEquals("P(x < 9) = 1 since both datapoints are less 9", 1.0, ed.getCummulativeProbability(9));
		assertEquals("P(x < 8.93) = 22/3", 2.0/3.0, ed.getCummulativeProbability(8.93));
		
	}
}
