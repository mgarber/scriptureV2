package Jama;

import broad.core.math.Statistics;
import junit.framework.TestCase;


public class QuantileNormalizationTest extends TestCase{
	
	public void testSimple() {
		Matrix m = new Matrix(5,3);
		m.set(0,0,2);
		m.set(1,0,1.1);
		m.set(2,0,5.5);
		m.set(3,0,1);
		m.set(4,0,5);
		
		m.set(0,1,3);
		m.set(1,1,5);
		m.set(2,1,7);
		m.set(3,1,6);
		m.set(4,1,6.3);
		
		m.set(0,2,2);
		m.set(1,2,3);
		m.set(2,2,9);
		m.set(3,2,6);
		m.set(4,2,2.5);

		Matrix msorted = new Matrix(5,3);
		
		msorted.set(0,0,1);
		msorted.set(1,0,1.1);
		msorted.set(2,0,2);
		msorted.set(3,0,5);
		msorted.set(4,0,5.5);
		
		msorted.set(0,1,3);
		msorted.set(1,1,5);
		msorted.set(2,1,6);
		msorted.set(3,1,6.3);
		msorted.set(4,1,7);
		
		msorted.set(0,2,2);
		msorted.set(1,2,2.5);
		msorted.set(2,2,3);
		msorted.set(3,2,6);
		msorted.set(4,2,9);
		
		double [] sortedRowMeans = new double[5];
		
		for(int i = 0; i < 5; i++) {	
			sortedRowMeans[i] = Statistics.mean(msorted.getRow(i));
		}
		
		Matrix mqnormalized = new Matrix(5,3);
		
		mqnormalized.set(0,0,sortedRowMeans[2]);
		mqnormalized.set(1,0,sortedRowMeans[1]);
		mqnormalized.set(2,0,sortedRowMeans[4]);
		mqnormalized.set(3,0,sortedRowMeans[0]);
		mqnormalized.set(4,0,sortedRowMeans[3]);
		
		mqnormalized.set(0,1,sortedRowMeans[0]);
		mqnormalized.set(1,1,sortedRowMeans[1]);
		mqnormalized.set(2,1,sortedRowMeans[4]);
		mqnormalized.set(3,1,sortedRowMeans[2]);
		mqnormalized.set(4,1,sortedRowMeans[3]);
		
		mqnormalized.set(0,2,sortedRowMeans[0]);
		mqnormalized.set(1,2,sortedRowMeans[2]);
		mqnormalized.set(2,2,sortedRowMeans[4]);
		mqnormalized.set(3,2,sortedRowMeans[3]);
		mqnormalized.set(4,2,sortedRowMeans[1]);
		
		m.quantileNormalizeColumns();
		
		for(int i = 0; i < m.getRowDimension(); i++) {
			for (int j = 0; j < m.getColumnDimension(); j++) {
				assertEquals("Value for ("+i+","+j+")", mqnormalized.get(i,j), m.get(i,j));
			}
		}
		
	}

}
