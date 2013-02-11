package broad.core.datastructures;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import junit.framework.TestCase;


public class MatrixWithHeadersTest  extends TestCase{
	public void testGetRowName() {
		int numRows = 200000;
		char first = 65;
		char last = 90;
		int wordSize =(int) (Math.log(numRows)/Math.log(last - first )) + 1;
		char [] counter = new char[wordSize];
		List<String> rowNames = new ArrayList<String>(numRows);
		for(int i = 0; i < numRows; i++) {
			StringBuilder word = new StringBuilder(wordSize);
			for(int j = 0; j < wordSize; j++) {
				word.append((char)(counter[j] +  first));
			}
			rowNames.add(word.toString());
			//System.out.println(word);
			increment(counter, last - first + 1);
		}
		
		ArrayList<String> columns = new ArrayList<String>(1);
		columns.add("Col1");
		MatrixWithHeaders mwh = new MatrixWithHeaders(rowNames, columns);
		
		long start = System.nanoTime();
		String row0 = mwh.getRowName(0);
		long end = System.nanoTime();
		System.out.println("getRowName took " + (end - start) + " nanoseconds");
		assertEquals(rowNames.get(0), row0);
		assertEquals(rowNames.get(wordSize - 1), mwh.getRowName(wordSize - 1));
		
		Random r = new Random();
		for(int i = 0; i < 50; i++) {
			int idx = r.nextInt(wordSize - 1 );
			start = System.nanoTime();
			assertEquals(rowNames.get(idx), mwh.getRowName(idx));
			end = System.nanoTime();
			System.out.println("getRowName took " + (end - start) + " nanoseconds");
		}
		
	}
	
	
	private void increment(char [] counter,int max) {
		int idx = counter.length - 1;
		while(idx >=0) {
			counter[idx] =(char)( (counter[idx] + 1) % max);
			if(counter[idx] > 0) {
				break;
			}
			idx--;
		}
	}
}
