package nextgen.core.scripture;

import java.util.Set;
import java.util.TreeSet;
import org.apache.log4j.Logger;


public class CombinationGenerator {

	static Logger logger = Logger.getLogger(CombinationGenerator.class.getName());
	int totalPossible;
	int[][] permutations;
	
	/**
	 * Choose a value from each of the possible values in the array until all permutations are exhausted
	 * @param values an array containing the number of allowable elements for each position
	 */
	public CombinationGenerator(int[] values){
		this.totalPossible=getTotalPossible(values);
		int[][] matrix=generatePermutations(values);
		Set<String> matrixNum=getStringRep(matrix);
		if(matrixNum.size()!=this.totalPossible){throw new IllegalStateException();}
		this.permutations=matrix;
	}
	
	public int[][] getPermutations(){return this.permutations;}
	
	private Set<String> getStringRep(int[][] matrix) {
		Set<String> set=new TreeSet<String>();
		
		
		for(int row=0; row<matrix.length; row++){
			String rtrn="";
			for(int column=0; column<matrix[row].length; column++){
				rtrn+=matrix[row][column]+",";
			}
			set.add(rtrn);
		}
		return set;
		
	}

	private String print(int[][] matrix) {
		String rtrn="\n";
		for(int row=0; row<matrix.length; row++){
			for(int column=0; column<matrix[row].length; column++){
				rtrn+=matrix[row][column]+",";
			}
			rtrn+="\n";
		}
		return rtrn;
	}

	private int[][] generatePermutations(int[] values) {
		//create matrix with values.length columns, and numPerm rows
		int[][] matrix=new int[this.totalPossible][values.length];
		
		//iterate through position in values and in order replace the values of that position in matrix
		int totalSoFar=1;
		for(int i=0; i<values.length; i++){
			matrix=modify(matrix, i, values[i], totalSoFar);
			totalSoFar=totalSoFar*values[i];
		}
		
		//return matrix
		return matrix;
		
		/*int[] seed=new int[values.length];
		Collection<>
		//shuffle one positions at a time
		for(int i=0; i<values.length; i++){
			//enumerate all possible at i
			Collection<int[]> extraSeeds=enumerate(seed, i);
			//take these and augment by enumerating all possible at i+1
			
			
			for(int j=0; j<values.length; j++){
			//at position i, j, enumerate all possible values
				 rtrn.addAll(enumerate(values, seed, i, j));
			}
		}
		logger.info(rtrn.size());
		for(String l: rtrn){
			logger.info(l);
		}*/
	}

	/**
	 * Augment the ith column of the matrix with possible values
	 * @param matrix to augment
	 * @param totalSoFar 
	 * @param i column to modify
	 * @param possible values to use
	 */
	private int[][] modify(int[][] matrix, int column, int possibleValues, int totalSoFar) {
		
		
		
		int row=0;	
		//TODO Need to iterate over all rows
		
		while(row<matrix.length){
		for(int k=0; k<possibleValues; k++){
			for(int i=0; i<totalSoFar; i++){
				//iterate through each number by the total number so far
				//logger.info((row)+" "+column+" "+k);
				matrix[row][column]=k;
				row++;
			}
		}
		}
		
		return matrix;
	}

	

	private int getTotalPossible(int[] values) {
		int rtrn=1;
		for(int i=0; i<values.length; i++){
			rtrn=rtrn*values[i];
		}
		return rtrn;
	}

	public static void main(String[] args){
		int[] array={2,3,2,2};
		new CombinationGenerator(array);
		/*MultiSetPermutations perm=new MultiSetPermutations(array);
		while(perm.hasNext()){
			logger.info(perm.nextPermutation());
		}*/
	}
	
}
