package broad.core.hmm;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import broad.core.hmm.MarkovModel;
import broad.core.hmm.MarkovState;
import broad.core.hmm.MarkovModel.BackwardResult;
import broad.core.hmm.MarkovModel.ForwardResult;


public class BlackBookDiceExampleTest extends junit.framework.TestCase{
	public static class RegularDice implements MarkovState<Integer> {
		static final double prob = 1/(double)6;
		static final double logProb = Math.log(prob);
		public double getEmissionLogProbability(Integer observation) {
			return logProb;
		}

		public double getEmissionProbability(Integer observation) {
			return 1/(double)6;
		}

		public String getName() {
			return "Regular Dice";
		}
		
		public Integer emitObservation() {
			double draw = new Random().nextDouble();			
			double cummProb =  prob;
			int ret = 1;
			while(cummProb < draw) {
				cummProb += prob;
				ret++;
			}
			return ret;
		}
		
	}
	
	public static class LoadedDice implements MarkovState<Integer> {
		static final double loaded6Prob = 1/(double)2; // Notice that in the black book this probability is 1/2.
		static final double loadedNon6Prob = 1/(double)10; //Notice that in the black book these probabilities are 1/10.
		static final double logLoaded6Prob = Math.log(loaded6Prob);
		static final double logLoadedNon6Prob = Math.log(loadedNon6Prob);
		public double getEmissionLogProbability(Integer observation) {
			return observation == 6 ? logLoaded6Prob : logLoadedNon6Prob;
		}

		public double getEmissionProbability(Integer observation) {
			return observation == 6 ? loaded6Prob: loadedNon6Prob;
		}

		public String getName() {
			return "Loaded Dice";
		}

		public Integer emitObservation() {
			double draw = new Random().nextDouble();
			int ret = 0;
			if(draw < 0.5) { ret = 6;}
			else {				
				double cummProb = 0.5 + loadedNon6Prob;
				ret = 1;
				while(cummProb < draw) {
					cummProb += loadedNon6Prob;
					ret++;
				}
			}
			return ret;
		}
		
	}
	
	public void testSmall() throws Exception {
		MarkovModel<Integer> mm = new MarkovModel<Integer>(2);  
		
		mm.addState(new RegularDice());
		mm.addState(new LoadedDice());
		
		mm.setStateTransitionProbability(0, 0, 0.95);
		mm.setStateTransitionProbability(0, 1, 0.05);
		mm.setStateTransitionProbability(1, 0, 0.1);
		mm.setStateTransitionProbability(1, 1, 0.9);
		
		//Make the starting and end state Loaded 
		mm.setInitialStateTransitionProbability(0, 0.1);
		mm.setInitialStateTransitionProbability(1, 0.9);
		
		mm.setEndStateTransitionProbability(0, 0.05);
		mm.setEndStateTransitionProbability(1, 0.9);
		
		List<Integer> rolls = new ArrayList<Integer>(12);
		rolls.add(4); rolls.add(5); rolls.add(3); rolls.add(1); rolls.add(3); rolls.add(2);
		rolls.add(6); rolls.add(5); rolls.add(1); rolls.add(2); rolls.add(4); rolls.add(5);
		rolls.add(6); //rolls.add(3); rolls.add(6); rolls.add(6); rolls.add(6); rolls.add(6); 
		
		short [] path = mm.viterbiMostLikelyEstimation(rolls);
		
		for(int i = 0; i < rolls.size(); i++) {
			System.out.print(rolls.get(i));
		}
		
		System.out.println("");
		for(int i = 0; i < rolls.size(); i++) {
			System.out.print(path[i] );
		}
		System.out.println("");
		
		short [] fairPath = new short[13];
		for(int i = 0; i < fairPath.length; i++) {fairPath[i] = 0;};
		fairPath[12] = 1;
		
		assertEquals("Likelihood of this fair path should be the same as Viterbi's",mm.computePathLogLikelihood(rolls, fairPath, 1, 1), mm.computePathLogLikelihood(rolls, path, 1, 1));
		for(int i =0; i < fairPath.length; i++) {
			assertEquals("fairPath step " + i + " should be the same as Viterbi's ", fairPath[i],path[i]);
		}


	}
	
	
	public void testViterbi() throws Exception {
		List<Integer> rolls = loadSampleRolls();
		
		MarkovModel<Integer> mm = new MarkovModel<Integer>(2);  
		
		mm.addState(new RegularDice());
		mm.addState(new LoadedDice());
		
		mm.setStateTransitionProbability(0, 0, 0.95);
		mm.setStateTransitionProbability(0, 1, 0.05);
		mm.setStateTransitionProbability(1, 0, 0.1);
		mm.setStateTransitionProbability(1, 1, 0.9);
		
		short [] path = mm.viterbiMostLikelyEstimation(rolls);
		/*for(int i = 0; i < rolls.size(); i++) {
			System.out.print(i + " ");
			if(i<10) {
				System.out.print("  ");
			} else if (i < 100) {
				System.out.print(" ");
			}
			
		}*/
		System.out.println("");
		for(int i = 0; i < rolls.size(); i++) {
			System.out.print(rolls.get(i) /*+ "   "*/);
		}
		
		System.out.println("");
		for(int i = 0; i < rolls.size(); i++) {
			if(i < 48) {
				assertEquals("Expected Fair but got loaded at roll " + i, 0, path[i] );
			} else if( i < 66) {
				assertEquals("Expected loaded but got fair at roll " + i, 1, path[i] );
			}  else if (i < 78) {
				assertEquals("Expected Fair but got loaded at roll " + i, 0, path[i] );
			} else if( i < 112) {
				assertEquals("Expected loaded but got fair at roll " + i, 1, path[i] );
			}  else if (i < 179) {
				assertEquals("Expected Fair but got loaded at roll " + i, 0, path[i] );
			} else if( i < 192) {
				assertEquals("Expected loaded but got fair at roll " + i, 1, path[i] );
			}  else if (i < 270) {
				assertEquals("Expected Fair but got loaded at roll " + i, 0, path[i] );
			} else if( i < 289) {
				assertEquals("Expected loaded but got fair at roll " + i, 1, path[i] );
			} else if (i < 300) {
				assertEquals("Expected Fair but got loaded at roll " + i, 0, path[i] );
			}   
			System.out.print(path[i] /*+ "   "*/);
		}
	}
	
	public void testBackwardAlgorithm() throws Exception {
		List<Integer> rolls = loadSampleRolls();
		
		MarkovModel<Integer> mm = new MarkovModel<Integer>(2);  
		
		mm.addState(new RegularDice());
		mm.addState(new LoadedDice());
		
		mm.setStateTransitionProbability(0, 0, 0.95);
		mm.setStateTransitionProbability(0, 1, 0.05);
		mm.setStateTransitionProbability(1, 0, 0.1);
		mm.setStateTransitionProbability(1, 1, 0.9);
		
		ForwardResult<Integer> forwardComputation = mm.runForwardAlgorithm(rolls);
		BackwardResult<Integer> backwardComputation = forwardComputation.runBackwardAlgorithm();
		System.out.println("\nRoll probability " + forwardComputation.getProbability());
		for(int i = 0; i < rolls.size(); i++) {

			System.out.println(i +"\t" + rolls.get(i) + "\t" + backwardComputation.getPosteriorProbability(0, i));
		}
		
	}
	
	

	private List<Integer> loadSampleRolls() throws FileNotFoundException, IOException {
		List<Integer> rolls = new ArrayList<Integer>(); 
		URL rollData = getClass().getResource("roll_sequence.txt");
		BufferedReader br = new BufferedReader(new FileReader(rollData.getFile()));
		String line = null;
		while((line = br.readLine()) != null) {
			String [] data = line.split("\\s");
			for(int i = 0; i < data.length; i++) {
				rolls.add(Integer.valueOf(data[i]));
			}
		}
		br.close();
		
		return rolls;
	}
}
