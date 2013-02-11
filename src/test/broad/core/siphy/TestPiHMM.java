package broad.core.siphy;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;


import Jama.Matrix;

import junit.framework.TestCase;
import broad.core.error.ParseException;
import broad.core.hmm.MarkovModel.BackwardResult;
import broad.core.hmm.MarkovModel.ForwardResult;
import broad.core.multiplealignment.MAFAlignment;
import broad.core.multiplealignment.MAFIO;
import broad.core.multiplealignment.MultipleAlignment;
import broad.core.siphy.EvolutionaryModel;
import broad.core.siphy.EvolutionaryModelParameters;
import broad.core.siphy.PiHMM;

public class TestPiHMM extends TestCase {

	public void testSegmentAlignment() throws Exception {
		EvolutionaryModel model = loadModel();
		MultipleAlignment aln = loadAlignment(model);
		double coverage = 0.08;
		PiHMM hmm = new PiHMM(model, 10, coverage);
		hmm.useDefaultConservedStates();
		
		short[] path = hmm.viterbiMostLikelyEstimation(aln);
		System.out.println("track name=\"Test Viterbi cov " + (100 * coverage) + "%\" color=200,0,100 visibility=2 type=wiggle_0");
		for(int i = 0; i < path.length; i++) {
			System.out.println("chr7\t"+(aln.getReferenceStart() + i) + "\t" + (aln.getReferenceStart() + i + 1) + "\t" + (path[i] > 0 ? 1 : 0));
		}
		
		ForwardResult<Map<String, Matrix>> forward = hmm.runForwardAlgorithm(aln);
		BackwardResult<Map<String,Matrix>> back    = forward.runBackwardAlgorithm();
		
		System.out.println("track name=\"Test Posterior cov " + (100 * coverage) + "%\" color=200,0,100 visibility=2 type=wiggle_0");
		for(int i = 0; i < aln.length(); i++) {
			System.out.println("chr7\t"+(aln.getReferenceStart() + i) + "\t" + (aln.getReferenceStart() + i + 1) + "\t"+ (1 - back.getPosteriorProbability(0, i)));
		}
		
	}
	private MultipleAlignment loadAlignment(EvolutionaryModel model) throws FileNotFoundException, IOException, ParseException {
		URL alnURL = getClass().getResource("ENm_26910761-26915760.maf");
		MAFIO mafio = new MAFIO();
		String [] seqs = model.getTree().getAllExternalSeqNames();
		List<String> seqsToLoad = new ArrayList<String>(seqs.length);
		for(int i =0; i < seqs.length; i++) {
			seqsToLoad.add(seqs[i]);
		}
		MAFAlignment rawAlignment =  mafio.load(alnURL.getFile(), seqsToLoad);
		rawAlignment.compress();
		MultipleAlignment alignment = rawAlignment.toMultipleAlignment();
		alignment.encodeAsMatrix();
		return alignment;
	}
	
	private EvolutionaryModel loadModel() throws IOException, ParseException {
		URL modelFileURL = getClass().getResource("combined.AR.new.out.mod");
		EvolutionaryModelParameters modelParams = new EvolutionaryModelParameters(new File(modelFileURL.getFile()));
		EvolutionaryModel model = new EvolutionaryModel(modelParams);
		return model;
	}
}
