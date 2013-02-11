package broad.core.siphy;

import java.io.BufferedWriter;
import java.io.File;
import java.io.PrintWriter;
import java.net.URL;
import java.util.List;

import org.forester.phylogeny.Phylogeny;

import Jama.Matrix;
import broad.core.multiplealignment.MultipleAlignment;
import broad.core.siphy.TreeScaler;

import junit.framework.TestCase;

public class TestSimulations extends TestCase {
	public void testSimulateEvolution() throws Exception {		
		URL modelURL = getClass().getResource("autosomal.mod");
		EvolutionaryModelParameters modelParams = new EvolutionaryModelParameters(new File (modelURL.getPath()));	
		EvolutionaryModel model = new EvolutionaryModel(modelParams);

		Matrix piMat = new Matrix(model.alphabetSize, 1);
		piMat.set(0, 0, 0.003);
		piMat.set(1, 0, 0.003);
		piMat.set(2, 0, 0.004);
		piMat.set(3, 0, 0.99);

		PiHMM hmm = PiHMM.createTwoStateChain(100, 0.9, model, piMat);
		List<Integer> path = hmm.emitPath(100);
		System.out.println(path);
		
		MultipleAlignment ma = hmm.generateAlignment(path, "PHYLIP");
		BufferedWriter outBW = new BufferedWriter(new PrintWriter(System.out));
		ma.write(outBW);
		outBW.flush();
	}
}
