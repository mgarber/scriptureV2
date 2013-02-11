/**
 * 
 */
package broad.pda.methylation;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.File;
import java.io.FileWriter;
import java.io.InputStream;
import java.util.*;

import broad.core.error.ParseException;
import broad.core.math.EmpiricalDistribution;
import broad.core.sequence.Sequence;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;
import broad.pda.methylation.PadlockProbe;
import broad.pda.methylation.OutputPadlock;
import broad.core.sequence.Sequence;

import net.sf.picard.util.TabbedInputParser;
import net.sf.samtools.util.StringUtil;
import net.sf.picard.fastq.*;

import org.apache.commons.lang3.StringUtils;

/**
 * @author engreitz
 *
 */
public class PadlockAnalysis {
	private Collection<PadlockProbe> probes;
	private FastqReader read1, read2;
	
	Map<String,PadlockProbe> probeMap, captureMap;
	Map<String,Integer> orphanTags, wrongCapture, tags;
	
	int circles, total;
	EmpiricalDistribution cPct, gPct;
	Map<PadlockProbe,Integer> probeCount, correctCaptureCount;
	
	String out;
	

	PadlockAnalysis(String probeFile, String read1, String read2, String out) throws FileNotFoundException {
		this.probes = loadPadlockProbes(probeFile);
		System.out.println("Loaded " + probes.size() + " unique padlock probes.");
		this.read1 = new FastqReader(new File(read1));
		this.read2 = new FastqReader(new File(read2));
		
		probeMap = initProbeMap(probes, 15);
		captureMap = initProbeMap(probes, 50);
		orphanTags = new HashMap<String,Integer>();
		wrongCapture = new HashMap<String,Integer>();
		
		cPct = new EmpiricalDistribution(100,0,1);
		gPct = new EmpiricalDistribution(100,0,1);
		
		circles = 0;
		total = 0;
		probeCount = initProbeCount(probes);
		correctCaptureCount = initProbeCount(probes);
		initTags();
		this.out = out;
		System.out.println(formatTabbedMap(tags));
	}
	
	
	private Map<PadlockProbe,Integer> initProbeCount(Collection<PadlockProbe> probes) {
		Map<PadlockProbe,Integer> map = new HashMap<PadlockProbe,Integer>();
		for (PadlockProbe p : probes) {
			map.put(p,0);
		}
		return map;
	}
	
	private void initTags() {
		tags = new HashMap<String,Integer>();
		tags.put(OutputPadlock.LEFT_AMP, 0);
		tags.put(OutputPadlock.RIGHT_AMP, 0);
		tags.put(Sequence.reverseSequence(OutputPadlock.LEFT_AMP), 0);
		tags.put(Sequence.reverseSequence(OutputPadlock.RIGHT_AMP), 0);
	}

	private Collection<PadlockProbe> loadPadlockProbes(String file) throws FileNotFoundException {
		InputStream is = new FileInputStream(file);
		TabbedInputParser ip = new TabbedInputParser(false, is);
		Iterator<String []> itr = ip.iterator();
		
		Collection<PadlockProbe> probes = new HashSet<PadlockProbe>();
		while (itr.hasNext()) {
			String[] fields = itr.next();
			PadlockProbe newProbe = new PadlockProbe(fields);
			boolean found = false;
			for (PadlockProbe p : probes) { 
				if (p.equals(newProbe)) {
					found = true;
					break;
				}
			}
			if (!found) probes.add(newProbe);
		}
		
		return probes;
	}

	
	/**
	 *  1. Probe representation
	 *  2. Whether probe captured the right sequence
	 *  3. Whether probe circularized and/or captured the tails
	 *  4. C:G ratio of captured sequence
	 */
	private void processRead1() {
		
		System.out.println("Processing Read 1");
		while (read1.hasNext()) {
			final FastqRecord seq = read1.next();
			total++;
			
			String head15 = seq.getReadString().substring(0,15);
			String head50 = seq.getReadString().substring(0,50);
			
			System.out.println(head15);
			System.out.println(head50);
			
			if (probeMap.containsKey(head15)) {
				PadlockProbe currProbe = probeMap.get(head15);
				probeCount.put(currProbe, probeCount.get(currProbe) + 1);
				//System.out.println("debug1");
				
				if (captureMap.containsKey(head50)) {
					correctCaptureCount.put(currProbe, correctCaptureCount.get(currProbe) + 1);
					//System.out.println("debug2");
				} else {
					String rest = head50.substring(currProbe.getRightCapture().length(), currProbe.getRightCapture().length() + 15);
					if (tags.containsKey(rest.substring(0,15))) {
						incrementKey(tags, rest.substring(0,15));
						//System.out.println("debug3");
					} else {
						incrementKey(wrongCapture, rest);
						//System.out.println("debug4");
					}
				}
				
			} else {
				incrementKey(orphanTags, head15);
				//System.out.println("debug5");
			}
			
			gPct.add(StringUtils.countMatches(seq.getReadString(), "G"));
			cPct.add(StringUtils.countMatches(seq.getReadString(), "C"));
		}
	}
	
	
	private <T> String formatTabbedMap(Map<T,Integer> map) {
		
		LinkedList list = new LinkedList(map.entrySet());
	    Collections.sort(list, new Comparator() {
	          public int compare(Object o1, Object o2) {
	               return ((Comparable) ((Map.Entry) (o1)).getValue())
	              .compareTo(((Map.Entry) (o2)).getValue());
	          }
	     });

	    StringBuilder sb = new StringBuilder("");
	    for (Iterator it = list.descendingIterator(); it.hasNext();) {
	        Map.Entry entry = (Map.Entry)it.next();
	        sb.append(entry.getKey() + "\t" + entry.getValue() + "\n");
	    }
		return sb.toString();
	}
	
	
	private <T> int sumCountMap(Map<T, Integer> map) {
		int count = 0;
		for (T key : map.keySet()) {
			count = count + map.get(key);
		}
		return count;
	}
	

	private void writeRead1(String out) throws IOException {
		System.out.println("Total reads examined: " + total);
		System.out.println("Correct captures: " + sumCountMap(correctCaptureCount));
		System.out.println("Tags: " + sumCountMap(tags));
		new FileWriter(new File(out + ".tags")).write(formatTabbedMap(tags));
		new FileWriter(new File(out + ".correctCaptureCount")).write(formatTabbedMap(correctCaptureCount));
		new FileWriter(new File(out + ".probeCount")).write(formatTabbedMap(probeCount));
		new FileWriter(new File(out + ".orphanTags")).write(formatTabbedMap(orphanTags));
		new FileWriter(new File(out + ".wrongCapture")).write(formatTabbedMap(wrongCapture));
	}
	
	
	private void incrementKey(Map<String,Integer> map, String key) {
		if (!map.containsKey(key))
			map.put(key, 1);
		else
			map.put(key, map.get(key) + 1);
	}
	
	
	private Map<String,PadlockProbe> initProbeMap(Collection<PadlockProbe> probes, int k) {
		Map<String,PadlockProbe> map = new HashMap<String,PadlockProbe>();
		for (PadlockProbe p : probes) {
			map.put(Sequence.reverseSequence(p.getCapturedSequence()).substring(0,k), p);
		}
		return map;
	}
	
	
	private static String USAGE = "PadlockAnalysis -in <probelist> -1 <read 1> -2 <read 2>\n";
	
	/**
	 * @param args
	 */
	public static void main(String[] args) throws FileNotFoundException, IOException, ParseException {
		ArgumentMap argmap = CLUtil.getParameters(args, USAGE, "write");
		
		final String in = argmap.getInput();
		final String read1 = argmap.getMandatory("1");
		final String read2 = argmap.getMandatory("2");
		final String out = argmap.getOutput();
		PadlockAnalysis p = new PadlockAnalysis(in, read1, read2, out);
		p.processRead1();
		p.writeRead1(out);
	}

}
