package nextgen.editing.crispr;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.concurrent.ConcurrentLinkedQueue;

import org.apache.log4j.Level;
import org.apache.log4j.Logger;

import broad.core.parser.CommandLineParser;
import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;
import broad.pda.annotation.BEDFileParser;

import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.Annotation.Strand;
import nextgen.core.annotation.Gene;
import nextgen.editing.RestrictionEnzyme;
import nextgen.editing.RestrictionEnzymeCutSite;
import nextgen.editing.RestrictionEnzymeCutSitePair;
import nextgen.editing.RestrictionEnzymeFactory;
import nextgen.editing.RestrictionEnzymePair;
import nextgen.editing.SingleCleavageTypeIIRestrictionEnzyme;
import nextgen.editing.TypeIISRestrictionEnzyme;
import nextgen.editing.crispr.predicate.GuideLacksEnzymeCutSite;
import nextgen.editing.crispr.predicate.GuidePairDoubleNickConfiguration;
import nextgen.editing.crispr.predicate.GuideProximityToNearestRegion;
import nextgen.editing.crispr.predicate.GuideSufficientEfficacy;
import nextgen.editing.crispr.predicate.GuideSufficientIsolation;
import nextgen.editing.crispr.score.GuideEfficacyScore;
import nextgen.editing.crispr.score.GuideOffTargetScore;
import nextgen.editing.crispr.score.GuidePairCombinedEfficacyDistanceScore;
import nextgen.editing.crispr.score.GuidePairInnerDistanceScore;

/**
 * Design guide RNA pairs for a double nick strategy of CRISPR editing
 * @author prussell
 */
public class DoubleNickCRISPRDesigner {
	
	private static Logger logger = Logger.getLogger(DoubleNickCRISPRDesigner.class.getName());
	private Map<String, Sequence> chrsByName;
	
	
	private static int MIN_UPSTREAM_DISTANCE = 1;
	private static int MAX_UPSTREAM_DISTANCE = 500;
	private static int MIN_DOWNSTREAM_DISTANCE = 500;
	private static int MAX_DOWNSTREAM_DISTANCE = 5000;
	private static int MIN_DIST_TO_NEAREST_GENE = 3000;
	private static int MAX_INNER_DIST_CUT_SITE_PAIRS = 150;
	private static int MAX_INNER_DIST_GUIDE_RNA_PAIRS = 40;
	private static int MAX_DIST_TO_RESTRICTION_SITE = 100;
	private static boolean ENFORCE_DOUBLE_NICK_CONFIGURATION = true;
	private static boolean ENFORCE_MIN_DIST_TO_NEAREST_GENE = true;
	private static boolean ENFORCE_DOWNSTREAM_PROXIMITY_TO_RESTRICTION_ENZYME = true;
	private static boolean ENFORCE_MIN_OFF_TARGET_SCORE = true;
	private static boolean ENFORCE_MAX_GUIDE_EFFICACY_SCORE = true;
	private static double MAX_GUIDE_EFFICACY_SCORE = 0.6;
	private static boolean WRITE_FAILED_PAIRS_FOR_MISSING_REGIONS = false;
	private static int MIN_OFF_TARGET_SCORE = 30;
	private static File OFF_TARGET_BITS = null;
	private static int NUM_BEST_PAIRS_TO_GET_PER_CUT = 20;
	private static boolean POOL_NON_OVERLAPPING = false;
	private static String LEFT_OLIGO_FLANKING_SEQUENCE;
	private static String RIGHT_OLIGO_FLANKING_SEQUENCE;
	private static Collection<RestrictionEnzyme> FORBIDDEN_ENZYMES = new ArrayList<RestrictionEnzyme>();
	private static int MAX_NUM_POOLS = 12;
	private static String PRIMER3_CORE;
	private FileWriter failedPairBedWriter;
	private FileWriter failedPairTableWriter;
	private Collection<NickingGuideRNAPair> allValidPairs;
	
	/**
	 * Full gene annotation
	 */
	private Map<String, Collection<Gene>> annotation;
	
	/**
	 * The genes to target with CRISPR system
	 */
	private Map<String, Collection<Gene>> targetGenes;

	/**
	 * Collection of restriction enzyme cut site pairs in window of interest downstream of each gene to edit
	 */
	private Map<Gene, Collection<RestrictionEnzymeCutSitePair>> downstreamRestrictionEnzymeCutSitePairsByGene;
	
	/**
	 * Collection of restriction enzyme cut sites in window of interest downstream of each gene to edit
	 */
	private Map<Gene, Collection<RestrictionEnzymeCutSite>> downstreamRestrictionEnzymeCutSitesByGene;

	/**
	 * @param genomeFasta Genome fasta file
	 * @param bedFileTargetGenes Bed file of target genes
	 * @param bedFileFullAnnotation Bed file of genome annotation
	 * @throws IOException
	 */
	private DoubleNickCRISPRDesigner(String genomeFasta, String bedFileTargetGenes, String bedFileFullAnnotation) throws IOException {
		chrsByName = FastaSequenceIO.getChrSequencesFromFasta(genomeFasta);
		targetGenes = BEDFileParser.loadDataByChr(new File(bedFileTargetGenes));
		annotation = BEDFileParser.loadDataByChr(new File(bedFileFullAnnotation));
	}
	
	/**
	 * Whether the guide RNA pair passes a set of filters common to upstream and downstream pairs
	 * @param target Target gene for the guide pair
	 * @param pair The guide RNA pair
	 * @param minDistToNearestGene Minimum required distance to nearest gene other than target gene
	 * @param guideEfficacy Guide efficacy object, preferably with score instantiated with all guide RNAs of interest (faster), or null if not using
	 * @return Whether the pair passes the filters
	 * @throws InterruptedException 
	 * @throws IOException 
	 */
	private boolean guidePairPassesAllBasicFilters(Gene target, NickingGuideRNAPair pair, GuideSufficientEfficacy guideEfficacy) throws IOException, InterruptedException {
		
		if(ENFORCE_DOUBLE_NICK_CONFIGURATION) {
			// Check that the pair are arranged correctly
			GuidePairDoubleNickConfiguration dnc = new GuidePairDoubleNickConfiguration(MAX_INNER_DIST_GUIDE_RNA_PAIRS);
			if(!dnc.evaluate(pair)) {
				//logger.debug("PAIR_FAILS_DOUBLE_NICK_CONFIGURATION\t" + pair.toString());
				return false;
			}
		}
				
		if(ENFORCE_MAX_GUIDE_EFFICACY_SCORE) {
			if(guideEfficacy == null) {
				throw new IllegalArgumentException("Must pass valid guide efficacy score object");
			}
			if(!guideEfficacy.evaluate(pair.getLeftGuideRNA())) {
				logger.debug("LEFT_GUIDE_FAILS_GUIDE_EFFICACY\t" + pair.toString());
				return false;
			}
			if(!guideEfficacy.evaluate(pair.getRightGuideRNA())) {
				logger.debug("RIGHT_GUIDE_FAILS_GUIDE_EFFICACY\t" + pair.toString());
				return false;
			}
		} else {
			if(guideEfficacy != null) {
				logger.warn("You provided a guide efficacy object but it is not being used because ENFORCE_MAX_GUIDE_EFFICACY_SCORE is false");
			}
		}
		
		if (ENFORCE_MIN_OFF_TARGET_SCORE) {
			if (pair.getLeftGuideRNA().getScore() < MIN_OFF_TARGET_SCORE || 
				pair.getRightGuideRNA().getScore() < MIN_OFF_TARGET_SCORE) {
				return false;
			}
		}
		
		if(!FORBIDDEN_ENZYMES.isEmpty()) {
			GuideLacksEnzymeCutSite g = new GuideLacksEnzymeCutSite(FORBIDDEN_ENZYMES);
			if(!g.evaluate(pair.getLeftGuideRNA())) {
				logger.debug("LEFT_GUIDE_CONTAINS_CUT_SITE\t" + pair.toString()); 
				return false;
			}
			if(!g.evaluate(pair.getRightGuideRNA())) {
				logger.debug("RIGHT_GUIDE_CONTAINS_CUT_SITE\t" + pair.toString()); 
				return false;
			}
		}
		
		logger.debug("PAIR_PASSES_BASIC_FILTERS\t" + pair.toString());
		return true;
		
	}
	
	private void writeValidGuidePairsAllGenes(String outFilePrefix) throws IOException, InterruptedException {
		
		if(allValidPairs.isEmpty()) {
			logger.warn("No valid guide pairs have been found. Make sure findValidPairs() was called.");
		}
		
		if(WRITE_FAILED_PAIRS_FOR_MISSING_REGIONS) {
			failedPairBedWriter = new FileWriter(outFilePrefix + "_guide_pairs_failing_filters.bed");
			failedPairTableWriter = new FileWriter(outFilePrefix + "_guide_pairs_failing_filters_oligos.out");
			String header = "target_gene\toligo_ID\t" + NickingGuideRNAPair.getOligoFieldNames();
			failedPairTableWriter.write(header + "\n");
		}
				
		String oligoTable = outFilePrefix + "_oligos.out";
		logger.info("Writing oligos to " + oligoTable);
		Collection<Collection<NickingGuideRNAPair>> pools = new ArrayList<Collection<NickingGuideRNAPair>>();
		if(POOL_NON_OVERLAPPING) {
			pools.addAll(NickingGuideRNAPair.poolNonOverlapping(allValidPairs, MAX_NUM_POOLS));
			String bedPrefix = outFilePrefix + "_guide_pairs";
			NickingGuideRNAPair.writeBED(pools, bedPrefix);

		} else {
			pools.add(allValidPairs);
			String bedFile = outFilePrefix + "_guide_pairs.bed";
			NickingGuideRNAPair.writeBED(allValidPairs, bedFile, false);
		}
		NickingGuideRNAPair.writeOligoTable(pools, LEFT_OLIGO_FLANKING_SEQUENCE, RIGHT_OLIGO_FLANKING_SEQUENCE, oligoTable, PRIMER3_CORE);
		
		if(WRITE_FAILED_PAIRS_FOR_MISSING_REGIONS) {
			failedPairBedWriter.close();
			failedPairTableWriter.close();
		}
		
	}
	
	/**
	 * Establish restriction enzymes for downstream proximity filter
	 * @param listFileSingleEnzymes List of restriction enzymes for single cut site (one enzyme name per line)
	 * @param listFileLeftPairedEnzymes List of left enzymes for paired cut sites (one enzyme name per line)
	 * @param listFileRightPairedEnzymes List of right enzymes for paired cut sites (one enzyme name per line)
	 * @param outFilePrefix Output file prefix for bed file of cut site locations
	 * @throws IOException
	 */
	private void useRestrictionEnzymes(String listFileSingleEnzymes, String listFileEnzymePairs, String outFilePrefix) throws IOException {
		Collection<RestrictionEnzyme> singleEnzymes = new ArrayList<RestrictionEnzyme>();
		Collection<RestrictionEnzymePair> pairedEnzymes = new ArrayList<RestrictionEnzymePair>();
		for(SingleCleavageTypeIIRestrictionEnzyme enzyme : RestrictionEnzymeFactory.readFromFile(listFileSingleEnzymes)) {
			singleEnzymes.add((RestrictionEnzyme) enzyme);
		}
		for(RestrictionEnzymePair pair : RestrictionEnzymeFactory.readPairsFromFile(listFileEnzymePairs)) {
			pairedEnzymes.add(pair);
		}
		findAndWriteDownstreamRestrictionEnzymeCutSitesAllGenes(singleEnzymes, pairedEnzymes, outFilePrefix);
	}
	
	/**
	 * Find all single cut sites and paired cut sites downstream of all target genes
	 * Save to data structure
	 * Also write to bed file
	 * @param singleEnzymes Enzymes for single cut sites
	 * @param leftPairedEnzymes Left enzymes for paired cut sites
	 * @param rightPairedEnzymes Right enzymes for paired cut sites
	 * @param outFilePrefix Output file prefix
	 * @throws IOException
	 */
	private void findAndWriteDownstreamRestrictionEnzymeCutSitesAllGenes(Collection<RestrictionEnzyme> singleEnzymes, Collection<RestrictionEnzymePair> pairedEnzymes, String outFilePrefix) throws IOException {
		String file = outFilePrefix + "_restriction_enzyme_sites.bed";
		logger.info("Finding restriction enzyme sites downstream of all target genes and writing to file " + file);
		downstreamRestrictionEnzymeCutSitePairsByGene = new TreeMap<Gene, Collection<RestrictionEnzymeCutSitePair>>();
		downstreamRestrictionEnzymeCutSitesByGene = new TreeMap<Gene, Collection<RestrictionEnzymeCutSite>>();
		FileWriter writer = new FileWriter(file);
		for(String chr : targetGenes.keySet()) {
			for(Gene gene : targetGenes.get(chr)) {
				Strand strand = gene.getOrientation();
				if(strand.equals(Strand.UNKNOWN)) {
					writer.close();
					throw new IllegalArgumentException("Strand must be known");
				}
				int end1 = strand.equals(Strand.POSITIVE) ? gene.getEnd() + MIN_DOWNSTREAM_DISTANCE : gene.getStart() - MIN_DOWNSTREAM_DISTANCE;
				int end2 = strand.equals(Strand.POSITIVE) ? gene.getEnd() + MAX_DOWNSTREAM_DISTANCE : gene.getStart() - MAX_DOWNSTREAM_DISTANCE;
				int windowStart = Math.min(end1, end2);
				int windowEnd = Math.max(end1, end2);
				// Get single enzyme cut sites for this gene
				Collection<RestrictionEnzymeCutSite> singleSites = new ArrayList<RestrictionEnzymeCutSite>();
				for(RestrictionEnzyme enzyme : singleEnzymes) {
					Collection<RestrictionEnzymeCutSite> sites = RestrictionEnzymeCutSite.getAllSites(enzyme, chrsByName.get(chr), windowStart, windowEnd);
					singleSites.addAll(sites);
					for(RestrictionEnzymeCutSite site : sites) {
						writer.write(site.toBED() + "\n");
					}
				}
				downstreamRestrictionEnzymeCutSitesByGene.put(gene, singleSites);
				// Get paired cut sites for this gene
				Collection<RestrictionEnzymeCutSitePair> pairedSites = new ArrayList<RestrictionEnzymeCutSitePair>();
				for(RestrictionEnzymePair pair : pairedEnzymes) {
					RestrictionEnzyme left = pair.getLeft();
					RestrictionEnzyme right = pair.getRight();
					Collection<RestrictionEnzymeCutSitePair> sites = RestrictionEnzymeCutSitePair.getAllSitePairs(left, right, chrsByName.get(chr), windowStart, windowEnd, MAX_INNER_DIST_CUT_SITE_PAIRS);
					pairedSites.addAll(sites);
					for(RestrictionEnzymeCutSitePair cutPair : sites) {
						writer.write(cutPair.toBED() + "\n");
					}
				}
				downstreamRestrictionEnzymeCutSitePairsByGene.put(gene, pairedSites);
			}
		}
		writer.close();
	}
	
	/**
	 * Pick out the guide RNA pairs with best scores (score hardcoded for now)
	 * Make sure they can be combined into pools of pairwise non-overlapping guide RNAs
	 * Throws exception if can't be placed into required number of pools
	 * @param guidePairs Full set of guide RNA pairs
	 * @return The best pairs, or the full set if already smaller than number to get
	 * @throws IOException
	 * @throws InterruptedException
	 */
	private static Collection<NickingGuideRNAPair> getPoolableBestGuidePairs(Collection<NickingGuideRNAPair> guidePairs) throws IOException, InterruptedException {
		//long start = System.currentTimeMillis();
		if(guidePairs.size() == 0) {
			return guidePairs;
		}
		if(NUM_BEST_PAIRS_TO_GET_PER_CUT < 1) {
			throw new IllegalArgumentException("Must get at least 1 guide RNA pair.");
		}
		if(MAX_NUM_POOLS < 1) {
			throw new IllegalArgumentException("Max number of pools must be at least 1.");
		}
		Collection<GuideRNA> indGuides = new ArrayList<GuideRNA>();
		for(NickingGuideRNAPair pair : guidePairs) {
			indGuides.add(pair.getLeftGuideRNA());
			indGuides.add(pair.getRightGuideRNA());
		}
		GuideEfficacyScore score = new GuideEfficacyScore(indGuides);
		GuidePairCombinedEfficacyDistanceScore combinedScore = new GuidePairCombinedEfficacyDistanceScore(score);
		//GuidePairInnerDistanceScore innerDistScore = new GuidePairInnerDistanceScore();
		for(NickingGuideRNAPair pair : guidePairs) {
			pair.setScore(combinedScore.getScore(pair));
		}
		ArrayList<NickingGuideRNAPair> sortedPairs = NickingGuideRNAPair.sortByScoreNondescending(guidePairs);
		// Lower scores are better
		Iterator<NickingGuideRNAPair> ascendingIter = sortedPairs.iterator();
		Collection<NickingGuideRNAPair> rtrn = new ArrayList<NickingGuideRNAPair>();
		// Try to pool
		Collection<Collection<NickingGuideRNAPair>> pools = new ArrayList<Collection<NickingGuideRNAPair>>();
		for(int i = 0; i < MAX_NUM_POOLS; i++) {
			pools.add(new TreeSet<NickingGuideRNAPair>());
		}
		int i = 0;
		while(i < NUM_BEST_PAIRS_TO_GET_PER_CUT) {
			try {
				NickingGuideRNAPair next = ascendingIter.next();
				boolean inPool = false;
				for(Collection<NickingGuideRNAPair> pool : pools) {
					boolean ok = true;
					for(NickingGuideRNAPair other : pool) {
						if(next.overlaps(other)) {
							ok = false;
							break;
						}
					} if(ok) {
						pool.add(next);
						inPool = true;
						break;
					}
				}
				if(!inPool) {
					// Can't be added to a pool; skip this one
					//logger.warn("Couldn't add " + next.toString() + " to a pool");
					continue;
				}
				//logger.info("Added pair. Guide efficacy scores " + score.getScore(next.getLeftGuideRNA()) + "," + score.getScore(next.getRightGuideRNA()) + "; inner distance " + next.getInnerDistance() + "; combined score " + combinedScore.getScore(next));
				i++;
				rtrn.add(next);
			} catch(NoSuchElementException e) {
				//long sec = (System.currentTimeMillis() - start) / 1000;
				//logger.info("Picked " + rtrn.size() + " poolable most effective and closest guide RNA pairs out of " + guidePairs.size() + " total pairs. Took " + sec + " seconds.");
				return rtrn;
			}
		}
		//long sec = (System.currentTimeMillis() - start) / 1000;
		//logger.info("Picked " + rtrn.size() + " poolable most effective and closest guide RNA pairs out of " + guidePairs.size() + " total pairs. Took " + sec + " seconds.");
		return rtrn;
	}
	
	
	private class ValidPairFinder implements Runnable {
		
		private ConcurrentLinkedQueue<Gene> geneQueue;
		
		public ValidPairFinder(ConcurrentLinkedQueue<Gene> genes) {
			geneQueue = genes;
		}
		
		@Override
		public void run() {
			while(!geneQueue.isEmpty()) {
				Gene gene = geneQueue.poll();
				try {
					synchronized(logger) {
						logger.info("Finding guide RNA pairs for gene " + gene.getName());
					}
					Collection<NickingGuideRNAPair> upstream = getPoolableBestGuidePairs(findAllValidPairsUpstreamOfTranscriptionStart(gene));
					if(upstream.size() != NUM_BEST_PAIRS_TO_GET_PER_CUT) {
						synchronized(logger) {
							logger.warn("Could only get " + upstream.size() + " good poolable pairs upstream of " + gene.getName());
						}
					}
					Collection<NickingGuideRNAPair> downstream = getPoolableBestGuidePairs(findAllValidPairsDownstreamOfTranscriptionStop(gene));
					if(downstream.size() != NUM_BEST_PAIRS_TO_GET_PER_CUT) {
						synchronized(logger) {
							logger.warn("Could only get " + downstream.size() + " good poolable pairs downstream of " + gene.getName());
						}
					}
					synchronized(allValidPairs) {
						allValidPairs.addAll(upstream);
						allValidPairs.addAll(downstream);
					}
				} catch (IOException e) {
					synchronized(logger) {
						logger.warn("Caught exception, skipping gene " + gene.getName());
					}
					e.printStackTrace();
				} catch (InterruptedException e) {
					synchronized(logger) {
						logger.warn("Caught exception, skipping gene " + gene.getName());
					}
					e.printStackTrace();
				}
			}
		}
		
	}
	
	
	private void findValidPairsAllGenes(int numThreads) throws IOException, InterruptedException {
		
		logger.info("");
		logger.info("Finding valid guide RNA pairs for all genes...");
		
		allValidPairs = new ArrayList<NickingGuideRNAPair>();

		Collection<Gene> genes = new TreeSet<Gene>();
		for(String chr : targetGenes.keySet()) {
			genes.addAll(targetGenes.get(chr));
		}
		
		ConcurrentLinkedQueue<Gene> geneQueue = new ConcurrentLinkedQueue<Gene>(genes);
		
		Collection<Thread> threads = new ArrayList<Thread>();
		
		for(int i = 0; i < numThreads; i++) {
			ValidPairFinder pf = new ValidPairFinder(geneQueue);
			Thread t = new Thread(pf);
			threads.add(t);
			t.start();
		}
		
		Thread.sleep(1000);
		for(Thread t : threads) {
			t.join();
		}
		
		logger.info("Found " + allValidPairs.size() + " total valid guide pairs.");
	
	}
	
	/**
	 * Find all guide RNA pairs in a window downstream of transcription stop that pass all the downstream pair filters
	 * @param gene The gene
	 * @return Collection of valid guide RNA pairs
	 * @throws InterruptedException 
	 * @throws IOException 
	 */
	private Collection<NickingGuideRNAPair> findAllValidPairsDownstreamOfTranscriptionStop(Gene gene) throws IOException, InterruptedException {		
		
		Collection<NickingGuideRNAPair> allPairs = findAllPossibleGuideRNAsDownstreamOfTranscriptionStop(gene);
		if(allPairs.isEmpty()) {
			return new ArrayList<NickingGuideRNAPair>();
		}
		
		GuideSufficientEfficacy ge = null;
		if(ENFORCE_MAX_GUIDE_EFFICACY_SCORE) {
			ge = new GuideSufficientEfficacy(new GuideEfficacyScore(NickingGuideRNAPair.getIndividualGuideRNAs(allPairs)), MAX_GUIDE_EFFICACY_SCORE);
		}
		
		GuideOffTargetScore scorer = null;
		if (ENFORCE_MIN_OFF_TARGET_SCORE) {
			scorer = new GuideOffTargetScore(OFF_TARGET_BITS);
			for (GuideRNA guide : NickingGuideRNAPair.getIndividualGuideRNAs(allPairs)) {
				guide.setScore(scorer.getScore(guide));
			}
		}
		
		synchronized (logger) {
			logger.debug("Before filters there are " + allPairs.size()
					+ " pairs downstream of transcription stop.");
		}
		Collection<NickingGuideRNAPair> rtrn = new ArrayList<NickingGuideRNAPair>();
		
		// Collect restriction enzyme sites downstream of the gene
		Collection<Annotation> restrictionSites = new ArrayList<Annotation>();
		if(ENFORCE_DOWNSTREAM_PROXIMITY_TO_RESTRICTION_ENZYME) {
			synchronized (downstreamRestrictionEnzymeCutSitesByGene) {
				restrictionSites
						.addAll(RestrictionEnzymeCutSite
								.asAnnotations(downstreamRestrictionEnzymeCutSitesByGene
										.get(gene)));
				restrictionSites
						.addAll(RestrictionEnzymeCutSitePair
								.asAnnotations(downstreamRestrictionEnzymeCutSitePairsByGene
										.get(gene)));
			}
		}
		
		for(NickingGuideRNAPair pair : allPairs) {
			boolean passes = true;
			if(!guidePairPassesAllBasicFilters(gene, pair, ge)) {
				passes = false;
			}
			
			if(passes) {
				if(ENFORCE_DOWNSTREAM_PROXIMITY_TO_RESTRICTION_ENZYME) {
					GuideProximityToNearestRegion p = new GuideProximityToNearestRegion(restrictionSites, MAX_DIST_TO_RESTRICTION_SITE, "restriction_enzyme");
					boolean ok = p.evaluate(pair.getLeftGuideRNA()) || p.evaluate(pair.getRightGuideRNA());
					if(!ok) {
						passes = false;
					}
				}
			}
			
			if(passes) {
				if(ENFORCE_MIN_DIST_TO_NEAREST_GENE) {
					// Check that the pair is not too close to another gene
					// Ignore all overlappers of the target gene
					GuideSufficientIsolation si = new GuideSufficientIsolation(annotation, gene, MIN_DIST_TO_NEAREST_GENE, true, "distance_from_nearest_downstream_gene");
					if(!si.evaluate(pair.getLeftGuideRNA())) {
						//logger.debug("LEFT_GUIDE_FAILS_MIN_DIST_TO_NEAREST_GENE\t" + pair.toString());
						passes = false;
					}
					if(!si.evaluate(pair.getRightGuideRNA())) {
						//logger.debug("RIGHT_GUIDE_FAILS_MIN_DIST_TO_NEAREST_GENE\t" + pair.toString());
						passes = false;
					}
				}
			}
			
			if(passes) {
				synchronized (logger) {
					logger.debug("PAIR_PASSES_ALL_DOWNSTREAM_FILTERS\t"
							+ pair.toString());
				}
				rtrn.add(pair);
			}
		}
		if(rtrn.isEmpty()) {
			logger.warn("NO_VALID_DOWNSTREAM_GUIDE_PAIRS\t" + gene.getName());
			
			if(WRITE_FAILED_PAIRS_FOR_MISSING_REGIONS) {
				
				synchronized (logger) {
					logger.warn("WRITING_FAILED_PAIRS\t" + gene.getName());
				}
				GuidePairDoubleNickConfiguration dnc = new GuidePairDoubleNickConfiguration(MAX_INNER_DIST_GUIDE_RNA_PAIRS);
				
				for(NickingGuideRNAPair pair : allPairs) {
					
					if(!dnc.evaluate(pair)) {
						continue;
					}
					
					String bedName = pair.toString() + ":" + getFailureMessageBasicFilters(pair, ge);
					
					if(ENFORCE_DOWNSTREAM_PROXIMITY_TO_RESTRICTION_ENZYME) {
						GuideProximityToNearestRegion p = new GuideProximityToNearestRegion(restrictionSites, MAX_DIST_TO_RESTRICTION_SITE, "proximity_to_restriction_enzyme_site");
						if(!p.evaluate(pair.getLeftGuideRNA())) {
							bedName += p.getShortFailureMessage(pair.getLeftGuideRNA()) + ":";
						}
						if(!p.evaluate(pair.getRightGuideRNA())) {
							bedName += p.getShortFailureMessage(pair.getRightGuideRNA()) + ":";
						}
					}

					if(ENFORCE_MIN_DIST_TO_NEAREST_GENE) {
						// Check that the pair is not too close to another gene
						// Ignore all overlappers of the target gene
						GuideSufficientIsolation si = new GuideSufficientIsolation(annotation, gene, MIN_DIST_TO_NEAREST_GENE, true, "distance_from_nearest_downstream_gene");
						if(!si.evaluate(pair.getLeftGuideRNA())) {
							bedName += si.getShortFailureMessage(pair.getLeftGuideRNA()) + ":";
						}
						if(!si.evaluate(pair.getRightGuideRNA())) {
							bedName += si.getShortFailureMessage(pair.getRightGuideRNA()) + ":";
						}
					}
					
					pair.setName(bedName);
					synchronized (failedPairBedWriter) {
						failedPairBedWriter.write(pair.toBED() + "\n");
					}
					synchronized (failedPairTableWriter) {
						failedPairTableWriter.write(gene.getName() + "\t"
								+ pair.toString() + "\t" + pair.getOligos()
								+ "\n");
					}

				}
			}
			

		}
		/*synchronized (logger) {
			logger.info("Found " + rtrn.size() + " pairs downstream of "
					+ gene.getName());
		}*/
		return rtrn;
	}
	
	/**
	 * Find all guide RNA pairs in a window upstream of transcription start that pass all the upstream pair filters
	 * @param gene The gene
	 * @return Collection of valid guide RNA pairs
	 * @throws InterruptedException 
	 * @throws IOException 
	 */
	private  Collection<NickingGuideRNAPair> findAllValidPairsUpstreamOfTranscriptionStart(Gene gene) throws IOException, InterruptedException {

		Collection<NickingGuideRNAPair> allPairs = findAllPossibleGuideRNAsUpstreamOfTranscriptionStart(gene);
		if(allPairs.isEmpty()) {
			return new ArrayList<NickingGuideRNAPair>();
		}
		
		GuideSufficientEfficacy ge = null;
		if(ENFORCE_MAX_GUIDE_EFFICACY_SCORE) {
			ge = new GuideSufficientEfficacy(new GuideEfficacyScore(NickingGuideRNAPair.getIndividualGuideRNAs(allPairs)), MAX_GUIDE_EFFICACY_SCORE);
		}
		
		GuideOffTargetScore scorer = null;
		if (ENFORCE_MIN_OFF_TARGET_SCORE) {
			scorer = new GuideOffTargetScore(OFF_TARGET_BITS);
			for (GuideRNA guide : NickingGuideRNAPair.getIndividualGuideRNAs(allPairs)) {
				guide.setScore(scorer.getScore(guide));
			}
		}
				

		
		Collection<NickingGuideRNAPair> rtrn = new ArrayList<NickingGuideRNAPair>();
		for(NickingGuideRNAPair pair : allPairs) {
			boolean passes = true;
			if(!guidePairPassesAllBasicFilters(gene, pair, ge)) {
				passes = false;
			}
			
			if(passes) {
				synchronized (logger) {
					logger.debug("PAIR_PASSES_ALL_UPSTREAM_FILTERS\t" + pair.toString());
				}
				rtrn.add(pair);
			}
		}
		if(rtrn.isEmpty()) {
			synchronized (logger) {
				logger.warn("NO_VALID_UPSTREAM_GUIDE_PAIRS\t" + gene.getName());
			}
			if(WRITE_FAILED_PAIRS_FOR_MISSING_REGIONS) {
				synchronized (logger) {
					logger.warn("WRITING_FAILED_PAIRS\t" + gene.getName());
				}
				GuidePairDoubleNickConfiguration dnc = new GuidePairDoubleNickConfiguration(MAX_INNER_DIST_GUIDE_RNA_PAIRS);
				for(NickingGuideRNAPair pair : allPairs) {
					if(!dnc.evaluate(pair)) {
						continue;
					}
					String bedName = pair.toString() + ":" + getFailureMessageBasicFilters(pair, ge);
					pair.setName(bedName);
					synchronized (failedPairBedWriter) {
						failedPairBedWriter.write(pair.toBED() + "\n");
					}
					synchronized (failedPairTableWriter) {
						failedPairTableWriter.write(gene.getName() + "\t"
								+ pair.toString() + "\t" + pair.getOligos()
								+ "\n");
					}
				}
			}
			
		}
		/*synchronized (logger) {
			logger.info("Found " + rtrn.size() + " pairs upstream of "
					+ gene.getName());
		}*/
		return rtrn;
	}
	
	/**
	 * Get a string describing which basic filters the guide pair fails
	 * @param pair Guide pair
	 * @param ge Guide efficacy object or null if not using
	 * @return String with all failed filters or empty string if passes all basic filters
	 * @throws InterruptedException 
	 * @throws IOException 
	 */
	private static String getFailureMessageBasicFilters(NickingGuideRNAPair pair, GuideSufficientEfficacy ge) throws IOException, InterruptedException {
		
		String rtrn = "";
		
		if(ENFORCE_DOUBLE_NICK_CONFIGURATION) {
			// Check that the pair are arranged correctly
			GuidePairDoubleNickConfiguration dnc = new GuidePairDoubleNickConfiguration(MAX_INNER_DIST_GUIDE_RNA_PAIRS);
			if(!dnc.evaluate(pair)) {
				rtrn += dnc.getShortFailureMessage(pair) + ":";
			}
		}
				
		if(ENFORCE_MAX_GUIDE_EFFICACY_SCORE) {
			if(ge == null) {
				throw new IllegalArgumentException("Must pass valid guide efficacy score object");
			}
			if(!ge.evaluate(pair.getLeftGuideRNA())) {
				rtrn += "left_" + ge.getShortFailureMessage(pair.getLeftGuideRNA()) + ":";
			}
			if(!ge.evaluate(pair.getRightGuideRNA())) {
				rtrn += "right_" + ge.getShortFailureMessage(pair.getRightGuideRNA()) + ":";
			}
		} else {
			if(ge == null) {
				logger.warn("You provided a guide efficacy object but it is not being used because ENFORCE_MAX_GUIDE_EFFICACY_SCORE is false");
			}
		}
		
		if(ENFORCE_MIN_OFF_TARGET_SCORE) {
			// TODO off-target score filter
		}

		return rtrn;
		
	}

	
	
	/**
	 * Find all guide RNA pairs in a window downstream of transcription stop
	 * @param gene Gene
	 * @return All guide RNA pairs fully contained in the window
	 */
	private Collection<NickingGuideRNAPair> findAllPossibleGuideRNAsDownstreamOfTranscriptionStop(Gene gene) {
		Sequence chr = chrsByName.get(gene.getChr());
		Strand strand = gene.getOrientation();
		int geneStart = gene.getStart();
		int geneEnd = gene.getEnd();
		if(strand.equals(Strand.UNKNOWN)) {
			throw new IllegalArgumentException("Gene strand must be known");
		}
		int windowStrandedBegin = strand.equals(Strand.POSITIVE) ? geneEnd + MIN_DOWNSTREAM_DISTANCE : geneStart - MIN_DOWNSTREAM_DISTANCE;
		int windowStrandedEnd = strand.equals(Strand.POSITIVE) ? geneEnd + MAX_DOWNSTREAM_DISTANCE : geneStart - MAX_DOWNSTREAM_DISTANCE;
		int windowStart = Math.min(windowStrandedBegin, windowStrandedEnd);
		int windowEnd = Math.max(windowStrandedBegin, windowStrandedEnd);
		logger.debug("DOWNSTREAM_WINDOW\t" + gene.getName() + "\t" + gene.toUCSC() + ":" + gene.getOrientation().toString() + "\tmin_dist=" + MIN_DOWNSTREAM_DISTANCE + "\tmax_dist=" + MAX_DOWNSTREAM_DISTANCE + "\t" + chr.getId() + ":" + windowStart + "-" + windowEnd);
		return NickingGuideRNAPair.findAll(chr, windowStart, windowEnd, gene);
	}
	
	/**
	 * Find all guide RNA pairs in a window upstream of transcription start
	 * @param gene Gene
	 * @return All guide RNA pairs fully contained in the window
	 */
	private Collection<NickingGuideRNAPair> findAllPossibleGuideRNAsUpstreamOfTranscriptionStart(Gene gene) {
		Sequence chr = chrsByName.get(gene.getChr());
		Strand strand = gene.getOrientation();
		int geneStart = gene.getStart();
		int geneEnd = gene.getEnd();
		if(strand.equals(Strand.UNKNOWN)) {
			throw new IllegalArgumentException("Gene strand must be known");
		}
		int windowStrandedBegin = strand.equals(Strand.POSITIVE) ? geneStart - MAX_UPSTREAM_DISTANCE : geneEnd + MAX_UPSTREAM_DISTANCE;
		int windowStrandedEnd = strand.equals(Strand.POSITIVE) ? geneStart - MIN_UPSTREAM_DISTANCE : geneEnd + MIN_UPSTREAM_DISTANCE;
		int windowStart = Math.min(windowStrandedBegin, windowStrandedEnd);
		int windowEnd = Math.max(windowStrandedBegin, windowStrandedEnd);
		logger.debug("UPSTREAM_WINDOW\t" + gene.getName() + "\t" + gene.toUCSC() + ":" + gene.getOrientation().toString() + "\tmin_dist=" + MIN_UPSTREAM_DISTANCE + "\tmax_dist=" + MAX_UPSTREAM_DISTANCE + "\t" + chr.getId() + ":" + windowStart + "-" + windowEnd);
		return NickingGuideRNAPair.findAll(chr, windowStart, windowEnd, gene);
	}
	
	
	

	
	private static void validateMinMaxDist(int minDistance, int maxDistance) {
		/*if(minDistance < 0 || maxDistance < 0) {
			throw new IllegalArgumentException("Min and max distances must be > 0");
		}*/
		if(minDistance >= maxDistance) {
			throw new IllegalArgumentException("Min distance must be < max distance");
		}
	}
	
	/**
	 * @param args
	 * @throws IOException 
	 * @throws InterruptedException 
	 */
	public static void main(String[] args) throws IOException, InterruptedException {
		
		CommandLineParser p = new CommandLineParser();
		
		p.addBooleanArg("--debug", "Debug logging", false, false);
		p.addStringArg("-g", "Genome fasta", true);
		p.addStringArg("-ba", "Full gene annotation in bed format", true);
		p.addStringArg("-bt", "Target genes in bed format", true);
		p.addIntArg("-mind", "Minimum distance downstream of transcription stop for downstream guide RNA pairs", false, MIN_DOWNSTREAM_DISTANCE);
		p.addIntArg("-maxd", "Maximum distance downstream of transcription stop for downstream guide RNA pairs", false, MAX_DOWNSTREAM_DISTANCE);
		p.addIntArg("-minu", "Minimum distance upstream of transcription start for upstream guide RNA pairs", false, MIN_UPSTREAM_DISTANCE);
		p.addIntArg("-maxu", "Maximum distance upstream of transcription start for upstream guide RNA pairs", false, MAX_UPSTREAM_DISTANCE);
		p.addIntArg("-minn", "Minimum distance between guide RNA and nearest gene other than target gene", false, MIN_DIST_TO_NEAREST_GENE);
		p.addIntArg("-maxp", "Max inner distance between paired restriction enzyme cut sites", false, MAX_INNER_DIST_CUT_SITE_PAIRS);
		p.addIntArg("-maxre", "Max distance to restriction enzyme cut site", false, MAX_DIST_TO_RESTRICTION_SITE);
		p.addIntArg("-maxgp", "Max inner distance between guide RNA pairs", false, MAX_INNER_DIST_GUIDE_RNA_PAIRS);
		p.addStringArg("-o", "Output file prefix", true);
		p.addBooleanArg("-dnc", "Enforce double nick configuration for paired guide RNAs", false, ENFORCE_DOUBLE_NICK_CONFIGURATION);
		p.addBooleanArg("-mdg", "Enforce minimum distance to nearest gene", false, ENFORCE_MIN_DIST_TO_NEAREST_GENE);
		p.addBooleanArg("-re", "Enforce proximity of downstream guides to restriction enzyme pair", false, ENFORCE_DOWNSTREAM_PROXIMITY_TO_RESTRICTION_ENZYME);
		p.addBooleanArg("-ot", "Enforce minimum off target score", false, ENFORCE_MIN_OFF_TARGET_SCORE);
		p.addBooleanArg("-ge", "Enforce maximum guide efficacy score", false, ENFORCE_MAX_GUIDE_EFFICACY_SCORE);
		p.addDoubleArg("-mge", "Max guide efficacy score", false, MAX_GUIDE_EFFICACY_SCORE);
		p.addStringArg("-se", "File containing list of restriction enzymes for single cut sites for downstream proximity filter", false, null);
		p.addStringArg("-pe", "File containing list of restriction enzyme pairs for paired cut sites for downstream proximity filter (line format: left_enzyme right_enzyme)", false, null);
		p.addBooleanArg("-fp", "For regions with no guide RNA pairs passing all filters, write all failed pairs to bed file", false, WRITE_FAILED_PAIRS_FOR_MISSING_REGIONS);
		p.addStringArg("-offTargetBits", "File containing bitpacked NGG sites", false, null);
		p.addIntArg("-minOffTargetScore", "Minimum score in the off target analysis to pass", false, MIN_OFF_TARGET_SCORE);
		p.addIntArg("-nbe", "Number of guide RNA pairs to get per cut, sorted by efficacy score", false, NUM_BEST_PAIRS_TO_GET_PER_CUT);
		p.addStringListArg("-fe", "Restriction enzyme whose recognition sequence cannot appear in any guide RNA sequences (repeatable)", false, null);	
		p.addIntArg("-nt", "Number of threads", false, 1);
		p.addBooleanArg("-pno", "Pool guide RNA pairs into as few pools of pairwise nonoverlapping guides as possible", false, POOL_NON_OVERLAPPING);
		p.addStringArg("-lfo", "Left flanking sequence for oligos", true);
		p.addStringArg("-rfo", "Right flanking sequence for oligos", true);
		p.addIntArg("-maxpool", "Max number of pools of pairwise non-overlapping guide pairs", false, MAX_NUM_POOLS);
		p.addStringArg("-p3", "primer3_core executable", true);
		p.addIntArg("-pl", "Primer length", false, NickingGuideRNAPair.PRIMER_LENGTH);
		p.addDoubleArg("-ptm", "Optimal primer Tm", false, NickingGuideRNAPair.OPTIMAL_PRIMER_TM);
		
		p.parse(args, true);
		
		if(p.getBooleanArg("--debug")) {
			RestrictionEnzymeCutSite.logger.setLevel(Level.DEBUG);
			RestrictionEnzymeFactory.logger.setLevel(Level.DEBUG);
			TypeIISRestrictionEnzyme.logger.setLevel(Level.DEBUG);
			GuidePairDoubleNickConfiguration.logger.setLevel(Level.DEBUG);
			DoubleNickCRISPRDesigner.logger.setLevel(Level.DEBUG);
			GuideEfficacyScore.logger.setLevel(Level.DEBUG);
			GuideOffTargetScore.log.setLevel(Level.DEBUG);
			GuideRNA.logger.setLevel(Level.DEBUG);
			NickingGuideRNAPair.logger.setLevel(Level.DEBUG);
			GuideSufficientIsolation.logger.setLevel(Level.DEBUG);
		}
		String genomeFasta = p.getStringArg("-g");
		String annotBed = p.getStringArg("-ba");
		String targetBed = p.getStringArg("-bt");
		PRIMER3_CORE = p.getStringArg("-p3");
		NickingGuideRNAPair.OPTIMAL_PRIMER_TM = p.getDoubleArg("-ptm");
		NickingGuideRNAPair.PRIMER_LENGTH = p.getIntArg("-pl");
		MAX_NUM_POOLS = p.getIntArg("-maxpool");
		LEFT_OLIGO_FLANKING_SEQUENCE = p.getStringArg("-lfo");
		RIGHT_OLIGO_FLANKING_SEQUENCE = p.getStringArg("-rfo");
		MIN_DOWNSTREAM_DISTANCE = p.getIntArg("-mind");
		MAX_DOWNSTREAM_DISTANCE = p.getIntArg("-maxd");
		MIN_UPSTREAM_DISTANCE = p.getIntArg("-minu");
		MAX_UPSTREAM_DISTANCE = p.getIntArg("-maxu");
		MIN_DIST_TO_NEAREST_GENE = p.getIntArg("-minn");
		ENFORCE_DOUBLE_NICK_CONFIGURATION = p.getBooleanArg("-dnc");
		ENFORCE_MIN_DIST_TO_NEAREST_GENE = p.getBooleanArg("-mdg");
		ENFORCE_DOWNSTREAM_PROXIMITY_TO_RESTRICTION_ENZYME = p.getBooleanArg("-re");
		ENFORCE_MIN_OFF_TARGET_SCORE = p.getBooleanArg("-ot");
		ENFORCE_MAX_GUIDE_EFFICACY_SCORE = p.getBooleanArg("-ge");
		MAX_INNER_DIST_CUT_SITE_PAIRS = p.getIntArg("-maxp");
		MAX_DIST_TO_RESTRICTION_SITE = p.getIntArg("-maxre");
		MAX_GUIDE_EFFICACY_SCORE = p.getDoubleArg("-mge");
		WRITE_FAILED_PAIRS_FOR_MISSING_REGIONS = p.getBooleanArg("-fp");
		MAX_INNER_DIST_GUIDE_RNA_PAIRS = p.getIntArg("-maxgp");
		NUM_BEST_PAIRS_TO_GET_PER_CUT = p.getIntArg("-nbe");
		POOL_NON_OVERLAPPING = p.getBooleanArg("-pno");
		OFF_TARGET_BITS = new File(p.getStringArg("-offTargetBits"));
		MIN_OFF_TARGET_SCORE = p.getIntArg("-minOffTargetScore");
		String outPrefix = p.getStringArg("-o");
		String listFileSingleEnzymes = p.getStringArg("-se");
		String listFilePairedEnzymes = p.getStringArg("-pe");
		Collection<String> enzymesToAvoid = p.getStringListArg("-fe");
		int numThreads = p.getIntArg("-nt");
		
		validateMinMaxDist(MIN_DOWNSTREAM_DISTANCE, MAX_DOWNSTREAM_DISTANCE);
		validateMinMaxDist(MIN_UPSTREAM_DISTANCE, MAX_UPSTREAM_DISTANCE);
		
		DoubleNickCRISPRDesigner dncd = new DoubleNickCRISPRDesigner(genomeFasta, targetBed, annotBed);
		
		if(enzymesToAvoid != null) {
			for(String enzymeName : enzymesToAvoid) {
				FORBIDDEN_ENZYMES.add(RestrictionEnzymeFactory.getRestrictionEnzyme(enzymeName));
			}
		}
		if(!FORBIDDEN_ENZYMES.isEmpty()) {
			for(RestrictionEnzyme e : FORBIDDEN_ENZYMES) {
				logger.info("Avoiding guide RNAs that contain recognition sequence for " + e.getName());
			}
		}
		
		if(ENFORCE_DOWNSTREAM_PROXIMITY_TO_RESTRICTION_ENZYME) {
			if(listFileSingleEnzymes == null || listFilePairedEnzymes == null) {
				throw new IllegalArgumentException("To use restriction enzyme filter must provide -se and -pe");
			}
			dncd.useRestrictionEnzymes(listFileSingleEnzymes, listFilePairedEnzymes, outPrefix);
		}

		dncd.findValidPairsAllGenes(numThreads);

		dncd.writeValidGuidePairsAllGenes(outPrefix);
		
		logger.info("");
		logger.info("All done.");
		
	}

}
