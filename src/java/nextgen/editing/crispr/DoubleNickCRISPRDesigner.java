package nextgen.editing.crispr;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;

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
import nextgen.editing.crispr.predicate.GuidePairDoubleNickConfiguration;
import nextgen.editing.crispr.predicate.GuideProximityToNearestRegion;
import nextgen.editing.crispr.predicate.GuideSufficientEfficacy;
import nextgen.editing.crispr.predicate.GuideSufficientIsolation;
import nextgen.editing.crispr.score.GuideEfficacyScore;
import nextgen.editing.crispr.score.GuideOffTargetScore;

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
	private static int MAX_DIST_TO_RESTRICTION_SITE = 100;
	private static boolean ENFORCE_DOUBLE_NICK_CONFIGURATION = true;
	private static boolean ENFORCE_MIN_DIST_TO_NEAREST_GENE = true;
	private static boolean ENFORCE_DOWNSTREAM_PROXIMITY_TO_RESTRICTION_ENZYME = true;
	private static boolean ENFORCE_MIN_OFF_TARGET_SCORE = true;
	private static boolean ENFORCE_MAX_GUIDE_EFFICACY_SCORE = true;
	private static double MAX_GUIDE_EFFICACY_SCORE = 0.6;
	private static boolean WRITE_FAILED_PAIRS_FOR_MISSING_REGIONS = false;
	private FileWriter failedPairBedWriter;
	private FileWriter failedPairTableWriter;
	
	
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
			GuidePairDoubleNickConfiguration dnc = new GuidePairDoubleNickConfiguration();
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
			if(guideEfficacy == null) {
				logger.warn("You provided a guide efficacy object but it is not being used because ENFORCE_MAX_GUIDE_EFFICACY_SCORE is false");
			}
		}
		
		if(ENFORCE_MIN_OFF_TARGET_SCORE) {
			// TODO off-target score filter
		}
		
		logger.debug("PAIR_PASSES_BASIC_FILTERS\t" + pair.toString());
		return true;
		
	}
	
	private void writeValidGuidePairsAllGenes(String outFilePrefix) throws IOException, InterruptedException {
		
		if(WRITE_FAILED_PAIRS_FOR_MISSING_REGIONS) {
			failedPairBedWriter = new FileWriter(outFilePrefix + "_guide_pairs_failing_filters.bed");
			failedPairTableWriter = new FileWriter(outFilePrefix + "_guide_pairs_failing_filters_oligos.out");
			String header = "target_gene\toligo_ID\t" + NickingGuideRNAPair.getOligoFieldNames();
			failedPairTableWriter.write(header + "\n");
		}
		
		// Get the valid pairs
		Collection<NickingGuideRNAPair> downstreamPairs = findValidDownstreamPairsAllGenes();
		Collection<NickingGuideRNAPair> upstreamPairs = findValidUpstreamPairsAllGenes();
		Collection<NickingGuideRNAPair> allPairs = new ArrayList<NickingGuideRNAPair>();
		allPairs.addAll(downstreamPairs);
		allPairs.addAll(upstreamPairs);
		
		// Write as bed file
		String bedFile = outFilePrefix + "_guide_pairs.bed";
		logger.debug("");
		logger.info("Writing all valid guide RNA pairs for all genes to file " + bedFile);
		NickingGuideRNAPair.writeBED(allPairs, bedFile);
		
		// Write as tables
		String oligoTable = outFilePrefix + "_oligos.out";
		logger.info("Writing oligos to " + oligoTable);
		NickingGuideRNAPair.writeOligoTable(allPairs, oligoTable);
		
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
	
	private Collection<NickingGuideRNAPair> findValidDownstreamPairsAllGenes() throws IOException, InterruptedException {
		Collection<NickingGuideRNAPair> rtrn = new ArrayList<NickingGuideRNAPair>();
		for(String chr : targetGenes.keySet()) {
			for(Gene gene : targetGenes.get(chr)) {
				logger.info("Finding downstream guide RNA pairs for gene " + gene.getName());
				rtrn.addAll(findAllValidPairsDownstreamOfTranscriptionStop(gene));
			}
		}
		return rtrn;
	}
	
	private Collection<NickingGuideRNAPair> findValidUpstreamPairsAllGenes() throws IOException, InterruptedException {
		Collection<NickingGuideRNAPair> rtrn = new ArrayList<NickingGuideRNAPair>();
		for(String chr : targetGenes.keySet()) {
			for(Gene gene : targetGenes.get(chr)) {
				logger.debug("");
				logger.info("Finding upstream guide RNA pairs for gene " + gene.getName());
				rtrn.addAll(findAllValidPairsUpstreamOfTranscriptionStart(gene));
			}
		}
		return rtrn;
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
		
		logger.debug("Before filters there are " + allPairs.size() + " pairs downstream of transcription stop.");
		Collection<NickingGuideRNAPair> rtrn = new ArrayList<NickingGuideRNAPair>();
		
		// Collect restriction enzyme sites downstream of the gene
		Collection<Annotation> restrictionSites = new ArrayList<Annotation>();
		if(ENFORCE_DOWNSTREAM_PROXIMITY_TO_RESTRICTION_ENZYME) {
			restrictionSites.addAll(RestrictionEnzymeCutSite.asAnnotations(downstreamRestrictionEnzymeCutSitesByGene.get(gene)));
			restrictionSites.addAll(RestrictionEnzymeCutSitePair.asAnnotations(downstreamRestrictionEnzymeCutSitePairsByGene.get(gene)));
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
				logger.debug("PAIR_PASSES_ALL_DOWNSTREAM_FILTERS\t" + pair.toString());
				rtrn.add(pair);
			}
		}
		if(rtrn.isEmpty()) {
			logger.warn("NO_VALID_DOWNSTREAM_GUIDE_PAIRS\t" + gene.getName());
			
			if(WRITE_FAILED_PAIRS_FOR_MISSING_REGIONS) {
				
				logger.warn("WRITING_FAILED_PAIRS\t" + gene.getName());
				GuidePairDoubleNickConfiguration dnc = new GuidePairDoubleNickConfiguration();
				
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
					
					failedPairBedWriter.write(pair.toBED(bedName) + "\n");
					failedPairTableWriter.write(gene.getName() + "\t" + pair.toString() + "\t" + pair.getOligos() + "\n");

				}
			}
			

		}
		return rtrn;
	}
	
	/**
	 * Find all guide RNA pairs in a window upstream of transcription start that pass all the upstream pair filters
	 * @param gene The gene
	 * @return Collection of valid guide RNA pairs
	 * @throws InterruptedException 
	 * @throws IOException 
	 */
	private Collection<NickingGuideRNAPair> findAllValidPairsUpstreamOfTranscriptionStart(Gene gene) throws IOException, InterruptedException {

		Collection<NickingGuideRNAPair> allPairs = findAllPossibleGuideRNAsUpstreamOfTranscriptionStart(gene);
		if(allPairs.isEmpty()) {
			return new ArrayList<NickingGuideRNAPair>();
		}
		
		GuideSufficientEfficacy ge = null;
		if(ENFORCE_MAX_GUIDE_EFFICACY_SCORE) {
			ge = new GuideSufficientEfficacy(new GuideEfficacyScore(NickingGuideRNAPair.getIndividualGuideRNAs(allPairs)), MAX_GUIDE_EFFICACY_SCORE);
		}

		
		Collection<NickingGuideRNAPair> rtrn = new ArrayList<NickingGuideRNAPair>();
		for(NickingGuideRNAPair pair : allPairs) {
			boolean passes = true;
			if(!guidePairPassesAllBasicFilters(gene, pair, ge)) {
				passes = false;
			}
			
			if(passes) {
				logger.debug("PAIR_PASSES_ALL_UPSTREAM_FILTERS\t" + pair.toString());
				rtrn.add(pair);
			}
		}
		if(rtrn.isEmpty()) {
			logger.warn("NO_VALID_UPSTREAM_GUIDE_PAIRS\t" + gene.getName());
			
			if(WRITE_FAILED_PAIRS_FOR_MISSING_REGIONS) {
				logger.warn("WRITING_FAILED_PAIRS\t" + gene.getName());
				GuidePairDoubleNickConfiguration dnc = new GuidePairDoubleNickConfiguration();
				for(NickingGuideRNAPair pair : allPairs) {
					if(!dnc.evaluate(pair)) {
						continue;
					}
					String bedName = pair.toString() + ":" + getFailureMessageBasicFilters(pair, ge);
					failedPairBedWriter.write(pair.toBED(bedName) + "\n");
					failedPairTableWriter.write(gene.getName() + "\t" + pair.toString() + "\t" + pair.getOligos() + "\n");
				}
			}
			
		}
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
			GuidePairDoubleNickConfiguration dnc = new GuidePairDoubleNickConfiguration();
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
		if(minDistance < 0 || maxDistance < 0) {
			throw new IllegalArgumentException("Min and max distances must be > 0");
		}
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
		
		p.parse(args);
		
		if(p.getBooleanArg("--debug")) {
			RestrictionEnzymeCutSite.logger.setLevel(Level.DEBUG);
			RestrictionEnzymeFactory.logger.setLevel(Level.DEBUG);
			TypeIISRestrictionEnzyme.logger.setLevel(Level.DEBUG);
			GuidePairDoubleNickConfiguration.logger.setLevel(Level.DEBUG);
			DoubleNickCRISPRDesigner.logger.setLevel(Level.DEBUG);
			GuideEfficacyScore.logger.setLevel(Level.DEBUG);
			GuideOffTargetScore.logger.setLevel(Level.DEBUG);
			GuideRNA.logger.setLevel(Level.DEBUG);
			NickingGuideRNAPair.logger.setLevel(Level.DEBUG);
			GuideSufficientIsolation.logger.setLevel(Level.DEBUG);
		}
		String genomeFasta = p.getStringArg("-g");
		String annotBed = p.getStringArg("-ba");
		String targetBed = p.getStringArg("-bt");
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
		String outPrefix = p.getStringArg("-o");
		String listFileSingleEnzymes = p.getStringArg("-se");
		String listFilePairedEnzymes = p.getStringArg("-pe");
		
		validateMinMaxDist(MIN_DOWNSTREAM_DISTANCE, MAX_DOWNSTREAM_DISTANCE);
		validateMinMaxDist(MIN_UPSTREAM_DISTANCE, MAX_UPSTREAM_DISTANCE);
		
		DoubleNickCRISPRDesigner dncd = new DoubleNickCRISPRDesigner(genomeFasta, targetBed, annotBed);
		
		if(ENFORCE_DOWNSTREAM_PROXIMITY_TO_RESTRICTION_ENZYME) {
			if(listFileSingleEnzymes == null || listFilePairedEnzymes == null) {
				throw new IllegalArgumentException("To use restriction enzyme filter must provide -se and -pe");
			}
			dncd.useRestrictionEnzymes(listFileSingleEnzymes, listFilePairedEnzymes, outPrefix);
		}

		dncd.writeValidGuidePairsAllGenes(outPrefix);
		
		logger.info("");
		logger.info("All done.");
		
	}

}
