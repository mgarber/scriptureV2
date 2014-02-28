package nextgen.core.capture.arrayscheme;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.TreeMap;

import org.apache.log4j.Logger;

import nextgen.core.capture.OligoPool;
import nextgen.core.capture.ProbeSet;
import nextgen.core.pipeline.ConfigFile;
import nextgen.core.pipeline.ConfigFileOptionValue;
import broad.core.parser.StringParser;
import broad.core.sequence.Sequence;

public class GroupedStackedPoolScheme implements PoolScheme {
	
	private Map<String, String> sequenceGroups;
	@SuppressWarnings("unused")
	private static Logger logger = Logger.getLogger(GroupedStackedPoolScheme.class.getName());
	private StackedSimplePoolScheme scheme;
	
	public GroupedStackedPoolScheme(ConfigFile configFile) {
		setFromConfigFile(configFile);
	}

	public GroupedStackedPoolScheme() {}

	@Override
	public String name() {
		return "grouped_stacked_pool_scheme";
	}

	@Override
	public String configFileLineDescription() {
		return OligoPool.poolSchemeOptionFlag + "\t" + name() + "\ttable_of_groups(format: seq_name   group_name)";
	}
	
	private static Map<String, String> toMap(String table) throws IOException {
		Map<String, String> rtrn = new TreeMap<String, String>();
		FileReader r = new FileReader(table);
		BufferedReader b = new BufferedReader(r);
		StringParser s = new StringParser();
		while(b.ready()) {
			s.parse(b.readLine());
			if(s.getFieldCount() == 0) continue;
			String seqName = s.asString(0);
			String groupName = s.asString(1);
			if(rtrn.containsKey(seqName)) {
				r.close();
				b.close();
				throw new IllegalArgumentException("Provided key " + seqName + " more than once in table");
			}
			rtrn.put(seqName, groupName);
		}
		r.close();
		b.close();
		return rtrn;
	}
	
	private void createGroups(String tableGroupNames) throws IOException {
		sequenceGroups = toMap(tableGroupNames);
	}

	@Override
	public boolean validConfigFileValue(ConfigFileOptionValue value) {
		// Third field is table with grouping information
		return value.getActualNumValues() == 3 && value.asString(0).equals(OligoPool.poolSchemeOptionFlag) && value.asString(1).equals(name());
	}

	@Override
	public void setParametersFromConfigFile(ConfigFileOptionValue value) {
		throw new UnsupportedOperationException("Method not applicable.");
	}

	@Override
	public Collection<ProbeSet> getProbes(Collection<Sequence> transcripts) {
		Map<String, Collection<Sequence>> seqsByGroup = new TreeMap<String, Collection<Sequence>>();
		for(Sequence seq : transcripts) {
			if(!sequenceGroups.containsKey(seq.getId())) {
				throw new IllegalArgumentException("Group for sequence " + seq.getId() + " was not specified in table.");
			}
			String groupName = sequenceGroups.get(seq.getId());
			if(!seqsByGroup.containsKey(groupName)) {
				seqsByGroup.put(groupName, new ArrayList<Sequence>());
			} 
			seqsByGroup.get(groupName).add(seq);
		}
		Collection<ProbeSet> rtrn = new ArrayList<ProbeSet>();
		for(String group : seqsByGroup.keySet()) {
			Collection<ProbeSet> fromSchemes = scheme.getProbes(seqsByGroup.get(group));
			for(ProbeSet ps : fromSchemes) {
				String layout;
				try {
					layout = ps.getProbes().iterator().next().getProbeLayout().toString();
				} catch (NoSuchElementException e) {
					continue;
				}
				ps.setName(group + "_" + layout);
				rtrn.add(ps);
			}
		}
		return rtrn;
	}

	@Override
	public void setFromConfigFile(ConfigFile file) {
		// Set layouts
		ConfigFileOptionValue poolSchemeVal = file.getSingleValue(OligoPool.arraySchemeSection, OligoPool.poolSchemeOption);
		if(!validConfigFileValue(poolSchemeVal)) {
			throw new IllegalArgumentException("File line not valid:\n" + poolSchemeVal.getFullOptionLine() + "\nlineformat:\n" + configFileLineDescription());
		}
		Collection<ProbeLayout> layouts = OligoPool.getProbeLayoutsFromConfigFile(file);
		Collection<SimplePoolScheme> simpleSchemes = new ArrayList<SimplePoolScheme>();
		for(ProbeLayout layout : layouts) {
			simpleSchemes.add(new SimplePoolScheme(layout));
		}
		scheme = new StackedSimplePoolScheme(simpleSchemes);
		// Set sequence groups
		try {
			createGroups(poolSchemeVal.asString(2));
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}

	}

}
