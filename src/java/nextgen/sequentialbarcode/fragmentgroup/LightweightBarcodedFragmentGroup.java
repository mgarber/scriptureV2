package nextgen.sequentialbarcode.fragmentgroup;

import java.util.Map;

import nextgen.core.annotation.Annotation;
import nextgen.sequentialbarcode.BarcodeSequence;
import nextgen.sequentialbarcode.BarcodedFragment;

/**
 * Fragment group that doesn't store read IDs, only the number of reads at each location
 * @author prussell
 *
 */
public class LightweightBarcodedFragmentGroup implements FragmentGroup {
	
	private BarcodeSequence barcodes;
	private Map<Annotation, Integer> locations;
	
	@Override
	public void addFragment(BarcodedFragment fragment) {
		Annotation a = fragment.getMappedLocation();
		a.setName(null);
		if(locations.containsKey(a)) {
			locations.put(a, Integer.valueOf(locations.get(a).intValue() + 1));
		} else {
			locations.put(a, Integer.valueOf(1));
		}
	}

	@Override
	public BarcodeSequence getBarcodes() {
		return barcodes;
	}

	@Override
	public String toSamAttributeString() {
		throw new UnsupportedOperationException();
	}
	
	
}
