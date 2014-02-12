package nextgen.sequentialbarcode.readlayout;

import java.util.ArrayList;
import java.util.Collection;
import java.util.TreeSet;

import org.apache.log4j.Logger;

/**
 * Static factory methods for read layouts
 * @author prussell
 *
 */
public class ReadLayoutFactory {

	public static Logger logger = Logger.getLogger(ReadLayoutFactory.class.getName());

	/**
	 * Create a read layout for sequential 3D barcoding project with ligations of an odd set and even set of barcodes
	 * @param evenBarcodes Even barcodes
	 * @param oddBarcodes Odd barcodes
	 * @param totalNumBarcodes Total number of barcode ligations
	 * @param rpm RPM sequence
	 * @param readLength Full length of sequencing reads
	 * @param maxMismatchBarcode Max number of mismatches in barcode sequence
	 * @param maxMismatchRpm Max number of mismatches in RPM sequence
	 * @param enforceOddEven Require odd and even barcodes to alternate in read sequences
	 * @return The read layout specified by the parameters
	 */
	public static ReadLayout getRead2LayoutRnaDna3D(Collection<Barcode> evenBarcodes, Collection<Barcode> oddBarcodes, int totalNumBarcodes, String rpm, int readLength, int maxMismatchBarcode, int maxMismatchRpm, boolean enforceOddEven) {
		if(enforceOddEven && totalNumBarcodes % 2 != 0) {
			throw new IllegalArgumentException("Total number of barcodes must be even if enforcing odd/even alternation");
		}
		ArrayList<ReadSequenceElement> eltsSequence = new ArrayList<ReadSequenceElement>();
		Collection<Barcode> allBarcodes = new TreeSet<Barcode>();
		allBarcodes.addAll(oddBarcodes);
		allBarcodes.addAll(evenBarcodes);
		
		if(enforceOddEven) {
			BarcodeSet oddBarcodesSet = new BarcodeSet("odd_barcodes", oddBarcodes, maxMismatchBarcode);
			BarcodeSet evenBarcodesSet = new BarcodeSet("even_barcodes", evenBarcodes, maxMismatchBarcode);
			for(int i = 0; i < totalNumBarcodes; i++) {
				if(i % 2 == 0) {
					eltsSequence.add(oddBarcodesSet);
				} else {
					eltsSequence.add(evenBarcodesSet);
				}
			}
		} else {
			BarcodeSet allBarcodesSet = new BarcodeSet("all_barcodes", allBarcodes, maxMismatchBarcode, true, rpm);
			eltsSequence.add(allBarcodesSet);
		}
		eltsSequence.add(new FixedSequence("rpm", rpm, maxMismatchRpm));
		return new ReadLayout(eltsSequence, readLength);
	}

}
