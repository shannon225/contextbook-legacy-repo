package edu.washington.gs.maccoss.encyclopedia.algorithms.library;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;

import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import junit.framework.TestCase;

public class ApexWindowsTest extends TestCase {
	public void testKeyComparesMassesAtTwoDecimals() {
		assertEquals(ApexWindows.key("AC[+57.0214635]K[+8.0142]", (byte)2), ApexWindows.key("AC[+57.02146]K[+8.014199]", (byte)2));
		assertFalse(ApexWindows.key("ACK", (byte)2).equals(ApexWindows.key("ACK", (byte)3)));
		assertFalse(ApexWindows.key("ACK[+8.0142]", (byte)2).equals(ApexWindows.key("ACK", (byte)2)));
	}

	public void testRead() throws IOException {
		File f=File.createTempFile("apex", ".tsv");
		f.deleteOnExit();
		Files.write(f.toPath(), ("PeptideModSeq\tPrecursorCharge\tMinApexRTInSeconds\tMaxApexRTInSeconds\n"
				+"PEPTIDEK[+8.0142]\t2\t100.5\t120.5\n"
				+"PEPTIDEK\t2\t100.5\t120.5\n").getBytes());
		ApexWindows w=ApexWindows.read(f);
		assertEquals(2, w.size());
		Range r=new Range(100.5f, 120.5f);
		assertEquals(r.getStart(), w.getWindow("PEPTIDEK[+8.0142]", (byte)2).getStart(), 1e-4);
		assertEquals(r.getStop(), w.getWindow("PEPTIDEK", (byte)2).getStop(), 1e-4);
		assertNull(w.getWindow("PEPTIDEK", (byte)3));
	}

	public void testReversedWindowFails() throws IOException {
		File f=File.createTempFile("apex", ".tsv");
		f.deleteOnExit();
		Files.write(f.toPath(), "PeptideModSeq\tPrecursorCharge\tMinApexRTInSeconds\tMaxApexRTInSeconds\nPEPTIDEK\t2\t120\t100\n".getBytes());
		try {
			ApexWindows.read(f);
			fail("a window that ends before it starts must be refused");
		} catch (EncyclopediaException e) {
			// expected
		}
	}
}
