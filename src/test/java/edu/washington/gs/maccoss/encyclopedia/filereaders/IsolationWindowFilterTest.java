package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.util.ArrayList;
import java.util.Arrays;

import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import junit.framework.TestCase;

public class IsolationWindowFilterTest extends TestCase {
	private static final Range SURVEY=new Range(450f, 500f);
	private static final Range OWN=new Range(474.5f, 475.5f);       // centred on 475.0
	private static final Range NEIGHBOUR=new Range(474.9f, 475.9f); // centred on 475.4

	private static IsolationWindowFilter.OwnWindows windows(Range... ranges) {
		return new IsolationWindowFilter.OwnWindows(new ArrayList<Range>(Arrays.asList(ranges)));
	}

	public void testPrecursorWithOwnWindowIsScoredOnlyThere() {
		IsolationWindowFilter.OwnWindows w=windows(SURVEY, OWN, NEIGHBOUR);
		assertEquals(OWN, w.get(475.003));
		assertTrue(w.isScoredIn(475.003, OWN));
		assertFalse(w.isScoredIn(475.003, NEIGHBOUR));
		assertFalse(w.isScoredIn(475.003, SURVEY));
	}

	public void testPrecursorWithoutOwnWindowIsScoredWhereverItIs() {
		IsolationWindowFilter.OwnWindows w=windows(SURVEY, OWN, NEIGHBOUR);
		assertNull(w.get(475.2));
		assertTrue(w.isScoredIn(475.2, OWN));
		assertTrue(w.isScoredIn(475.2, NEIGHBOUR));
		assertTrue(w.isScoredIn(475.2, SURVEY));
	}

	public void testNearbyCentreIsNotAnOwnWindow() {
		IsolationWindowFilter.OwnWindows w=windows(SURVEY, OWN, NEIGHBOUR);
		assertNull(w.get(475.03));
		assertTrue(w.isScoredIn(475.03, NEIGHBOUR));
	}

	public void testClosestCentreWinsForNearIsobaricTargets() {
		Range a=new Range(499.5f, 500.5f);
		Range b=new Range(499.51f, 500.51f);
		IsolationWindowFilter.OwnWindows w=windows(a, b);
		assertEquals(a, w.get(500.004));
		assertEquals(b, w.get(500.008));
	}

	public void testEqualCentresPreferTheNarrowerWindow() {
		Range narrow=new Range(600f, 601f);
		Range wide=new Range(599f, 602f);
		assertEquals(narrow, windows(wide, narrow).get(600.5));
	}

	public void testTiledWindowsAreUnaffected() {
		Range low=new Range(400f, 425f);
		Range high=new Range(425f, 450f);
		IsolationWindowFilter.OwnWindows w=windows(low, high);
		assertTrue(w.isScoredIn(412.5, low));
		assertTrue(w.isScoredIn(425.0, low));
		assertTrue(w.isScoredIn(425.0, high));
	}

	private static FragmentScan scan(int index, double lower, double upper) {
		return new FragmentScan("s"+index, "p", index, index, 0, null, lower, upper, new double[0], new float[0], new float[0]);
	}

	public void testNearIsobaricWindowsAreNotMerged() {
		Range own=new Range(715.86243f, 716.86243f);
		ArrayList<FragmentScan> stripes=new ArrayList<FragmentScan>();
		stripes.add(scan(1, own.getStart(), own.getStop()));
		stripes.add(scan(2, 715.87012f, 716.87012f)); // a target 0.008 Th away, scheduled elsewhere
		ArrayList<FragmentScan> kept=IsolationWindowFilter.restrictToWindow(stripes, own);
		assertEquals(1, kept.size());
		assertEquals(1, kept.get(0).getSpectrumIndex());
	}

	public void testPropertyDisablesTheRule() {
		IsolationWindowFilter.OwnWindows w=windows(SURVEY, OWN, NEIGHBOUR);
		System.setProperty("encyclopedia.exactIsolationWindow", "false");
		try {
			assertTrue(w.isScoredIn(475.003, SURVEY));
		} finally {
			System.clearProperty("encyclopedia.exactIsolationWindow");
		}
	}
}
