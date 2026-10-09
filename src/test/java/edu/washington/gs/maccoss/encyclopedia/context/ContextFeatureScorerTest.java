package edu.washington.gs.maccoss.encyclopedia.context;

import java.util.ArrayList;

import junit.framework.TestCase;

public class ContextFeatureScorerTest extends TestCase {
	private static ArrayList<IsolationWindow> list() {
		ArrayList<IsolationWindow> windows=new ArrayList<IsolationWindow>();
		windows.add(new IsolationWindow("GEHPGLSIGDVAK", 427.2262, (byte)3, 0f, 0f, false));
		return windows;
	}

	private static ScoredFeature feature(double mz, int charge) {
		return new ScoredFeature(mz, (byte)charge, false, 1f, 1300f, "-.GEHPGLSIGDVAK.-", "P1", "");
	}

	public void testListedChargeIsAReference() {
		assertNotNull(ContextFeatureScorer.findMatchingMassListWindow(feature(427.2262, 3), list()));
	}

	public void testOtherChargeOfAListedPeptideIsNot() {
		assertNull(ContextFeatureScorer.findMatchingMassListWindow(feature(640.3357, 2), list()));
	}

	public void testUnknownChargeFallsBackToTheSequence() {
		assertNotNull(ContextFeatureScorer.findMatchingMassListWindow(feature(342.0, 0), list()));
	}
}
