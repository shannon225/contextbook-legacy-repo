package edu.washington.gs.maccoss.encyclopedia.filereaders;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;
import java.util.Set;

import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;

// Keeps only the scans acquired for the isolation window being scored, and
// scores a precursor that has a window of its own only in that window
public class IsolationWindowFilter {
	// a window's bounds repeat exactly in each of its scans, while the windows of
	// near-isobaric targets can differ by a thousandth of a Th
	private static final double EPSILON=0.0005;
	// scheduled methods centre the window on the target m/z, at most rounded to
	// 0.01; wider would hand background precursors a window by coincidence
	private static final float OWN_WINDOW_TOLERANCE=0.01f;
	private static final String PROPERTY="encyclopedia.exactIsolationWindow";

	public static boolean isEnabled() {
		return !"false".equalsIgnoreCase(System.getProperty(PROPERTY, "true"));
	}

	public static ArrayList<FragmentScan> restrictToWindow(ArrayList<FragmentScan> stripes, Range window) {
		if (!isEnabled()||stripes==null||stripes.isEmpty()||window==null) return stripes;

		if (countDistinctWindows(stripes)<=1) return stripes;

		ArrayList<FragmentScan> kept=new ArrayList<FragmentScan>(stripes.size());
		for (FragmentScan scan : stripes) {
			if (isFromWindow(scan, window)) {
				kept.add(scan);
			}
		}

		if (kept.isEmpty()) {
			Logger.logLine("No scan exactly matches "+window+" m/z; scoring against all "+stripes.size()
					+" overlapping scans.");
			return stripes;
		}

		Logger.logLine("Scheduled acquisition: "+window+" m/z is covered by scans from more than one isolation "
				+"window; keeping the "+kept.size()+" of "+stripes.size()+" acquired for this one. Disable with -D"
				+PROPERTY+"=false.");
		return kept;
	}

	private static int countDistinctWindows(ArrayList<FragmentScan> stripes) {
		Set<String> distinct=new HashSet<String>();
		for (FragmentScan scan : stripes) {
			distinct.add(windowKey(scan.getIsolationWindowLower(), scan.getIsolationWindowUpper()));
			if (distinct.size()>1) return distinct.size();
		}
		return distinct.size();
	}

	private static boolean isFromWindow(FragmentScan scan, Range window) {
		return Math.abs(scan.getIsolationWindowLower()-window.getStart())<EPSILON
				&&Math.abs(scan.getIsolationWindowUpper()-window.getStop())<EPSILON;
	}

	private static String windowKey(double lower, double upper) {
		return Math.round(lower/EPSILON)+"_"+Math.round(upper/EPSILON);
	}

	// A precursor's own window is the one centred on it (closest centre within
	// OWN_WINDOW_TOLERANCE, then the narrowest). Other windows that contain the
	// precursor were acquired for other targets, at other retention times.
	public static class OwnWindows {
		private final Range[] byMiddle;
		private final float[] middles;
		private int skipped=0;

		public OwnWindows(Collection<Range> ranges) {
			byMiddle=ranges.toArray(new Range[0]);
			Arrays.sort(byMiddle, (a, b) -> Float.compare(a.getMiddle(), b.getMiddle()));
			middles=new float[byMiddle.length];
			for (int i=0; i<byMiddle.length; i++) {
				middles[i]=byMiddle[i].getMiddle();
			}
		}

		public Range get(double precursorMz) {
			int i=Arrays.binarySearch(middles, (float)(precursorMz-OWN_WINDOW_TOLERANCE));
			if (i<0) i=-i-1;
			Range best=null;
			float bestOffset=Float.MAX_VALUE;
			for (; i<middles.length&&middles[i]<=precursorMz+OWN_WINDOW_TOLERANCE; i++) {
				Range range=byMiddle[i];
				float offset=(float)Math.abs(middles[i]-precursorMz);
				if (!range.contains(precursorMz)) continue;
				if (best==null||offset<bestOffset||(offset==bestOffset&&range.getRange()<best.getRange())) {
					best=range;
					bestOffset=offset;
				}
			}
			return best;
		}

		// true unless the precursor has an own window and it is not this one
		public boolean isScoredIn(double precursorMz, Range range) {
			if (!isEnabled()) return true;
			Range own=get(precursorMz);
			if (own==null||own.equals(range)) return true;
			synchronized (this) {
				skipped++;
			}
			return false;
		}

		public void log() {
			if (skipped>0) {
				Logger.logLine("Scheduled acquisition: skipped "+skipped+" precursor/window pairs outside the precursor's own "
						+"isolation window. Disable with -D"+PROPERTY+"=false.");
			}
		}
	}
}
