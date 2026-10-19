package edu.washington.gs.maccoss.encyclopedia.algorithms.library;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Locale;
import java.util.Optional;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;

public class ApexWindows {
	public static final String ARGUMENT="-apexWindows";
	private static final Pattern MASS=Pattern.compile("\\[([+-]?\\d+(?:\\.\\d+)?)\\]");

	private final HashMap<String, Range> windows;

	private ApexWindows(HashMap<String, Range> windows) {
		this.windows=windows;
	}

	public static ApexWindows read(File file) {
		HashMap<String, Range> windows=new HashMap<String, Range>();
		try (BufferedReader reader=new BufferedReader(new FileReader(file))) {
			String header=reader.readLine();
			if (header==null) {
				throw new EncyclopediaException(ARGUMENT+" "+file+" is empty");
			}
			String[] names=header.split("\t", -1);
			int seq=column(names, "PeptideModSeq", file), charge=column(names, "PrecursorCharge", file);
			int min=column(names, "MinApexRTInSeconds", file), max=column(names, "MaxApexRTInSeconds", file);
			String line;
			while ((line=reader.readLine())!=null) {
				if (line.trim().isEmpty()) continue;
				String[] cells=line.split("\t", -1);
				float start=Float.parseFloat(cells[min]), stop=Float.parseFloat(cells[max]);
				if (!(start<=stop)) {
					throw new EncyclopediaException(ARGUMENT+" "+file+": window ends before it starts: "+line);
				}
				windows.put(key(cells[seq].trim(), Byte.parseByte(cells[charge].trim())), new Range(start, stop));
			}
		} catch (IOException e) {
			throw new EncyclopediaException("Error reading "+ARGUMENT+" "+file, e);
		} catch (RuntimeException e) {
			if (e instanceof EncyclopediaException) throw e;
			throw new EncyclopediaException("Error parsing "+ARGUMENT+" "+file+": "+e.getMessage(), e);
		}
		Logger.logLine("Apex windows: read "+windows.size()+" precursors from "+file.getName());
		return new ApexWindows(windows);
	}

	private static int column(String[] names, String name, File file) {
		for (int i=0; i<names.length; i++) {
			if (names[i].trim().equalsIgnoreCase(name)) return i;
		}
		throw new EncyclopediaException(ARGUMENT+" "+file+" has no column "+name);
	}

	static String key(String peptideModSeq, byte charge) {
		Matcher m=MASS.matcher(peptideModSeq);
		StringBuffer sb=new StringBuffer();
		while (m.find()) {
			m.appendReplacement(sb, Matcher.quoteReplacement(String.format(Locale.US, "[%+.2f]", Double.parseDouble(m.group(1)))));
		}
		m.appendTail(sb);
		return sb.toString()+"+"+charge;
	}

	public Optional<Range> get(LibraryEntry entry) {
		return Optional.ofNullable(windows.get(key(entry.getPeptideModSeq(), entry.getPrecursorCharge())));
	}

	Range getWindow(String peptideModSeq, byte charge) {
		return windows.get(key(peptideModSeq, charge));
	}

	public int size() {
		return windows.size();
	}
}
