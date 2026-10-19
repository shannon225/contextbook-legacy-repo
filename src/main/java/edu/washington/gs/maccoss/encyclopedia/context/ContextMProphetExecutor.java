package edu.washington.gs.maccoss.encyclopedia.context;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Optional;

import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.MProphetExecutionData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.ApexWindows;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.MProphetReiter;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.MProphetResult;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.utils.CommandLineParser;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.math.LinearDiscriminantAnalysis;

public class ContextMProphetExecutor {

	private static final float  DEFAULT_FDR  = 0.01f;
	private static final int    DEFAULT_SEED = 1;
	private static final String DEFAULT_PLOTS_DIR = "mprophet_plots";

	public static void main(String[] args) throws IOException {
		HashMap<String, String> arguments = CommandLineParser.parseArguments(args);

		if (args.length == 0
				|| arguments.containsKey("-h")
				|| arguments.containsKey("-help")
				|| arguments.containsKey("--help")) {
			printHelp();
			System.exit(args.length == 0 ? 1 : 0);
		}

		// Feature-file mode: skip the search and run only the train-on-background /
		// apply-to-reference step on feature files that already exist
		if (arguments.containsKey("-background") && arguments.containsKey("-reference")) {
			try {
				File backgroundFeatures = requiredFile(arguments, "-background");
				File referenceFeatures  = requiredFile(arguments, "-reference");
				File fasta              = requiredFile(arguments, "-f");
				float fdr               = Float.parseFloat(arguments.getOrDefault("-fdr", Float.toString(DEFAULT_FDR)));
				int   seed              = Integer.parseInt(arguments.getOrDefault("-seed", Integer.toString(DEFAULT_SEED)));
				File  plotsDir          = plotsDirectory(arguments, referenceFeatures);
				runTrainApply(referenceFeatures, backgroundFeatures, fasta, fdr, seed, plotsDir);
				return;
			} catch (Exception e) {
				Logger.errorLine("ContextMProphetExecutor (feature-file mode) failed: " + e.getMessage());
				e.printStackTrace();
				System.exit(2);
			}
		}

		File dia       = requiredFile(arguments, "-i");
		File library   = requiredFile(arguments, "-l");
		File fasta     = requiredFile(arguments, "-f");
		File massList  = requiredFile(arguments, "-massList");
		float fdr      = Float.parseFloat(arguments.getOrDefault("-fdr",  Float.toString(DEFAULT_FDR)));
		int   seed     = Integer.parseInt(arguments.getOrDefault("-seed", Integer.toString(DEFAULT_SEED)));
		File  plotsDir = plotsDirectory(arguments, dia);
		Optional<File> apexWindows = arguments.containsKey(ApexWindows.ARGUMENT)
				? Optional.of(requiredFile(arguments, ApexWindows.ARGUMENT)) : Optional.empty();

		try {
			// EncyclopeDIA's defaults, overridden by anything the caller passed
			HashMap<String, String> searchArgs = SearchParameterParser.getDefaultParameters();
			searchArgs.putAll(arguments);
			SearchParameters searchParameters = SearchParameterParser.parseParameters(searchArgs);
			Logger.logLine("Fragment tolerance: " + searchParameters.getFragmentTolerance()
					+ ", fragmentation: " + searchParameters.getFragType());

			runEndToEnd(library, fasta, dia, massList, fdr, seed, plotsDir, searchParameters, apexWindows);
		} catch (Exception e) {
			Logger.errorLine("ContextMProphetExecutor failed: " + e.getMessage());
			e.printStackTrace();
			System.exit(2);
		}
	}

	// Full pipeline: score -> split -> train(bg) -> apply(ref)
	public static void runEndToEnd(File library, File fasta, File dia, File massList,
	                               float fdr, int seed, File plotsDir) throws Exception {
		runEndToEnd(library, fasta, dia, massList, fdr, seed, plotsDir,
				SearchParameterParser.getDefaultParametersObject());
	}

	public static void runEndToEnd(File library, File fasta, File dia, File massList,
	                               float fdr, int seed, File plotsDir,
	                               SearchParameters searchParameters) throws Exception {
		runEndToEnd(library, fasta, dia, massList, fdr, seed, plotsDir, searchParameters, Optional.empty());
	}

	public static void runEndToEnd(File library, File fasta, File dia, File massList,
	                               float fdr, int seed, File plotsDir,
	                               SearchParameters searchParameters, Optional<File> apexWindows) throws Exception {

		if (!dia.exists())      throw new IOException("Input file not found: " + dia);
		if (!library.exists())  throw new IOException("Library file not found: " + library);
		if (!fasta.exists())    throw new IOException("FASTA file not found: " + fasta);
		if (!massList.exists()) throw new IOException("Mass list file not found: " + massList);

		// Strip .dia or .mzML (case-insensitive). EncyclopeDIA accepts both.
		String baseName = dia.getAbsolutePath().replaceFirst("(?i)\\.(dia|mzML)$", "");

		Logger.logLine("ContextFeatureScorer: scoring " + dia.getName());
		ContextFeatureScorer.scoreFeatures(library, dia, fasta, baseName, massList.getAbsolutePath(), searchParameters, apexWindows);

		File referenceFeatures  = new File(baseName + "_reference.features.txt");
		File backgroundFeatures = new File(baseName + "_background.features.txt");

		runTrainApply(referenceFeatures, backgroundFeatures, fasta, fdr, seed, plotsDir, searchParameters);
	}

	public static void runTrainApply(File referenceFeatures, File backgroundFeatures,
	                                 File fasta, float fdr, int seed, File plotsDir) throws Exception {
		runTrainApply(referenceFeatures, backgroundFeatures, fasta, fdr, seed, plotsDir,
				SearchParameterParser.getDefaultParametersObject());
	}

	public static void runTrainApply(File referenceFeatures, File backgroundFeatures,
	                                 File fasta, float fdr, int seed, File plotsDir,
	                                 SearchParameters searchParameters) throws Exception {

		SearchParameters params = searchParameters;

		MProphetExecutionData backgroundData = buildExecutionData(backgroundFeatures, fasta, params);
		MProphetExecutionData referenceData  = buildExecutionData(referenceFeatures,  fasta, params);

		final int round = 1;

		Logger.logLine("Training mProphet LDA on background features: " + backgroundFeatures.getName());
		MProphetResult backgroundResult = MProphetReiter.executeMProphetTSV(
				backgroundData, fdr, seed, params.getAAConstants(), round);
		LinearDiscriminantAnalysis backgroundLDA = backgroundResult.getLDA();

		Logger.logLine("Applying background-trained LDA to reference features: " + referenceFeatures.getName());
		MProphetResult referenceResult = MProphetReiter.executeMProphetTSVWithModel(
				referenceData, fdr, backgroundLDA, params.getAAConstants());

		Logger.logLine("Reference passing peptides: " + referenceResult.getPassingPeptides().size());

		if (!plotsDir.exists() && !plotsDir.mkdirs()) {
			Logger.errorLine("Could not create plots directory: " + plotsDir);
		}
		ContextMProphetPlotter.plotContextMProphetResults(
				backgroundData.getPeptideOutputFile(),
				backgroundData.getPeptideDecoyFile(),
				referenceData.getPeptideOutputFile(),
				referenceData.getPeptideDecoyFile(),
				plotsDir);
	}

	private static MProphetExecutionData buildExecutionData(File inputFeatures, File fasta, SearchParameters params) {
		String base = inputFeatures.getAbsolutePath().replaceAll("\\.txt$", "");
		File peptideOutputFile = new File(base + ".pep.output.txt");
		File peptideDecoyFile  = new File(base + ".pep.decoy.txt");
		return new MProphetExecutionData(inputFeatures, fasta, peptideOutputFile, peptideDecoyFile, params);
	}

	// Default the plots dir beside the input
	private static File plotsDirectory(HashMap<String, String> arguments, File anchor) {
		String configured=arguments.get("-plotsdir");
		if (configured!=null) return new File(configured);

		File parent=anchor.getAbsoluteFile().getParentFile();
		return parent==null?new File(DEFAULT_PLOTS_DIR):new File(parent, DEFAULT_PLOTS_DIR);
	}

	private static File requiredFile(HashMap<String, String> args, String flag) {
		String value = args.get(flag);
		if (value == null) {
			Logger.errorLine("Missing required argument: " + flag);
			printHelp();
			System.exit(1);
		}
		return new File(value);
	}

	private static void printHelp() {
		Logger.timelessLogLine("ContextMProphetExecutor");
		Logger.timelessLogLine("Context-mode mProphet: score one PRM file, split by mass-list membership,");
		Logger.timelessLogLine("train LDA on background features, apply to reference features.");
		Logger.timelessLogLine("");
		Logger.timelessLogLine("Feature-file mode (skips the search):");
		Logger.timelessLogLine("  -background <file> background feature TSV");
		Logger.timelessLogLine("  -reference  <file> reference feature TSV");
		Logger.timelessLogLine("  -f          <file> FASTA");
		Logger.timelessLogLine("");
		Logger.timelessLogLine("Required:");
		Logger.timelessLogLine("  -i        <file>   input .dia (or .mzML)");
		Logger.timelessLogLine("  -l        <file>   library (.elib preferred, .dlib accepted)");
		Logger.timelessLogLine("  -f        <file>   FASTA protein database");
		Logger.timelessLogLine("  -massList <file>   assay / mass-list .txt");
		Logger.timelessLogLine("Optional:");
		Logger.timelessLogLine("  -fdr      <float>  peptide FDR threshold (default: " + DEFAULT_FDR + ")");
		Logger.timelessLogLine("  -seed     <int>    random seed for LDA training (default: " + DEFAULT_SEED + ")");
		Logger.timelessLogLine("  -plotsdir <dir>    diagnostic plot directory (default: " + DEFAULT_PLOTS_DIR + "/ beside the input)");
		Logger.timelessLogLine("  " + ApexWindows.ARGUMENT + " <file> where each listed precursor's apex may lie, in seconds (tab-separated:");
		Logger.timelessLogLine("                     PeptideModSeq, PrecursorCharge, MinApexRTInSeconds, MaxApexRTInSeconds); its decoys");
		Logger.timelessLogLine("                     likewise; scores still use every scan of its window");
		Logger.timelessLogLine("");
		Logger.timelessLogLine("Search parameters (passed through to EncyclopeDIA; defaults are ORBITRAP):");
		Logger.timelessLogLine("  -ftol/-ftolunits   fragment tolerance, e.g. -ftol 0.4 -ftolunits AMU");
		Logger.timelessLogLine("  -lftol/-lftolunits library fragment tolerance");
		Logger.timelessLogLine("  -frag              CID | HCD | ...");
		Logger.timelessLogLine("  ION TRAP data (Stellar, LTQ) needs -ftol ~0.4 AMU; the 10 ppm default");
		Logger.timelessLogLine("  yields ZERO identifications on it.");
	}
}
