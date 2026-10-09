package edu.washington.gs.maccoss.encyclopedia.context;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.zip.DataFormatException;

import edu.washington.gs.maccoss.encyclopedia.EncyclopediaTwo;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.EncyclopediaTwoJobData;
import edu.washington.gs.maccoss.encyclopedia.algorithms.library.LibraryScoringFactory;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.BlibToLibraryConverter;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface;
import edu.washington.gs.maccoss.encyclopedia.filereaders.SearchParameterParser;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.EmptyProgressIndicator;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.ProgressIndicator;

public class ContextFeatureScorer {

	public static void main(String[] args) throws IOException, SQLException, InterruptedException, DataFormatException {
		
		if (args.length != 4) {
			System.err.println("Usage: ");
			System.err.println("java edu.washington.gs.maccoss.encyclopedia.context.ContextFeatureScorer " + 
			"<rawFilePath> <libraryFilePath> <fastaPath> <massListPath");
			}
		
		String rawFilePath = args[0];
		String libraryFilePath = args[1];
		String fastaPath = args[2];
		String massListPath = args[3];
		
		// Inputs for Search

//		String rawFilePath = "C:/Users/m334793/Documents/targeted_bootstrapper_eval_20260522/varying_number_of_peptides/100_pep/IL2A_GPFDIA_0combined_masked0_assay.dia"; // Raw file to search
//		String libraryFilePath = "C:/Users/m334793/Documents/targeted_bootstrapper_eval_20260522/varying_number_of_peptides/100_pep/IL2_and_IL15_Combo.elib"; // Library to

//		String fastaPath = "C:/Users/m334793/Documents/targeted_bootstrapper_eval_20260522/varying_number_of_peptides/100_pep/mus_musculus_reviewed_uniprot.fasta"; // fasta file for serach
		String baseName = rawFilePath.replaceFirst("\\.dia$", "");

		// Mass list file
//		String massListPath = "C:/Users/m334793/Documents/targeted_bootstrapper_eval_20260522/varying_number_of_peptides/100_pep/IL2A_GPFDIA_0combined_masked0_assay.txt";
//		String rawFilePath = "C:/Users/m334793/Documents/Library/for_context_50perCycle/IT_100ngCurve_100p.dia"; // Raw file to search
//		String libraryFilePath = "C:/Users/m334793/Documents/Library/for_context_50perCycle/IL2_and_IL15_Combo.elib"; // Library to

		final File fasta = new File(fastaPath);
		File rawFile = new File(rawFilePath);
		File library = new File(libraryFilePath);

		try {
			ArrayList<ScoredFeature> partitionedFeatures = scoreFeatures(library, rawFile, fasta, baseName, massListPath);

		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	// a reference is the listed precursor: sequence and charge (0 = charge unknown)
	static IsolationWindow findMatchingMassListWindow(ScoredFeature feature, ArrayList<IsolationWindow> targetWindows) {
		String sequence = cleanPeptideSequence(feature.getSequence());

		for (IsolationWindow window : targetWindows) {
			String compound = cleanPeptideSequence(window.getCompound());

			if (compound.equals(sequence)
					&& (feature.getCharge() == 0 || window.getCharge() == feature.getCharge())) {
				return window;
			}
		}

		return null;
	}
	
//	private static String cleanPeptideSequence(String sequence) {
//		if (sequence == null) return "";
//		sequence = sequence.trim();
//		String[] parts = sequence.split("\\.");
//		if (parts.length == 3) {
//			return parts[1].trim();
//		}
//		return sequence;
//	}

	private static String cleanPeptideSequence(String sequence) {
		if (sequence == null) return "";

		return sequence
				.trim()
				.replaceFirst("^[A-Za-z-]?\\.", "")
				.replaceFirst("\\.[A-Za-z-]?$", "");
	}
	
	private static void writeScoredFeatures(File outputFile, ArrayList<ScoredFeature> features, String header)
			throws IOException {
		try (BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile))) {
			writer.write(header);
			writer.newLine();

			for (ScoredFeature feature : features) {
				writer.write(feature.getOriginalLine());
				writer.newLine();
			}
		}
	}

	private static byte parseCharge(String[] columns) {
		if (Integer.parseInt(columns[23])==1) return 1;
		if (Integer.parseInt(columns[24])==1) return 2;
		if (Integer.parseInt(columns[25])==1) return 3;
		if (Integer.parseInt(columns[26])==1) return 4;

		return 0;
	}



	// Uses EncyclopeDIA's default search parameters: 10 ppm fragment tolerance, CID
	public static ArrayList<ScoredFeature> scoreFeatures(File library, File rawFile, File fasta, String baseName, String massListPath) throws IOException, SQLException, DataFormatException, InterruptedException {
		return scoreFeatures(library, rawFile, fasta, baseName, massListPath,
				SearchParameterParser.getDefaultParametersObject());
	}

	// Score with user-supplied search parameters
	public static ArrayList<ScoredFeature> scoreFeatures(File library, File rawFile, File fasta, String baseName, String massListPath, SearchParameters searchParameters) throws IOException, SQLException, DataFormatException, InterruptedException {

		// Run an Encyclopedia job
		SearchParameters params = searchParameters;
		LibraryScoringFactory scoringForLibrary = EncyclopediaScoringFactory.getDefaultScoringFactory(params);
		LibraryInterface interfaceForLibrary = BlibToLibraryConverter.getFile(library, fasta, params);

		EncyclopediaTwoJobData job = new EncyclopediaTwoJobData(rawFile, fasta, interfaceForLibrary,
				interfaceForLibrary, rawFile, scoringForLibrary);

		ProgressIndicator progress = new EmptyProgressIndicator(true);
		StripeFileInterface interfaceForStripeFile = job.getDiaFileReader();

		// Run Encyclopedia job to get the feature file
		File featuresToSplit = job.getPercolatorFiles().getInputTSV();

		if (featuresToSplit.exists() && featuresToSplit.canRead()) {
			System.out.println("Feature file already exists, skipping feature calculation!");
			System.out.println(featuresToSplit.getAbsolutePath());
		} else {
			System.out.println("Calculating features...");
			EncyclopediaTwo.generateFeatureFile(progress, interfaceForLibrary, job,
					interfaceForStripeFile, java.util.Optional.empty());
		}

		ArrayList<ScoredFeature> uniqueFeatures = new ArrayList<>();
		ArrayList<ScoredFeature> uniqueFeaturesList = uniqueFeatures;
		HashMap<String, ScoredFeature> bestFeatureByPrecursor = new HashMap<>();
		String header;

		// Read all rows the feature file
		try (BufferedReader br = new BufferedReader(new FileReader(featuresToSplit))) {
			header = br.readLine();
			if (header == null) {
				throw new IOException("Feature file is empty, so no header could be read.");
			}

			String line;
			while ((line = br.readLine()) != null) {
				String columns[] = line.split("\t", -1);

				double mz = Double.parseDouble(columns[27]);
				byte featureCharge = parseCharge(columns);
				boolean isDecoy = Integer.parseInt(columns[1]) == -1;
				float primary = Float.parseFloat(columns[3]);
				float retentionTime = Float.parseFloat(columns[29]);
				String sequence = columns[30];
				String protein = columns[31];

				ScoredFeature feature = new ScoredFeature(mz, featureCharge, isDecoy, primary, retentionTime, sequence, protein, line);

				uniqueFeaturesList.add(feature);

				// best feature per precursor, so a charge state that was not
				// targeted cannot stand in for the one that was
				String precursor = sequence + "+" + featureCharge;
				ScoredFeature currentBest = bestFeatureByPrecursor.get(precursor);

				if (currentBest == null || feature.getPrimary() > currentBest.getPrimary()) {
					bestFeatureByPrecursor.put(precursor, feature);
				}

			}
		}
		ArrayList<ScoredFeature> bestFeatures = new ArrayList<>(bestFeatureByPrecursor.values());
		bestFeatures.sort(Comparator.comparing(ScoredFeature::getPrimary).reversed());

		System.out.println("Unique scored features have been found " + bestFeatures.size());


		// Output Paths
		String referenceOutputPath = baseName + "_reference.features.txt";
		String backgroundOutputPath = baseName + "_background.features.txt";
		
		// Output Files
		File referenceOutput = new File(referenceOutputPath);
		File backgroundOutput = new File(backgroundOutputPath);

		// Target mass list
		ArrayList<IsolationWindow> targetWindows = IsolationWindowReader.parseMassList(massListPath);
		System.out.println(targetWindows.size() + " windows cataloged from the mass list");

		ArrayList<ScoredFeature> referenceFeatures = new ArrayList<>();
		ArrayList<ScoredFeature> backgroundFeatures = new ArrayList<>();

		ArrayList<ScoredFeature> partitionedFeatures = new ArrayList<>(); // so that the return is all of the features

		for (ScoredFeature feature : bestFeatures) {
			IsolationWindow matchingWindow = findMatchingMassListWindow(feature, targetWindows);
			
			boolean isOnMassList = matchingWindow != null;
			boolean isMassListDecoy = isOnMassList && matchingWindow.isDecoy();
			boolean isBackground = !isOnMassList;

			ScoredFeature annotatedFeature = new ScoredFeature(feature.getMz(), feature.isDecoy(), feature.getPrimary(), feature.getRetentionTime(), cleanPeptideSequence(feature.getSequence()), feature.getProtein(), feature.getOriginalLine(), isBackground);
			partitionedFeatures.add(annotatedFeature);
			
			System.out.println("Features are being read, sequence for this feature is " + cleanPeptideSequence(feature.getSequence()));
			
			if (isOnMassList) {
				referenceFeatures.add(feature); 
			} else {
				backgroundFeatures.add(feature);
			}
//			if (!feature.isDecoy() && isOnMassList) {
//				referenceFeatures.add(feature);
//			} else if (!feature.isDecoy() && !isOnMassList) {
//				backgroundFeatures.add(feature);
//			} else if (feature.isDecoy() && isOnMassList) {
//				referenceFeatures.add(feature);
//			} else {
//				backgroundFeatures.add(feature);
//			}
		}
		writeScoredFeatures(referenceOutput, referenceFeatures, header);
		writeScoredFeatures(backgroundOutput, backgroundFeatures, header);

		System.out.println("Reference target features: " + referenceFeatures.size());
//		System.out.println("Background target features: " + backgroundFeatures.size());
//		System.out.println("Reference decoy features: " + referenceDecoyFeatures.size());
//		System.out.println("Background decoy features: " + backgroundDecoyFeatures.size());
	return partitionedFeatures;
}
}
