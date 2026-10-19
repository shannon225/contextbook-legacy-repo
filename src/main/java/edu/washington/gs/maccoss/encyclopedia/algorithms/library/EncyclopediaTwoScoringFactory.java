package edu.washington.gs.maccoss.encyclopedia.algorithms.library;

import java.io.File;
import java.util.ArrayList;
import java.util.Optional;
import java.util.concurrent.BlockingQueue;

import edu.washington.gs.maccoss.encyclopedia.algorithms.AbstractLibraryScoringTask;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PSMPeakScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.PSMScorer;
import edu.washington.gs.maccoss.encyclopedia.algorithms.AbstractScoringResult;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.LibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.filewriters.PeptideScoringResultsConsumer;
import edu.washington.gs.maccoss.encyclopedia.filewriters.ScoringResultsToTSVConsumer;

public class EncyclopediaTwoScoringFactory implements LibraryScoringFactory {
	private final SearchParameters parameters;

	public EncyclopediaTwoScoringFactory(SearchParameters parameters) {
		this.parameters=parameters;
	}
	
	@Override
	public String getPrimaryScoreName() {
		return EncyclopediaTwoScorer.getPrimaryScoreName();
	}
	
	@Override
	public String getName() {
		return "EncyclopeDIA 2.X Scoring System";
	}

	@Override
	public PSMPeakScorer getLibraryScorer(LibraryBackgroundInterface background) {
		return new EncyclopediaTwoScorer(parameters); 
	}

	@Override
	public PeptideScoringResultsConsumer getResultsConsumer(File outputFile, BlockingQueue<AbstractScoringResult> resultsQueue, StripeFileInterface diaFile, LibraryInterface library) {
		return new ScoringResultsToTSVConsumer(outputFile, diaFile, EncyclopediaTwoPointOneAuxillaryPSMScorer.getScoreNames(), resultsQueue, parameters);
	}

	public PeptideScoringResultsConsumer getResultsConsumer(File outputFile, String[] scoreNames, BlockingQueue<AbstractScoringResult> resultsQueue, StripeFileInterface diaFile, LibraryInterface library) {
		return new ScoringResultsToTSVConsumer(outputFile, diaFile, scoreNames, resultsQueue, parameters);
	}

	@Override
	public AbstractLibraryScoringTask getScoringTask(PSMScorer scorer, ArrayList<LibraryEntry> entries, ArrayList<FragmentScan> stripes, Range precursorIsolationRange, float dutyCycle, PrecursorScanMap precursors, BlockingQueue<AbstractScoringResult> resultsQueue) {
		//return new EncyclopediaTwoScoringTask(scorer, entries, stripes, precursorIsolationRange, dutyCycle, precursors, resultsQueue, parameters);
		return new EncyclopediaTwoPointOneScoringTask(scorer, entries, stripes, precursorIsolationRange, dutyCycle, precursors, resultsQueue, parameters);
		//return new IsomerAwareEncyclopediaTwoPointOneScoringTask(scorer, entries, stripes, precursorIsolationRange, dutyCycle, precursors, resultsQueue, parameters);
		
	}

	@Override
	public AbstractLibraryScoringTask getScoringTask(PSMScorer scorer, ArrayList<LibraryEntry> entries, ArrayList<FragmentScan> stripes, Range precursorIsolationRange, float dutyCycle, PrecursorScanMap precursors, BlockingQueue<AbstractScoringResult> resultsQueue, Optional<Range> apexWindow) {
		return new EncyclopediaTwoPointOneScoringTask(scorer, entries, stripes, precursorIsolationRange, dutyCycle, precursors, resultsQueue, parameters, apexWindow);
	}
	
	@Override
	public AbstractLibraryScoringTask getDDAScoringTask(PSMScorer scorer, ArrayList<LibraryEntry> entries, ArrayList<FragmentScan> stripes, PrecursorScanMap precursors, BlockingQueue<AbstractScoringResult> resultsQueue) {
		return new EncyclopediaDDAScoringTask(scorer, entries, stripes, precursors, resultsQueue, parameters);
	}

	@Override
	public SearchParameters getParameters() {
		return parameters;
	}
}
