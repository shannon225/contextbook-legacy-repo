package edu.washington.gs.maccoss.encyclopedia.algorithms.quantitation;

import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Optional;
import java.util.Vector;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadFactory;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.zip.DataFormatException;

import com.google.common.util.concurrent.ThreadFactoryBuilder;

import edu.washington.gs.maccoss.encyclopedia.algorithms.alignment.PeakLocationInferrerInterface;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorPeptide;
import edu.washington.gs.maccoss.encyclopedia.datastructures.FragmentScan;
import edu.washington.gs.maccoss.encyclopedia.datastructures.IntegratedLibraryEntry;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PSMData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.PrecursorScanMap;
import edu.washington.gs.maccoss.encyclopedia.datastructures.Range;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchJobData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.SearchParameters;
import edu.washington.gs.maccoss.encyclopedia.filereaders.IsolationWindowFilter;
import edu.washington.gs.maccoss.encyclopedia.filereaders.LibraryInterface;
import edu.washington.gs.maccoss.encyclopedia.filereaders.StripeFileInterface;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParser;
import edu.washington.gs.maccoss.encyclopedia.utils.io.TableParserMuscle;
import edu.washington.gs.maccoss.encyclopedia.utils.massspec.PeptideUtils;
import edu.washington.gs.maccoss.encyclopedia.utils.threading.ProgressIndicator;
import gnu.trove.map.hash.TObjectFloatHashMap;

public class PeptideQuantExtractor {
	private static final byte NO_ENTRY_CHARGE = 0;

	private final ProgressIndicator progress;
	private final StripeFileInterface stripefile;
	//private final PhosphoLocalizer localizer;
	private final SearchParameters parameters;
	
	public PeptideQuantExtractor(ProgressIndicator progress, StripeFileInterface stripefile, SearchParameters parameters) {
		this.progress=progress;
		this.stripefile=stripefile;
		this.parameters=parameters;
		
		/*PhosphoLocalizer localizer;
		if (parameters.getLocalizingModification().isPresent()&&searchedLibrary!=null) {
			try {
				localizer=new PhosphoLocalizer(stripefile, parameters.getLocalizingModification().get(), searchedLibrary, parameters);
			} catch (DataFormatException dfe) {
				localizer=null;
			} catch (IOException ioe) {
				localizer=null;
			} catch (SQLException sqle) {
				localizer=null;
			}
		} else {
			localizer=null;
		}
		this.localizer=localizer;*/
	}

	public static ArrayList<IntegratedLibraryEntry> parseSearchFeatures(ProgressIndicator progress, final SearchJobData job, boolean limitToQuantifiable, ArrayList<PercolatorPeptide> globalPassingPSMIDs, ArrayList<PercolatorPeptide> localPassingPSMIDs, final Optional<PeakLocationInferrerInterface> inferrer, StripeFileInterface stripeFile, LibraryInterface searchedLibrary, final SearchParameters parameters) {
		HashMap<String, PSMData> uniquedData = findTargetPSMData(job, globalPassingPSMIDs, localPassingPSMIDs, inferrer,
				parameters);
		
		try {
			PeptideQuantExtractor extractor=new PeptideQuantExtractor(progress, stripeFile, parameters);
			ArrayList<IntegratedLibraryEntry> extractPeptides=extractor.extractPeptides(uniquedData.values(), inferrer, limitToQuantifiable);
			
			Logger.logLine("Attempted extraction for: "+uniquedData.size()+", found "+extractPeptides.size());
			return extractPeptides;
		} catch (IOException ioe) {
			Logger.errorLine("Error processing "+stripeFile.getFile().getName());
			throw new EncyclopediaException("Error parsing Stripe file", ioe);
		} catch (SQLException sqle) {
			Logger.errorLine("Error processing "+stripeFile.getFile().getName());
			throw new EncyclopediaException("Error parsing Stripe file", sqle);
		} catch (DataFormatException dfe) {
			Logger.errorLine("Error processing "+stripeFile.getFile().getName());
			throw new EncyclopediaException("Error parsing Stripe file", dfe);
		} catch (InterruptedException ie) {
			Logger.errorLine("Error processing "+stripeFile.getFile().getName());
			throw new EncyclopediaException("Error parsing Stripe file", ie);
		}
	}

	public static HashMap<String, PSMData> findTargetPSMData(final SearchJobData job, ArrayList<PercolatorPeptide> globalPassingPSMIDs, ArrayList<PercolatorPeptide> localPassingPSMIDs, final Optional<PeakLocationInferrerInterface> inferrer, final SearchParameters parameters) {
		HashSet<String> passingPeptideSequences=new HashSet<String>();
		final HashMap<String, PercolatorPeptide> savedPeptides = new HashMap<>();
		for (PercolatorPeptide psm : globalPassingPSMIDs) {
			String psmID=psm.getPsmID();
			boolean isDecoy=PercolatorPeptide.isPSMIDDecoy(psmID);
			
			if (!isDecoy) {
				String peptideModSeq=PeptideUtils.getCorrectedMasses(PercolatorPeptide.getPeptideSequence(psmID), job.getParameters().getAAConstants());
				passingPeptideSequences.add(peptideModSeq);
				savedPeptides.put(peptideModSeq, psm);
			}
		}
		Logger.logLine("Number of global peptides: "+savedPeptides.size()+" vs local peptides: "+localPassingPSMIDs.size());

		final TObjectFloatHashMap<String> localSavedIDs=new TObjectFloatHashMap<String>();
		for (PercolatorPeptide psm : localPassingPSMIDs) {
			String peptideModSeq=PeptideUtils.getCorrectedMasses(PercolatorPeptide.getPeptideSequence(psm.getPsmID()), job.getParameters().getAAConstants());
			if (passingPeptideSequences.contains(peptideModSeq)) {
				localSavedIDs.put(psm.getPsmID(), psm.getPosteriorErrorProb());
			}
		}

		final Vector<PSMData> data=new Vector<PSMData>();
		
		TableParserMuscle muscle=new TableParserMuscle() {
			@Override
			public void processRow(Map<String, String> row) {
				final String psmID=row.get("id");
				final String peptideModSeq=PeptideUtils.getCorrectedMasses(PercolatorPeptide.getPeptideSequence(psmID), job.getParameters().getAAConstants());

				// FIXME need to get peptide charge from window
				final byte precursorCharge=PercolatorPeptide.getCharge(psmID);

				// It shouldn't be possible for the PSM to have zero charge. This assertion won't run unless the JVM was run with -ea
				assert precursorCharge != NO_ENTRY_CHARGE;

				// Only accepted IDs for the correct charge state
				PercolatorPeptide targetPeptide=savedPeptides.get(peptideModSeq);
				if (targetPeptide!=null&&targetPeptide.getPrecursorCharge()==precursorCharge) {
					boolean isDecoy=PercolatorPeptide.isPSMIDDecoy(psmID);
					if (!isDecoy) {
						boolean wasInferred;
						float retentionTime;// in seconds
						int scanID;
						
						// prefer actual identification, fall back on RT inference
						if (localSavedIDs.contains(psmID)) {
							String rtString=row.get("midTime"); // in seconds
							if (rtString!=null) {
								retentionTime=Float.parseFloat(rtString);
								wasInferred=false;
							} else {
								rtString=row.get("RTinMin"); // in minutes so *60
								retentionTime=Float.parseFloat(rtString)*60f;
								wasInferred=false;
							}
							scanID=Integer.parseInt(row.get("ScanNr"));
						} else if (inferrer.isPresent()) {
							if (localSavedIDs.contains(psmID)) {
								String rtString=row.get("midTime"); // in seconds
								if (rtString!=null) {
									retentionTime=Float.parseFloat(rtString);
									wasInferred=false;
								} else {
									rtString=row.get("RTinMin"); // in minutes so *60
									retentionTime=Float.parseFloat(rtString)*60f;
									wasInferred=false;
								}
								
								float warpedRT=inferrer.get().getPreciseRTInSec(job, peptideModSeq, retentionTime);
								if (warpedRT!=retentionTime) {
									Logger.errorLine("Don't trust ID for "+peptideModSeq+" (global RT:"+warpedRT+", local RT:"+retentionTime+"). Using the warped RT!");
									// warping is better (original is way outside the warping margins)
									retentionTime=warpedRT;
									wasInferred=true;
									scanID=-1; // negative scan ID for inferred IDs
								} else {
									// original detection is better (within the warping margins)
									scanID=Integer.parseInt(row.get("ScanNr"));
								}
									
							} else {
								// no detect, so use warped retention time
								retentionTime=inferrer.get().getWarpedRTInSec(job, peptideModSeq);
								wasInferred=true;
								scanID=-1; // negative scan ID for inferred IDs
							}
						} else {
							if (localSavedIDs.contains(psmID)) {
								// no warping, so use RT
								String rtString=row.get("midTime"); // in seconds
								if (rtString!=null) {
									retentionTime=Float.parseFloat(rtString);
									wasInferred=false;
								} else {
									rtString=row.get("RTinMin"); // in minutes so *60
									retentionTime=Float.parseFloat(rtString)*60f;
									wasInferred=false;
								}
								scanID=Integer.parseInt(row.get("ScanNr"));
							} else {
								// not in local search and no warping available
								return;
							}
						}

						double precursorMZ=Double.parseDouble(row.get("precursorMz"));

						float score=localSavedIDs.contains(psmID)?localSavedIDs.get(psmID):1.0f;

						float sortingScore;
						String sortingScoreString=row.get("primary"); // Encyclopedia/XCordia
						if (sortingScoreString==null) {
							sortingScoreString=row.get("xTandem"); // old Encyclopedia
						}
						if (sortingScoreString==null) {
							sortingScoreString=row.get("peakZScore"); // Pecan
						}
						if (sortingScoreString==null) {
							sortingScoreString=row.get("peakBGScore"); // Pecan
						}
						if (sortingScoreString==null) {
							Logger.errorLine("Can't parse score from header from ["+row.keySet()+"]");
							throw new EncyclopediaException("Can't parse score from header from ["+row.keySet()+"]");
						}
						sortingScore=Float.parseFloat(sortingScoreString);

						//String samplingTimeString=row.get("sampledTimes");
						//float duration=samplingTimeString==null?(parameters.getExpectedPeakWidth()):Float.parseFloat(samplingTimeString);

						String proteinString=row.get("protein");
						HashSet<String> accessions=PSMData.stringToAccessions(proteinString);

						savedPeptides.remove(peptideModSeq);
						data.add(new PSMData(accessions, scanID, precursorMZ, precursorCharge, peptideModSeq, retentionTime, score, sortingScore, parameters.getExpectedPeakWidth(), wasInferred, parameters.getAAConstants()));
					}
				}
			}
			
			@Override
			public void cleanup() {
			}
		};
		
		TableParser.parseTSV(job.getPercolatorFiles().getInputTSV(), muscle);

		Logger.logLine("Found additional "+savedPeptides.size()+" globally detected peptides without scores in the individual run. Adding those into the list for re-extraction.");
		for (Entry<String, PercolatorPeptide> entry : savedPeptides.entrySet()) {
			String peptideModSeq=entry.getKey();
			PercolatorPeptide percolatorPeptide=entry.getValue();
			HashSet<String> accessions=PSMData.stringToAccessions(percolatorPeptide.getProteinIDs());
			float retentionTime;
			if (inferrer.isPresent()) {
				retentionTime=inferrer.get().getWarpedRTInSec(job, peptideModSeq);
			} else {
				retentionTime=percolatorPeptide.getRetentionTimeInSec();
			}
			boolean wasInferred=true;
			int scanID=-1; // negative scan ID for inferred IDs
			byte precursorCharge=percolatorPeptide.getPrecursorCharge();
			double precursorMZ=parameters.getAAConstants().getChargedMass(peptideModSeq, precursorCharge);
			data.add(new PSMData(accessions, scanID, precursorMZ, precursorCharge, peptideModSeq, retentionTime, 1.0f, -1.0f, parameters.getExpectedPeakWidth(), wasInferred, parameters.getAAConstants()));
		}
		Logger.logLine("Parsed features and scores for "+data.size()+" peptides.");
		HashMap<String, PSMData> uniquedData=new HashMap<String, PSMData>();
		for (PSMData psmData : data) {
			String key=generateKey(psmData.getPeptideModSeq(), psmData.getPrecursorCharge());
			PSMData prev=uniquedData.get(key);
			if (prev!=null) {
				if (prev.getSortingScore()<psmData.getSortingScore()) {
					// good scores are high
					uniquedData.put(key, psmData);
				}
			} else {
				uniquedData.put(key, psmData);
			}
		}
		return uniquedData;
	}
	
	public static String generateKey(String peptideModSeq, byte precursorCharge) {
		return peptideModSeq+"+"+precursorCharge;
	}
	
	public ArrayList<IntegratedLibraryEntry> extractPeptides(Collection<PSMData> data, final Optional<PeakLocationInferrerInterface> inferrer, boolean limitToQuantifiable) throws IOException, SQLException, DataFormatException, InterruptedException {
		ConcurrentLinkedQueue<IntegratedLibraryEntry> savedEntries=new ConcurrentLinkedQueue<IntegratedLibraryEntry>();
		int cores=parameters.getNumberOfThreadsUsed();
		Logger.logLine("Extracting "+data.size()+" peptides...");
		
		String filename=stripefile.getOriginalFileName();
		
		// get targeted ranges
		ArrayList<Range> ranges=new ArrayList<Range>();
		for (Range range : stripefile.getRanges().keySet()) {
			if (!parameters.useTargetWindowCenter()||range.contains(parameters.getTargetWindowCenter())) {
				ranges.add(range);
			}
		}
		Collections.sort(ranges);
		IsolationWindowFilter.OwnWindows ownWindows=new IsolationWindowFilter.OwnWindows(ranges);
		
		Logger.logLine("Processing precursors scans...");
		PrecursorScanMap precursors=new PrecursorScanMap(stripefile.getPrecursors(-Float.MAX_VALUE, Float.MAX_VALUE));
		
		// get stripes
		int rangesFinished=0;
		float numberOfTasks=2.0f+ranges.size();
		for (Range range : ranges) {
			ArrayList<PSMData> inRange=new ArrayList<PSMData>();
			float minRetentionTime=Float.MAX_VALUE;
			float maxRetentionTime=-Float.MAX_VALUE;
			for (PSMData psm : data) {
				if (range.contains((float)psm.getPrecursorMZ())&&ownWindows.isScoredIn(psm.getPrecursorMZ(), range)) {
					minRetentionTime=Math.min(minRetentionTime, psm.getRetentionTime()-10*psm.getDuration());
					maxRetentionTime=Math.max(maxRetentionTime, psm.getRetentionTime()+10*psm.getDuration());
					inRange.add(psm);
				}
			}
			boolean used=!inRange.isEmpty();

			float baseProgress=(1.0f+rangesFinished)/numberOfTasks;
			String baseMessage="Extracting "+range+" m/z ("+Math.max(0.0f, minRetentionTime/60f)+" to "+Math.max(0.0f, maxRetentionTime/60f)+" min)";
			progress.update(baseMessage, baseProgress);
			if (!used) continue;
			
			Logger.logLine("Quant "+baseMessage);

			ArrayList<FragmentScan> stripes=IsolationWindowFilter.restrictToWindow(
					stripefile.getStripes(range.getMiddle(), minRetentionTime, maxRetentionTime, false), range);
			Collections.sort(stripes);

			// prepare executor for background
			ThreadFactory threadFactory=new ThreadFactoryBuilder().setNameFormat("STRIPE_"+range.getStart()+"to"+range.getStop()+"-%d").setDaemon(true).build();
			LinkedBlockingQueue<Runnable> workQueue=new LinkedBlockingQueue<Runnable>();
			ExecutorService executor=new ThreadPoolExecutor(cores, cores, Long.MAX_VALUE, TimeUnit.NANOSECONDS, workQueue, threadFactory); 


			for (PSMData psm : inRange) {
				executor.submit(new PeptideQuantExtractorTask(filename, psm, inferrer, Optional.ofNullable(null), stripes, Optional.ofNullable(precursors), parameters, savedEntries, limitToQuantifiable));
			}

			executor.shutdown();
			executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
			
			rangesFinished++;
		}
		ownWindows.log();

		ArrayList<IntegratedLibraryEntry> entryList=new ArrayList<IntegratedLibraryEntry>();
		for (IntegratedLibraryEntry entry : savedEntries) {
			entryList.add(entry);
		}

		return entryList;
	}
}
