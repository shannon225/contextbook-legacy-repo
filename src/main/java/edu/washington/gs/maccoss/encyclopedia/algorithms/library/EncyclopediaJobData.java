package edu.washington.gs.maccoss.encyclopedia.algorithms.library;

import edu.washington.gs.maccoss.encyclopedia.ProgramType;
import edu.washington.gs.maccoss.encyclopedia.algorithms.percolator.PercolatorExecutionData;
import edu.washington.gs.maccoss.encyclopedia.datastructures.*;
import edu.washington.gs.maccoss.encyclopedia.filereaders.*;
import edu.washington.gs.maccoss.encyclopedia.utils.EncyclopediaException;
import edu.washington.gs.maccoss.encyclopedia.utils.Logger;
import edu.washington.gs.maccoss.encyclopedia.utils.io.XMLObject;
import edu.washington.gs.maccoss.encyclopedia.utils.io.XMLUtils;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.sql.*;
import java.util.ArrayList;
import java.util.Map;
import java.util.Optional;
import java.util.zip.DataFormatException;

public class EncyclopediaJobData extends QuantitativeSearchJobData implements LibrarySearchJobData, XMLObject {
	public static final String LOG_FILE_SUFFIX=".log";
	public static final String DECOY_PROTEIN_FILE_SUFFIX=".encyclopedia.protein_decoy.txt";
	public static final String OUTPUT_PROTEIN_FILE_SUFFIX=".encyclopedia.protein.txt";
	public static final String DECOY_FILE_SUFFIX=".encyclopedia.decoy.txt";
	public static final String OUTPUT_FILE_SUFFIX=".encyclopedia.txt";
	public static final String FEATURE_FILE_SUFFIX=".features.txt";

	private final LibraryInterface library;
	private final LibraryScoringFactory taskFactory;
	private ApexWindows apexWindows=null;

	public void setApexWindows(ApexWindows apexWindows) {
		this.apexWindows=apexWindows;
	}

	public boolean hasApexWindows() {
		return apexWindows!=null;
	}

	public Optional<Range> getApexWindow(LibraryEntry entry) {
		return apexWindows==null?Optional.empty():apexWindows.get(entry);
	}

	public EncyclopediaJobData(File diaFile, File fastaFile, LibraryInterface library, LibraryScoringFactory taskFactory) {
		this(diaFile, null, getPercolatorExecutionData(diaFile, fastaFile, taskFactory.getParameters()), taskFactory.getParameters(), ProgramType.getGlobalVersion().toString(), library, taskFactory);
	}

	public EncyclopediaJobData(File diaFile, File fastaFile, LibraryInterface library, File outputFile, LibraryScoringFactory taskFactory) {
		this(diaFile, null, getPercolatorExecutionData(outputFile, fastaFile, taskFactory.getParameters()), taskFactory.getParameters(), ProgramType.getGlobalVersion().toString(), library, taskFactory);
	}

	public EncyclopediaJobData(File diaFile, PercolatorExecutionData percolatorFiles, SearchParameters parameters, String version, LibraryInterface library, LibraryScoringFactory taskFactory) {
		this(diaFile, null, percolatorFiles, parameters, version, library, taskFactory);
	}

	public EncyclopediaJobData(File diaFile, StripeFileInterface diaFileReader, PercolatorExecutionData percolatorFiles, SearchParameters parameters, String version, LibraryInterface library, LibraryScoringFactory taskFactory) {
		super(diaFile, diaFileReader, percolatorFiles, parameters, version);

		this.library = library;
		this.taskFactory = taskFactory;
	}

	@Override
	public void writeToXML(Document doc, Element parentElement) {
		Element rootElement=doc.createElement(getClass().getSimpleName());
		parentElement.appendChild(rootElement);

		XMLUtils.writeTag(doc, rootElement, "diaFile", getDiaFile().getAbsolutePath());
		XMLUtils.writeTag(doc, rootElement, "version", getVersion());
		if (library instanceof LibraryFile) {
			XMLUtils.writeTag(doc, rootElement, "library", ((LibraryFile) library).getFile().getAbsolutePath());
		}

		getPercolatorFiles().writeToXML(doc, rootElement);
		getParameters().writeToXML(doc, rootElement);
	}


	public static EncyclopediaJobData readFromXML(Document doc, Element rootElement) {
		if (!rootElement.getTagName().equals(EncyclopediaJobData.class.getSimpleName())) {
			throw new EncyclopediaException("Unexpected XML parsing element, found ["+rootElement.getTagName()+"] when expecting ["+EncyclopediaJobData.class.getSimpleName()+"]");
		}
		File diaFile=null;
		File library=null;
		String version=null;
		PercolatorExecutionData percolatorData=null;
		SearchParameters readParams=null;

		NodeList nodes=rootElement.getChildNodes();

		// read params first
		for (int i = 0; i < nodes.getLength(); i++) {
			Node node = nodes.item(i);
            if (node.getNodeType() == Node.ELEMENT_NODE) {
                Element element = (Element) node;
                if (element.getTagName().equals(SearchParameters.class.getSimpleName())) {
                	readParams=SearchParameters.readFromXML(doc, element);
                }
            }
		}
		if (readParams==null) throw new EncyclopediaException("Found null readParams in "+rootElement.getTagName());

		for (int i = 0; i < nodes.getLength(); i++) {
			Node node = nodes.item(i);
            if (node.getNodeType() == Node.ELEMENT_NODE) {
                Element element = (Element) node;
                if ("diaFile".equals(element.getTagName())) {
                	diaFile=new File(element.getTextContent());
                } else if ("library".equals(element.getTagName())) {
                	library=new File(element.getTextContent());
                } else if ("version".equals(element.getTagName())) {
                	version=element.getTextContent();
                } else if (element.getTagName().equals(PercolatorExecutionData.class.getSimpleName())) {
                	percolatorData=PercolatorExecutionData.readFromXML(doc, element, readParams);
                }
            }
		}

		if (diaFile==null) throw new EncyclopediaException("Found null diaFile in "+rootElement.getTagName());
		if (library==null) throw new EncyclopediaException("Found null library in "+rootElement.getTagName());
		if (version==null) throw new EncyclopediaException("Found null version in "+rootElement.getTagName());
		if (percolatorData==null) throw new EncyclopediaException("Found null percolatorData in "+rootElement.getTagName());

		LibraryInterface libraryObject=BlibToLibraryConverter.getFile(library, percolatorData.getFastaFile(), readParams);

		LibraryScoringFactory factory=EncyclopediaScoringFactory.getDefaultScoringFactory(readParams);
		return new EncyclopediaJobData(diaFile,  percolatorData, readParams,  version, libraryObject, factory);
	}

	@Override
	public SearchJobData updateQuantFile(File f) {
		return new EncyclopediaJobData(f, getPercolatorFiles(), getParameters(), getVersion(), getLibrary(), getTaskFactory());
	}

	protected static PercolatorExecutionData getPercolatorExecutionData(File referenceFileLocation, File fastaFile, SearchParameters parameters) {
		return new PercolatorExecutionData(new File(getPrefixFromOutput(referenceFileLocation) + FEATURE_FILE_SUFFIX), fastaFile,
				new File(getPrefixFromOutput(referenceFileLocation) + OUTPUT_FILE_SUFFIX), new File(getPrefixFromOutput(referenceFileLocation) + DECOY_FILE_SUFFIX), 
				new File(getPrefixFromOutput(referenceFileLocation) + OUTPUT_PROTEIN_FILE_SUFFIX), new File(getPrefixFromOutput(referenceFileLocation) + DECOY_PROTEIN_FILE_SUFFIX), parameters);
	}

	static String getPrefixFromOutput(File outputFile) {
		final String absolutePath = outputFile.getAbsolutePath();

		if (absolutePath.endsWith(OUTPUT_FILE_SUFFIX)) {
			return absolutePath.substring(0, absolutePath.length() - OUTPUT_FILE_SUFFIX.length());
		} else {
			return absolutePath;
		}
	}

	public EncyclopediaJobData updateTaskFactory(LibraryScoringFactory taskFactory) {
		return new EncyclopediaJobData(getDiaFile(), diaFileReader, getPercolatorFiles(), getParameters(), getVersion(), getLibrary(), taskFactory);
	}

	public LibraryInterface getLibrary() {
		return library;
	}

	public LibraryScoringFactory getTaskFactory() {
		return taskFactory;
	}

	public File getResultLibrary() {
		String absolutePath = getPrefixFromOutput(getPercolatorFiles().getPeptideOutputFile());
		return new File(absolutePath + LibraryFile.ELIB);
	}

	@Override
	public String getSearchType() {
		return "EncyclopeDIA";
	}
	
	@Override
	public String getPrimaryScoreName() {
		return taskFactory.getPrimaryScoreName();
	}

	/**
	 * Return an {@code EncyclopeDIAJobData} instance for {@code diaFile} that will
	 * allow processing results for jobs even if the .DIA file doesn't exist. This allows
	 * SearchToBLIB's {@code -alignOnly} option to run without access to raw data, which
	 * is useful e.g. for large experiments where collecting all the raw files is costly.
	 * <p>
	 * Instances will return {@code true} from {@link #hasBeenRun()} even if the .DIA
	 * doesn't exist. It will also return a {@code StripeFileInterface} that provides the
	 * appropriate original file name and TIC (read from the file's results ELIB) but
	 * otherwise throws on attempts to read the .DIA's contents.
	 */
	public static EncyclopediaJobData getDummyFor(File diaFile, File fastaFile, LibraryInterface library, LibraryScoringFactory factory) {
		return new DummyEncyclopediaJobData(diaFile, fastaFile, library, factory);
	}

	/**
	 * Special class overriding key methods to allow processing results for jobs
	 * even if the .DIA file doesn't exist. This allows SearchToBLIB's {@code -alignOnly}
	 * option to run without access to raw data, which is useful e.g. for large
	 * experiments where collecting all the raw files is costly.
	 * <p>
	 * Instances will return {@code true} from {@link #hasBeenRun()} even if the .DIA
	 * doesn't exist. It will also return a {@code StripeFileInterface} that provides the
	 * appropriate original file name and TIC (read from the file's results ELIB) but
	 * otherwise throws on attempts to read the .DIA's contents.
	 */
	private static class DummyEncyclopediaJobData extends EncyclopediaJobData {
		private String originalFileName = null;
		private Float tic = null;

		private DummyEncyclopediaJobData(File diaFile, File fastaFile, LibraryInterface library, LibraryScoringFactory factory) {
			super(diaFile, fastaFile, library, factory);

			if (Files.exists(diaFile.toPath())) {
				Logger.errorLine("Creating a dummy job datum for " + diaFile.getName() + " even though it exists!");
			}
		}

		@Override
		public boolean hasBeenRun() {
			final PercolatorExecutionData percolatorFiles = getPercolatorFiles();
			if (!percolatorFiles.getInputTSV().exists()) {
				Logger.errorLine("Missing feature file: " + percolatorFiles.getInputTSV().getName());
				return false;
			}
			if (!percolatorFiles.getPeptideOutputFile().exists()) {
				Logger.errorLine("Missing output file: " + percolatorFiles.getPeptideOutputFile().getName());
				return false;
			}
			if (!getResultLibrary().exists()) {
				Logger.errorLine("Missing output library: " + getResultLibrary().getName());
				return false;
			}
			return true;
		}

		@Override
		public StripeFileInterface getDiaFileReader() {
			return new StripeFileInterface() {
				@Override
				public Map<Range, WindowData> getRanges() {
					throw new UnsupportedOperationException("File not found: " + getDiaFile().getAbsolutePath());
				}

				@Override
				public Map<String, String> getMetadata() throws IOException, SQLException {
					throw new UnsupportedOperationException("File not found: " + getDiaFile().getAbsolutePath());
				}

				@Override
				public void openFile(File userFile) throws IOException, SQLException {
					throw new UnsupportedOperationException();
				}

				@Override
				public ArrayList<PrecursorScan> getPrecursors(float minRT, float maxRT) throws IOException, SQLException, DataFormatException {
					throw new UnsupportedOperationException("File not found: " + getDiaFile().getAbsolutePath());
				}

				@Override
				public ArrayList<FragmentScan> getStripes(double targetMz, float minRT, float maxRT, boolean sqrt) throws IOException, SQLException {
					throw new UnsupportedOperationException("File not found: " + getDiaFile().getAbsolutePath());
				}

				@Override
				public ArrayList<FragmentScan> getStripes(Range targetMzRange, float minRT, float maxRT, boolean sqrt) throws IOException, SQLException {
					throw new UnsupportedOperationException("File not found: " + getDiaFile().getAbsolutePath());
				}

				@Override
				public float getGradientLength() throws IOException, SQLException {
					throw new UnsupportedOperationException("File not found: " + getDiaFile().getAbsolutePath());
				}

				@Override
				public void close() {
					// no-op
				}

				@Override
				public boolean isOpen() {
					return false;
				}

				@Override
				public File getFile() {
					return getDiaFile();
				}

				@Override
				public String getOriginalFileName() {
					if (null != originalFileName) {
						return originalFileName;
					}

					synchronized (DummyEncyclopediaJobData.this) {
						SET_NAME: if (null == originalFileName) {
							// Workaround: if the DIA file is missing we can't read the
							// original file name, which likely has an .mzML extension.
							// Instead, get the name used in this job's results ELIB.
							try (Connection c = new SQLFile() {}.getConnection(getResultLibrary())) {
								try (Statement s = c.createStatement()) {
									try (ResultSet rs = s.executeQuery(
											"SELECT Value" +
													" FROM Metadata" +
													" WHERE Key = 'sourcefile'" +
													" LIMIT 1;"
									)) {
										if (rs.next()) {
											originalFileName = rs.getString(1);
											break SET_NAME;
										} else {
											throw new SQLException("No `sourcefile` in results ELIB's `metadata`!");
										}
									}
								}
							} catch (IOException | SQLException e) {
								Logger.errorLine("Unable to read from results ELIB for job " + getDiaFile().getName());
								Logger.errorException(e);
							}

							originalFileName = getDiaFile().getName();
						}
					}

					return originalFileName;
				}

				@Override
				public float getTIC() throws IOException, SQLException {
					if (null != tic) {
						return tic;
					}

					synchronized (this) {
						SET_TIC: if (null == tic) {
							// TIC is saved in this file's results ELIB; read directly to avoid copying the ELIB to temp
							try (Connection c = new SQLFile() {}.getConnection(getResultLibrary())) {
								try (PreparedStatement ps = c.prepareStatement(
										"SELECT Value" +
										" FROM Metadata" +
										" WHERE KEY=?" +
										" LIMIT 1;"
								)) {
									ps.setString(1, LibraryFile.SOURCEFILE_TIC_PREFIX + getOriginalFileName());
									try (ResultSet rs = ps.executeQuery()) {
										if (rs.next()) {
											tic = Float.parseFloat(rs.getString(1));
											break SET_TIC;
										} else {
											throw new SQLException("No recorded TIC in results ELIB!");
										}
									}
								}
							} catch (IOException | SQLException e) {
								Logger.errorLine("Unable to read from results ELIB for job " + getOriginalFileName());
								Logger.errorException(e);
							}

							tic = 0.0f;
						}
					}

					return tic;
				}
			};
		}
	}
}