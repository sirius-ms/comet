package de.unijena.bioinf.ms.io.projectspace;

import de.unijena.bioinf.ms.properties.PropertyManager;
import de.unijena.bioinf.sirius.ExperimentResult;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import net.lingala.zip4j.core.ZipFile;
import net.lingala.zip4j.exception.ZipException;
import net.lingala.zip4j.model.ZipParameters;
import net.lingala.zip4j.util.Zip4jConstants;
import org.apache.commons.io.FileUtils;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.text.ParseException;
import java.util.*;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;

/**
 * This is the SiriusProjectSpace. It can only write to uncompressed/
 * directory-based project-spaces. For zip based projects: They will be
 * decompressed in some temp dir and copied/compressed back after the
 * project-space has been successfully closed.
 *
 * @author Markus Fleischauer
 * */
public class SiriusProjectSpace implements ProjectSpace {
    protected static final Logger LOG = LoggerFactory.getLogger(SiriusProjectSpace.class);

    //region Static Helper Methods
    public static @NotNull SiriusProjectSpace create(@Nullable FilenameFormatter filenameFormatter, @NotNull final File projectSpaceRoot, MetaDataSerializer... metaDataSerializers) throws IOException {
        return create(filenameFormatter, projectSpaceRoot, (currentProgress, maxProgress, Message) -> {
        }, metaDataSerializers);
    }

    public static @NotNull SiriusProjectSpace create(@Nullable FilenameFormatter filenameFormatter, @NotNull final File projectSpaceRoot, @NotNull ProgressListener progress, MetaDataSerializer... metaDataSerializers) throws IOException {
        SiriusProjectSpace space = new SiriusProjectSpace(projectSpaceRoot, filenameFormatter, metaDataSerializers);
        space.loadIntoProjectSpace(progress);
        return space;
    }

    public static SiriusProjectSpace create(@NotNull final File rootOutPath, @NotNull final Collection<File> rootInputPaths, @Nullable final FilenameFormatter filenameFormatter, MetaDataSerializer... metaDataSerializers) throws IOException {
        return create(rootOutPath, rootInputPaths, filenameFormatter, (currentProgress, maxProgress, Message) -> {
        }, metaDataSerializers);
    }

    public static SiriusProjectSpace create(@NotNull final File rootOutPath, @NotNull final Collection<File> rootInputPaths, @Nullable final FilenameFormatter filenameFormatter, @NotNull ProgressListener progress, MetaDataSerializer... metaDataSerializers) throws IOException {
        final SiriusProjectSpace merged = new SiriusProjectSpace(rootOutPath, filenameFormatter, metaDataSerializers);
        final TIntSet ids = new TIntHashSet();
        merged.loadIntoProjectSpace(merged.reader, ids, false, progress);//todo correct progress stack
        rootInputPaths.remove(rootOutPath);
        merged.load(progress, ids, rootInputPaths);
        return merged;
    }

    public static void exortToZip(SiriusProjectSpace projectSpace, File zipFile) throws IOException {
        if (!isCompressedProjectSpaceName(zipFile.getName()))
            throw new IllegalArgumentException("Filename needs to be a valid zipped ProjectSpace output.");
        if (zipFile.exists()) {
            if (zipFile.isFile())
                zipFile.delete();
            else throw new IllegalArgumentException("Output path is not a file!");
        }

        projectSpace.copyToZip(zipFile);
    }

    /**
     * Check for a compressed project-space by file ending
     */
    public static boolean isCompressedProjectSpaceName(String fileName) {
        final String lowercaseName = fileName.toLowerCase();
        return lowercaseName.endsWith(".workspace") || lowercaseName.endsWith(".zip") || lowercaseName.endsWith(".sirius");
    }

    /**
     * Just a quick check to discriminate a project-space for an arbitrary folder
     * */
    public static boolean isSiriusWorkspaceDirectory(File f) {
        final File fv = new File(f, "version.txt");
        if (!fv.exists()) return false;
        try (final BufferedReader br = new BufferedReader(new FileReader(fv), 512)) {
            String line = br.readLine();
            if (line == null) return false;
            line = line.toUpperCase();
            if (line.startsWith("SIRIUS")) return true;
            else return false;
        } catch (IOException e) {
            // not critical: if file cannot be read, it is not a valid workspace
            LOG.error(e.getMessage(), e);
            return false;
        }
    }
    //endregion

    //region Internal Fields
    private final LinkedHashMap<String, ExperimentDirectory> experimentIDs = new LinkedHashMap<>();
    protected final List<SummaryWriter> summaryWriters = new ArrayList<>();
    protected final Map<String, String> versionInfo = new HashMap<>();
    protected final AtomicBoolean changed = new AtomicBoolean(false);

    protected DirectoryWriter writer;
    protected DirectoryReader reader;
    protected FilenameFormatter filenameFormatter;

    protected int currentMaxIndex = ExperimentDirectory.NO_INDEX;

    protected File zipRoot = null;
    protected File rootPath;
    //endregion

    //loads existing project-space from reader and uses given writer
    protected SiriusProjectSpace(@NotNull File root, @Nullable FilenameFormatter filenameFormatter, MetaDataSerializer... metaDataSerializers) throws IOException {
        rootPath = root;
        if (rootPath.isFile()) {
            if (isCompressedProjectSpaceName(rootPath.getName().toLowerCase())) {
                zipRoot = rootPath;
                rootPath = de.unijena.bioinf.ChemistryBase.utils.FileUtils.newTempFile(".", "-" + zipRoot.getName() + "-").toFile();
                LOG.info("Zipped workspace found! Unpacking it to temp directory: " + rootPath.getAbsolutePath());
                extractZip();
            }
        } else {
            rootPath.mkdirs();
        }


        if (!rootPath.isDirectory())
            throw new IOException("Project-Space path is not a Directory but has to be: " + rootPath.getAbsolutePath());
        if (!rootPath.canWrite())
            throw new IOException("No writing permission for Project-Space makePath: " + rootPath.getAbsolutePath());

        this.reader = new DirectoryReader(new SiriusFileReader(rootPath), metaDataSerializers);
        this.writer = new DirectoryWriter(new SiriusFileWriter(rootPath), metaDataSerializers);

        this.filenameFormatter = filenameFormatter != null ? filenameFormatter : readFormatter(rootPath);

        addVersionInfo("sirius", PropertyManager.getProperty("de.unijena.bioinf.sirius.versionString", null, "Unknown"));
        registerSummaryWriter(Arrays.stream(metaDataSerializers)
                .filter(v -> v instanceof SummaryWriter).map(v -> (SummaryWriter) v).collect(Collectors.toList()));

    }

    //region Internal Methods
    @NotNull
    private FilenameFormatter readFormatter(File rootPath) {
        FilenameFormatter filenameFormatter = null;
        Path formatFile = rootPath.toPath().resolve(SiriusLocations.SIRIUS_FORMATTER_FILE.fileName());
        if (Files.exists(formatFile)) {
            try {
                String l = Files.readAllLines(formatFile).stream().findFirst().orElse(null);
                if (l != null && !l.isEmpty())
                    try {
                        filenameFormatter = new StandardMSFilenameFormatter(l);
                    } catch (ParseException e) {
                        LOG.warn("Could not parse FileFormatter string: " + l);
                    }
            } catch (IOException ignored) {
                LOG.warn("Could not read FileFormatter file: " + formatFile.toAbsolutePath().toString());
            }
        }
        return filenameFormatter != null ? filenameFormatter : new StandardMSFilenameFormatter();
    }

    private void load(@Nullable final ProgressListener progress, @NotNull final TIntSet ids, @NotNull final Collection<File> toLoad) {
        final List<DirectoryReader.ReadingEnvironment> envs = toLoad.stream().map(file -> {
            try {
                // check for zip stream
                return file.isDirectory() ? new SiriusFileReader(file) : new SiriusZipFileReader(file);
            } catch (IOException e) {
                LOG.error("Cannot read this zipFile: " + file.getPath());
            }
            return null;
        }).filter(Objects::nonNull).collect(Collectors.toList());

        final ProgressListener listener;
        if (progress == null) {
            listener = (a, b, c) -> {
            };
        } else if (envs.size() > 1) {
            final long c = envs.stream().flatMap(env -> env.list().stream().filter(env::isDirectory)).count();
            listener = new StackedProgressListener(3 * (int) c, progress);
        } else {
            listener = progress;
        }

        envs.forEach(env ->
                loadIntoProjectSpace(new DirectoryReader(env), ids, true, listener));
    }

    private ArrayList<ExperimentDirectory> loadIntoProjectSpace(@NotNull ProgressListener listener) {
        experimentIDs.clear();
        return loadIntoProjectSpace(reader, new TIntHashSet(), false, listener);
    }

    private ArrayList<ExperimentDirectory> loadIntoProjectSpace(@NotNull final DirectoryReader r, @NotNull final TIntSet ids, final boolean forceRewrite, @NotNull ProgressListener progress) {
        DirectoryReader.DirectoryReaderIterator it = r.iterator();
        AtomicInteger curProgress = new AtomicInteger(0);
        AtomicInteger maxProgress = new AtomicInteger(it.getMaxPossibleSize() * 3);
        final ArrayList<ExperimentDirectory> expDirs = new ArrayList<>(maxProgress.get());

        progress.doOnProgress(curProgress.get(), maxProgress.get());
        it.forEachRemaining(expDir -> {
            expDirs.add(expDir);
            progress.doOnProgress(curProgress.incrementAndGet(), maxProgress.get());

            //resolve index conflicts and specify rewrite behaviour
            if (expDir.hasIndex() && !ids.add(expDir.getIndex())) {
                expDir.setIndex(ExperimentDirectory.NO_INDEX);
                expDir.setRewrite(true);
            }

            if (forceRewrite) {
                expDir.setRewrite(true);
            } else {
                try {
                    addID(expDir);
                } catch (IOException e) {//should not be possible in practice
                    throw new RuntimeException("Duplicate experiment ID/Name when reading Workspace. This looks like a BUG!", e);
                }
            }

            currentMaxIndex = Math.max(currentMaxIndex, expDir.getIndex());

            if (expDir.hasIndex())
                progress.doOnProgress(curProgress.incrementAndGet(), maxProgress.get());
            if (!expDir.isRewrite())
                progress.doOnProgress(curProgress.incrementAndGet(), maxProgress.get());

        });

        progress.doOnProgress(curProgress.get(), maxProgress.updateAndGet((i) -> (3 * expDirs.size())));

        expDirs.stream().filter(ExperimentDirectory::hasNoIndex).forEach(expDir -> {
            final int index = ++currentMaxIndex;
            if (!ids.add(index))
                throw new IllegalArgumentException("Index is already assigned " + index + ". Assigned indices: " + Arrays.toString(ids.toArray())); //todo special index exception
            expDir.setIndex(index);

            progress.doOnProgress(curProgress.incrementAndGet(), maxProgress.get());
        });

        /*
         * rewrite Experiments if necessary. Either to convert the project space to a new format
         * or to move the Experiments to a new writing location or apply a new naming convention
         */
        expDirs.stream().filter(ExperimentDirectory::isRewrite).forEach(expDir -> {
            try {
                LOG.info("Writing " + expDir.getDirectoryName() + "to update your project-space to the current version");
                final ExperimentResult expResult = r.parseExperiment(expDir);
                if (forceRewrite)
                    expDir.setDirectoryName(null);
                writeExperiment(expResult);
                expDir.setRewrite(false);
            } catch (IOException e) {
                LOG.error("Could not update you project-space to the new format", e);
            } finally {
                progress.doOnProgress(curProgress.incrementAndGet(), maxProgress.get());
            }
        });

        return expDirs;
    }

    private ExperimentDirectory createID(@NotNull ExperimentResult result) {
        final int index = ++currentMaxIndex;
        final ExperimentDirectory dir = new ExperimentDirectory(filenameFormatter.formatName(result, index));
        dir.setIndex(index);
        return dir;
    }

    // return true if id was newly  added to the map
    private boolean addID(ExperimentDirectory expDir) throws IOException {
        ExperimentDirectory old = experimentIDs.putIfAbsent(expDir.getDirectoryName(), expDir);
        if (old != null && old != expDir)
            throw new IOException(
                    "Duplicate Experiment ID/Name: " + old.getDirectoryName() +
                            " Either your naming scheme doe not create unique experiment names " +
                            "or some project merging or renaming results in a naming conflict");
        return old == null;
    }

    private ExperimentDirectory removeID(ExperimentDirectory expDir) {
        return experimentIDs.remove(expDir.getDirectoryName());
    }

    protected List<SummaryWriter> makeBasicSummaries() {
        List<SummaryWriter> sums = new ArrayList<>();

        sums.add((experiments, writer) -> {
            try {
                writer.write(SiriusLocations.SIRIUS_VERSION_FILE.fileName(), w -> {
                    for (Map.Entry<String, String> entry : versionInfo.entrySet()) {
                        if (entry.getValue() != null) {
                            w.write(entry.getKey() + "\t" + entry.getValue());
                            w.write(System.lineSeparator());
                        } else {
                            LOG.error("NO version info for key: " + entry.getKey());
                        }
                    }
                });
            } catch (IOException e) {
                LOG.error("Could not write Version info", e);
            }
        });


        sums.add((experiments, writer) -> {
            final String v = PropertyManager.getProperty("de.unijena.bioinf.sirius.cite");
            if (v != null) {
                try {
                    writer.write(SiriusLocations.SIRIUS_CITATION_FILE.fileName(), w -> w.write(v));
                } catch (IOException e) {
                    LOG.error("Could not write Citation info", e);
                }
            } else {
                LOG.error("NO Citation info found!");
            }
        });

        sums.add((experiments, writer) -> {
            try {
                writer.write(SiriusLocations.SIRIUS_FORMATTER_FILE.fileName(), w -> w.write(filenameFormatter.getFormatExpression()));
            } catch (IOException e) {
                LOG.error("Could not write file formatter info", e);
            }
        });

        return sums;
    }

    private void extractZip() throws IOException {
        try {
            ZipFile zipFile = new ZipFile(zipRoot.getPath());
            if (zipFile.isEncrypted())
                throw new IllegalArgumentException("Encrypted Workspaces are not Supported!");
            zipFile.extractAll(rootPath.getPath());
        } catch (ZipException e) {
            throw new IOException(e);
        }
    }

    private void moveBackToZip() throws IOException {
        Path zipRootNew = de.unijena.bioinf.ChemistryBase.utils.FileUtils.newTempFile(zipRoot.toPath().getParent().toString(), ".", "-" + zipRoot.getName());
        copyToZip(zipRootNew.toFile());
        Files.deleteIfExists(zipRoot.toPath());
        Files.move(zipRootNew, zipRoot.toPath());
        if (rootPath.exists()) FileUtils.deleteDirectory(rootPath);
    }

    private void copyToZip(@NotNull File rootPathFrom, @NotNull File zipRootTo) throws IOException {
        try {
            ZipFile zipFile = new ZipFile(zipRootTo);
            ZipParameters p = new ZipParameters();
            p.setIncludeRootFolder(false);
            p.setCompressionLevel(Zip4jConstants.DEFLATE_LEVEL_ULTRA);
            zipFile.createZipFileFromFolder(rootPathFrom, p, false, 65536);
        } catch (ZipException e) {
            throw new IOException("Error during compression. Your compressed Workspace is incomplete an can be found in: "
                    + zipRootTo.toString() + " The uncompressed version can be found at: " + rootPathFrom.toString(), e);
        }
    }
    //endregion


    //region API Methods
    public synchronized void copyToZip(@NotNull File zipRootTo) throws IOException {
        copyToZip(rootPath, zipRootTo);
    }

    public synchronized void load(@Nullable ProgressListener progress, @NotNull File... toLoad) {
        load(progress, Arrays.asList(toLoad));
    }

    public synchronized void load(@Nullable ProgressListener progress, @NotNull Collection<File> toLoad) {
        load(progress, new TIntHashSet(experimentIDs.values().stream().mapToInt(ExperimentDirectory::getIndex).toArray()), toLoad);
    }

    /**
     * Parses an experiment with the given ID.
     * Override this method to implement caching for the project-space
     *
     * @param id ExperimentDirectory that identifies the experiment to parse.
     *           The id needs to be part of the project-space
     * @return The parsed ExperimentResult annotated with the given id
     *
     * @throws IllegalArgumentException if the project-space does not contain the given id/path
     * @throws IOException            if the given path cannot be parsed from disk.
     *
     */
    @Override
    public synchronized ExperimentResult parseExperiment(ExperimentDirectory id) throws IOException, IllegalArgumentException {
        if (!experimentIDs.containsKey(id.getDirectoryName()))
            throw new IllegalArgumentException("The project-space does not contain the given ID: " + id.getDirectoryName());

        return reader.parseExperiment(id);
    }

    @Override
    public synchronized boolean deleteExperiment(ExperimentDirectory id) throws IOException {
        if (removeID(id) != null) {
            writer.deleteExperiment(id);
            return true;
        }
        return false;

    }

    @Override
    public synchronized void writeExperiment(final @NotNull ExperimentResult result) throws IOException {
        final ExperimentDirectory expDir = result.computeAnnotationIfAbsent(ExperimentDirectory.class, () -> createID(result));
        final String nuName = filenameFormatter.formatName(result, expDir.getIndex());

        if (!nuName.equals(expDir.getDirectoryName())) { //rewrite with new name and delete old
            //check if new name conflicts with old name
            final ExperimentDirectory deleteKey = new ExperimentDirectory(expDir.getDirectoryName());
            expDir.setDirectoryName(nuName);
            addID(expDir);
            writer.writeExperiment(result);
            deleteExperiment(deleteKey);
        } else { //override old
            addID(expDir);
            writer.writeExperiment(result);
        }
        changed.set(true);
    }

    @Override
    public synchronized void writeSummaries(@NotNull final Iterable<ExperimentResult> resultsToSummarize) {
        writeSummaries(resultsToSummarize, (a, b, c) -> {
        });
    }

    public synchronized void writeSummaries(@NotNull final Iterable<ExperimentResult> resultsToSummarize, ProgressListener progress) {
        if (!changed.get())
            return;

        final List<SummaryWriter> basicSummaryWriter = makeBasicSummaries();

        int maxProgress = basicSummaryWriter.size() + summaryWriters.size();
        int currentProgress = 0;
        progress.doOnProgress(currentProgress, maxProgress, "Summarizing Project-Space");

        for (SummaryWriter summaryWriter : summaryWriters) {
            summaryWriter.writeSummary(resultsToSummarize, writer);
            progress.doOnProgress(++currentProgress, maxProgress, summaryWriter.getClass().getSimpleName() + " DONE!");
        }

        for (SummaryWriter summaryWriter : basicSummaryWriter) {
            summaryWriter.writeSummary(resultsToSummarize, writer);
            progress.doOnProgress(++currentProgress, maxProgress, summaryWriter.getClass().getSimpleName() + " DONE!");
        }

        changed.set(false);
    }

    @Override
    public synchronized void registerSummaryWriter(List<SummaryWriter> writers) {
        summaryWriters.addAll(writers);
        writers.forEach(sw -> sw.getVersionInfo().forEach(this::addVersionInfo));
    }

    @Override
    public synchronized boolean removeSummaryWriter(List<SummaryWriter> writers) {
        return summaryWriters.removeAll(writers);
    }

    @Override
    public synchronized int getNumberOfWrittenExperiments() {
        return experimentIDs.size();
    }

    public String addVersionInfo(@NotNull final String key, @NotNull final String value) {
        if (key.contains("\t"))
            throw new IllegalArgumentException("TAB is not allowed in the key");
        return versionInfo.put(key, value.replaceAll("\t", " "));
    }


    @Override
    public void close() throws IOException {
        reader.close();
        writer.close();
        if (zipRoot != null)
            moveBackToZip();
    }


    @NotNull
    @Override
    public Iterator<ExperimentResult> parseExperimentIterator() {
        return new ExperimentIterator();
    }

    @NotNull
    @Override
    public Iterator<ExperimentDirectory> iterator() {
        return new IdIterator();
    }
    //endregion

    //internal classes
    public class IdIterator implements Iterator<ExperimentDirectory> {
        private final Iterator<ExperimentDirectory> baseIter = experimentIDs.values().iterator();
        private ExperimentDirectory current;

        @Override
        public boolean hasNext() {
            return baseIter.hasNext();
        }

        @Override
        public ExperimentDirectory next() {
            if (!hasNext()) throw new NoSuchElementException();
            current = baseIter.next();
            return current;
        }

        @Override
        public void remove() {
            try {
                deleteExperiment(current);
            } catch (IOException e) {
                throw new RuntimeException(e);
            }
        }
    }


    public class ExperimentIterator implements Iterator<ExperimentResult> {
        private final Iterator<ExperimentDirectory> baseIter = experimentIDs.values().iterator();
        private ExperimentDirectory current;

        @Override
        public boolean hasNext() {
            return baseIter.hasNext();
        }

        @Override
        public ExperimentResult next() {
            if (!hasNext()) throw new NoSuchElementException();
            try {
                current = baseIter.next();
                return parseExperiment(current);
            } catch (IOException e) {
                throw new RuntimeException(e);
            }
        }

        @Override
        public void remove() {
            try {
                deleteExperiment(current);
            } catch (IOException e) {
                throw new RuntimeException(e);
            }
        }
    }
}
