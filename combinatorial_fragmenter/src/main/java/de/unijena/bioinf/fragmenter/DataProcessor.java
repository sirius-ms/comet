package de.unijena.bioinf.fragmenter;

import de.unijena.bioinf.ChemistryBase.chem.MolecularFormula;
import de.unijena.bioinf.ChemistryBase.chem.utils.UnknownElementException;
import de.unijena.bioinf.ChemistryBase.ms.ft.FTree;
import de.unijena.bioinf.babelms.MsIO;
import gnu.trove.map.hash.TObjectIntHashMap;
import gurobi.GRBException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.StandardOpenOption;
import java.util.*;
import java.util.concurrent.*;
import java.util.stream.Collectors;

public class DataProcessor {

    private File spectraDir;
    private File predictionsDir;

    private final File fTreeDir;
    private final File outputDir;
    private final String[] fileNames;

    private static final long SHUFFLE_SEED = 42;

    public DataProcessor(File spectraDir, File predictionsDir, File fTreeDir, File outputDir){
        this(spectraDir, predictionsDir, fTreeDir, outputDir, 1, 0);
    }

    public DataProcessor(File spectraDir, File predictionsDir, File fTreeDir, File outputDir, int numPartitions, int idxPartition){
        this(spectraDir, predictionsDir, fTreeDir, outputDir,
                Arrays.stream(Objects.requireNonNull(outputDir.list())).
                        map(processedInstanceFileName -> processedInstanceFileName.replaceFirst("\\..+$", "")).
                        collect(Collectors.toList()),
                numPartitions, idxPartition);
    }

    public DataProcessor(File spectraDir, File predictionsDir, File fTreeDir, File outputDir, Collection<String> processedInstances){
        this(spectraDir, predictionsDir, fTreeDir, outputDir, processedInstances, 1, 0);
    }

    public DataProcessor(File spectraDir, File predictionsDir, File fTreeDir, File outputDir, Collection<String> processedInstances, int numPartitions, int idxPartition){
        this(fTreeDir, outputDir, processedInstances, numPartitions, idxPartition);

        if(spectraDir != null){
            if(spectraDir.isDirectory()){
                this.spectraDir = spectraDir;
            }else{
                throw new RuntimeException("The given abstract path name denoting the directory containing the mass spectra " +
                        "does not exist or is not a directory.");
            }
        }
        if(predictionsDir != null){
            if(predictionsDir.isDirectory()){
                this.predictionsDir = predictionsDir;
            }else{
                throw new RuntimeException("The given abstract path name denoting the directory containing the predicted compounds" +
                        "does not exist or is not a directory.");
            }
        }
    }

    private DataProcessor(File fTreeDir, File outputDir, Collection<String> processedInstances, int numPartitions, int idxPartition) {
        if(fTreeDir.isDirectory() && outputDir.isDirectory()) {
            this.fTreeDir = fTreeDir;
            this.outputDir = outputDir;

            System.out.println("Filter out already processed instances.");
            List<String> filteredInstanceFileNames = this.filterOutProcessedInstances(processedInstances);
            System.out.println(filteredInstanceFileNames.size() + " instances remain after filtering.");

            System.out.println("Sort the remaining instances lexicographically and shuffle them with " +
                    "seed " + SHUFFLE_SEED + ".");
            Collections.sort(filteredInstanceFileNames);
            Collections.shuffle(filteredInstanceFileNames, new Random(SHUFFLE_SEED));

            System.out.println("Partition the instances into " + numPartitions + ".");
            this.fileNames = this.getPartitionOfInstances(filteredInstanceFileNames, numPartitions, idxPartition);
            System.out.println("The partition with index " + idxPartition + " was created and contains " +
                    this.fileNames.length + " instances.");
        }else{
            throw new RuntimeException("Whether the given abstract path name for the fragmentation tree" +
                    "directory or the output directory does not exist or is not a directory.");
        }
    }

    private List<String> filterOutProcessedInstances(Collection<String> processedInstanceFileNames) {
        // For each method (computation of subtrees, comparison or structure ranking), FTrees are always needed.
        return Arrays.stream(Objects.requireNonNull(this.fTreeDir.list())).
                map(fTreeFileName -> fTreeFileName.replaceFirst("\\.json", "")).
                filter(fileName -> {
                    for (String processedFileName : processedInstanceFileNames) {
                        if (fileName.equals(processedFileName)) return false;
                    }
                    return true;
                }).collect(Collectors.toList());
    }

    private String[] getPartitionOfInstances(List<String> instanceFileNames, int numPartitions, int idxPartition) {
        int numberOfInstances = instanceFileNames.size();
        if (numberOfInstances <= numPartitions) {
            if (idxPartition <= numberOfInstances - 1) {
                return new String[]{instanceFileNames.get(idxPartition)};
            } else {
                return new String[0];
            }
        } else {
            int lengthPartition = numberOfInstances / numPartitions;
            int rest = numberOfInstances - numPartitions * lengthPartition;

            int startIndex, endIndex;// 'startIndex' is inclusive, 'endIndex' is exclusive
            if (idxPartition < rest) {
                // the partition will take one additional element and
                // all successor partitions have one additional element from the rest as well
                startIndex = idxPartition * (lengthPartition + 1);
                endIndex = startIndex + lengthPartition + 1;
            } else {
                startIndex = idxPartition * lengthPartition + rest;
                endIndex = startIndex + lengthPartition;
            }

            String[] fileNames = new String[endIndex - startIndex];
            for (int i = 0; i < fileNames.length; i++) {
                fileNames[i] = instanceFileNames.get(startIndex + i);
            }
            return fileNames;
        }
    }

    private MolecularGraph readMoleculeFromMsFile(String fileName) throws IOException, InvalidSmilesException, UnknownElementException {
        File file = new File(this.spectraDir, fileName);
        try (BufferedReader fileReader = new BufferedReader(new FileReader(file))) {

            String currentLine = fileReader.readLine();
            String molecularFormula = null, smiles = null;
            boolean mfWasRead = false;
            boolean smilesWasRead = false;

            while (currentLine != null) {
                if (currentLine.startsWith(">formula")) {
                    molecularFormula = currentLine.split(" ")[1];
                    mfWasRead = true;
                } else if (currentLine.startsWith(">smiles")) {
                    smiles = currentLine.split(" ")[1];
                    smilesWasRead = true;
                }
                if (mfWasRead && smilesWasRead) break;
                currentLine = fileReader.readLine();
            }
            // Assumption: In any case, 'file' contains two lines that start with ">formula" and ">smiles".
            // --> the molecular formula and the smiles string have been read
            SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
            return new MolecularGraph(MolecularFormula.parse(molecularFormula), smilesParser.parseSmiles(smiles));
        }
    }

    private FTree readFTree(String fileName) throws IOException {
        File file = new File(this.fTreeDir, fileName);
        FTree fTree = MsIO.readTreeFromFile(file);
        return fTree;
    }

    private CSIPredictionData readPredictionDataFromTSV(String fileName) throws IOException, UnknownElementException {
        File file = new File(this.predictionsDir, fileName);
        try(BufferedReader fileReader = Files.newBufferedReader(file.toPath())){
            // 1.: Get the index of the columns of interest:
            String[] columnNames = fileReader.readLine().split("\\t");
            int rankIdx = -1, csiScoreIdx = -1, confidenceScoreIdx = -1, mfIdx = -1, smilesIdx = -1;

            for(int i = 0; i < columnNames.length; i++){
                String columnName = columnNames[i];
                switch(columnName){
                    case "rank":
                        rankIdx = i;
                        break;
                    case "ConfidenceScore":
                        confidenceScoreIdx = i;
                        break;
                    case "CSI:FingerIDScore":
                        csiScoreIdx = i;
                        break;
                    case "molecularFormula":
                        mfIdx = i;
                        break;
                    case "smiles":
                        smilesIdx = i;
                        break;
                }
            }

            // 2.: Read all lines in the file and return a list of String arrays:
            ArrayList<String[]> lines = this.readAllRemainingLines(fileReader);

            // 3.: Now parse these lines and create a CSIPredictionData object:
            int numberOfDataLines = lines.size();
            int[] ranks = new int[numberOfDataLines];
            double[] csiScores = new double[numberOfDataLines];
            double[] confidenceScores = new double[numberOfDataLines];
            MolecularFormula[] molecularFormulas = new MolecularFormula[numberOfDataLines];
            String[] smilesStrings = new String[numberOfDataLines];

            for(int i = 0; i < numberOfDataLines; i++){
                String[] currentLine = lines.get(i);
                ranks[i] = Integer.parseInt(currentLine[rankIdx]);
                csiScores[i] = Double.parseDouble(currentLine[csiScoreIdx]);
                confidenceScores[i] = Double.parseDouble(currentLine[confidenceScoreIdx]);
                molecularFormulas[i] = MolecularFormula.parse(currentLine[mfIdx]);
                smilesStrings[i] = currentLine[smilesIdx];
            }

            return new CSIPredictionData(ranks, csiScores, confidenceScores, molecularFormulas, smilesStrings);
        }
    }

    private ArrayList<String[]> readAllRemainingLines(BufferedReader fileReader) throws IOException{
        ArrayList<String[]> lines = new ArrayList<>();
        String currentLine = fileReader.readLine();
        while(currentLine != null){
            String[] data = currentLine.split("\\t");
            lines.add(data);

            currentLine = fileReader.readLine();
        }

        return lines;
    }

    private synchronized void appendStringToFile(File file, String str) throws IOException{
        try(BufferedWriter fileWriter = Files.newBufferedWriter(file.toPath(), StandardOpenOption.APPEND)){
            fileWriter.newLine();
            fileWriter.write(str);
        }
    }

    private void unreference(MolecularGraph molecule, FTree fTree, CombinatorialFragmenterScoring scoring, CombinatorialSubtreeCalculator subtreeCalc) {
        molecule = null;
        fTree = null;
        scoring = null;
        subtreeCalc = null;
    }

    public void computeCombinatorialSubtrees(CombinatorialFragmenter.Callback2 fragmentationConstraint, SubtreeComputationMethod subtreeCompMethod) throws InterruptedException {
        // INITIALISATION:
        // Initialise the ExecutorService:
        System.out.println("Initialise ExecutorService...");
        final int NUM_CONCURRENT_THREADS = Runtime.getRuntime().availableProcessors();
        ExecutorService executor = Executors.newFixedThreadPool(NUM_CONCURRENT_THREADS);

        // LOOP - CREATE TASKS FOR EACH INSTANCE:
        // For each instance, create a task which computes the in silico fragmentation graph, computes the subtree using
        // CriticalPath1 and saves the result into a JSON file.
        System.out.println("Each instance corresponds to one task. Collect all tasks...");
        ArrayList<Callable<Object>> tasks = new ArrayList<>(this.fileNames.length);
        for (String fileName : this.fileNames) {
            Callable<Object> task = Executors.callable(() -> {
                try {
                    MolecularGraph molecule = this.readMoleculeFromMsFile(fileName + ".ms");
                    FTree fTree = this.readFTree(fileName + ".json");
                    DirectedBondTypeScoring scoring = new DirectedBondTypeScoring(molecule);

                    CombinatorialSubtreeCalculator subtreeCalc = SubtreeComputationMethod.getComputedSubtreeCalculator(fTree, molecule, scoring, fragmentationConstraint, subtreeCompMethod);

                    CombinatorialSubtreeCalculatorJsonWriter.writeResultsToFile(subtreeCalc, new File(this.outputDir, fileName + ".json"));
                    this.unreference(molecule, fTree, scoring, subtreeCalc);
                } catch (Exception e) {
                    System.out.println("An error occurred during processing instance " + fileName);
                    File resultFile = new File(this.outputDir, fileName + ".json");
                    if (resultFile.exists()) {
                        boolean wasDeleted = resultFile.delete();
                        if (wasDeleted) {
                            System.out.println(fileName + ".json was successfully deleted.");
                        } else {
                            System.out.println("Could not delete " + fileName + ".json.");
                        }
                    }
                    e.printStackTrace();
                }
            });
            tasks.add(task);
        }

        // COMPUTING THE TASKS:
        // Wait, until all tasks terminated.
        System.out.println("Execute all tasks...");
        Collection<Future<Object>> futures = executor.invokeAll(tasks);

        System.out.println(futures.size() + " tasks of " + tasks.size() + " have been processed and the executor service will be shutdown.");
        executor.shutdown();
    }

    public void compareSubtreeComputationMethods(CombinatorialFragmenter.Callback2 fragmentationConstraint, String outputFileName) throws InterruptedException, IOException, ExecutionException {
        // We know: this.fileNames contains the filenames of all instances that have to be processed.
        // For each instance, we want to compare all subtree computation methods.
        // That means, we want to compute the fragmentation graph and the subtrees, their score and running times and
        // the tanimoto coefficient between each heuristic subtree and the optimal subtree (ILP).
        // We store these results in the specified CSV output file 'this.outputDir/outputFileName'.

        // INITIALISATION:
        // Initialise the ExecutorService:
        final int NUM_AVAILABLE_PROCESSORS = Runtime.getRuntime().availableProcessors();
        ExecutorService executor = Executors.newFixedThreadPool(NUM_AVAILABLE_PROCESSORS);

        // Create an array with all SubtreeComputationMethods and determine the index of the ILP method:
        final SubtreeComputationMethod[] methods = SubtreeComputationMethod.values(); // array of the methods in the order they're declared
        final int ilpIdx = SubtreeComputationMethod.ILP.ordinal();

        // Test if 'outputFileName' is an existing file in this.outputDir
        // and if not, create a new file and write the starting string to it:
        File outputFile = new File(this.outputDir, outputFileName);
        if(outputFile.createNewFile()){
            StringBuilder startingString = new StringBuilder("instance_name,construction_runtime");
            for (int i = 0; i < methods.length; i++) startingString.append("," + methods[i].name() + "_runtime");
            for (int i = 0; i < methods.length; i++) startingString.append("," + methods[i].name() + "_score");
            for (int i = 0; i < methods.length; i++) startingString.append("," + methods[i].name() + "_tanimoto");

            try(BufferedWriter fileWriter = Files.newBufferedWriter(outputFile.toPath())){
                fileWriter.write(startingString.toString());
            }
        }

        // CREATING ALL TASKS:
        // Each instance in this.fileNames corresponds to one task.
        // Each task does the following:
        // - load the molecule and the FTree from their files in 'spectraDir' and 'fTreeDir', and create a scoring object
        // - create a CombinatorialGraph with terminal nodes and measure the running time of construction
        // - compute the CombinatorialSubtree for each method (ILP, Insertion, Prim, Critical Path1-3)
        // - measure the running the times, the score and the tanimoto coefficient against the ILP solution
        // - save all data in a string, matching the CSV ordering, and write it to the specified output file.
        System.out.println("Collect all tasks...");
        ArrayList<Callable<Object>> tasks = new ArrayList<>(this.fileNames.length);
        for (String fileName : this.fileNames) {
            Callable<Object> task = Executors.callable(() -> {
                try {
                    // 1.) Initialise the molecule, the FTree and the scoring object:
                    MolecularGraph molecule = this.readMoleculeFromMsFile(fileName + ".ms");
                    FTree fTree = this.readFTree(fileName + ".json");
                    DirectedBondTypeScoring scoring = new DirectedBondTypeScoring(molecule);

                    // 2.) Create the CombinatorialGraph and add the terminal nodes to it:
                    // Measure the running time for constructing such fragmentation graph.
                    long timeStamp, constructionRuntime;
                    CombinatorialFragmenter fragmenter = new CombinatorialFragmenter(molecule, scoring);
                    timeStamp = System.nanoTime();
                    CombinatorialGraph graph = fragmenter.createCombinatorialFragmentationGraph(fragmentationConstraint);
                    CombinatorialGraphManipulator.addTerminalNodes(graph, scoring, fTree);
                    constructionRuntime = System.nanoTime() - timeStamp;

                    // 3.) Compute the CombinatorialSubtree of 'graph' with each SubtreeComputationMethod:
                    // Save the score, the running time and the tanimoto coefficient between the ILP solution for each method.
                    CombinatorialSubtree[] subtrees = new CombinatorialSubtree[methods.length];
                    long[] runningTimes = new long[methods.length];
                    double[] scores = new double[methods.length];
                    double[] tanimotoScores = new double[methods.length];

                    // 3.1: For each method, compute the subtree, measure the runtime and score:
                    for (int i = 0; i < methods.length; i++) {
                        SubtreeComputationMethod method = methods[i];
                        timeStamp = System.nanoTime();
                        CombinatorialSubtreeCalculator subtreeCalc = SubtreeComputationMethod.getComputedSubtreeCalculator(fTree, graph, scoring, method);

                        runningTimes[i] = System.nanoTime() - timeStamp;
                        subtrees[i] = subtreeCalc.getSubtree();
                        scores[i] = subtreeCalc.getScore();
                    }

                    // 3.2: For each method, compute the tanimoto coefficient between the subtree and the ILP subtree:
                    CombinatorialSubtree ilpSubtree = subtrees[ilpIdx];
                    TObjectIntHashMap<BitSet> mergedEdgeBitSet2Index = graph.mergedEdgeBitSet2Index();
                    int maxBitSetLength = graph.maximalBitSetLength();

                    for (int i = 0; i < methods.length; i++) {
                        CombinatorialSubtree subtree = subtrees[i];
                        tanimotoScores[i] = CombinatorialSubtreeManipulator.tanimoto(subtree, ilpSubtree, mergedEdgeBitSet2Index, maxBitSetLength);
                    }

                    // 4.) Save all data into one string and write this string into the output file:
                    // The string should look like this:
                    // "<instance name>,<constrRunningtime>,{Running times},{Scores},{Tanimoto-Scores}"
                    StringBuilder strBuilder = new StringBuilder(fileName + "," + constructionRuntime);
                    for (int i = 0; i < methods.length; i++) strBuilder.append("," + runningTimes[i]);
                    for (int i = 0; i < methods.length; i++) strBuilder.append("," + scores[i]);
                    for (int i = 0; i < methods.length; i++) strBuilder.append("," + tanimotoScores[i]);

                    this.appendStringToFile(outputFile, strBuilder.toString());
                }catch(IOException | InvalidSmilesException | UnknownElementException | GRBException e){
                    System.out.println("An error occurred during computing instance"+fileName);
                    e.printStackTrace();
                }
            });
            tasks.add(task);
        }
        System.out.println("All tasks were collected and will be executed now...");
        executor.invokeAll(tasks);

        System.out.println("The executor service will be shutdown.");
        executor.shutdown();
    }

    public void runStructureRanking(CombinatorialFragmenter.Callback2 fragmentationConstraint, SubtreeComputationMethod subtreeCompMethod){
        throw new UnsupportedOperationException("This method is currently not supported.");
    }

    private class CSIPredictionData{

        protected final int[] ranks;
        protected final double[] csiScores;
        protected final double[] confidenceScores;
        protected final MolecularFormula[] molecularFormulas;
        protected final String[] smilesStrings;

        public CSIPredictionData(int[] ranks, double[] csiScores, double[] confidenceScores, MolecularFormula[] molecularFormulas, String[] smilesStrings){
            this.ranks = ranks;
            this.csiScores = csiScores;
            this.confidenceScores = confidenceScores;
            this.molecularFormulas = molecularFormulas;
            this.smilesStrings = smilesStrings;
        }

    }

}
