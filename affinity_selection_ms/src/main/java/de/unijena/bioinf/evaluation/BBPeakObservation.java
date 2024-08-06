package de.unijena.bioinf.evaluation;

import de.unijena.bioinf.ChemistryBase.chem.MolecularFormula;
import de.unijena.bioinf.ChemistryBase.chem.PrecursorIonType;
import de.unijena.bioinf.ChemistryBase.chem.Smiles;
import de.unijena.bioinf.ChemistryBase.chem.utils.UnknownElementException;
import de.unijena.bioinf.ChemistryBase.ms.Deviation;
import de.unijena.bioinf.ChemistryBase.ms.Ms2Experiment;
import de.unijena.bioinf.ChemistryBase.ms.Peak;
import de.unijena.bioinf.ChemistryBase.ms.Spectrum;
import de.unijena.bioinf.ChemistryBase.ms.utils.SimpleMutableSpectrum;
import de.unijena.bioinf.ChemistryBase.ms.utils.Spectrums;
import de.unijena.bioinf.babelms.MsIO;
import de.unijena.bioinf.cmlSpectrumPrediction.BBBarcodeSpectrumPredictor;
import de.unijena.bioinf.fragmenter.MolecularGraph;
import de.unijena.bioinf.sirius.ProcessedInput;
import de.unijena.bioinf.sirius.ProcessedPeak;
import de.unijena.bioinf.sirius.Sirius;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.*;

public class BBPeakObservation {

    private static BitSet incrementBitSet(BitSet bitSet){
        BitSet newBitSet = (BitSet) bitSet.clone();
        int i = 0;
        while(newBitSet.get(i)){
            newBitSet.set(i, false);
            i++;
        }
        newBitSet.set(i, true);
        return newBitSet;
    }

    public static void main(String[] args){
        try {
            final File spectraDir = new File(args[0]);
            final File outputFile = new File(args[1]);
            final File bbFile = new File(args[2]);
            final double scaffoldMass = Double.parseDouble(args[3]);
            final Deviation deviation = new Deviation(Double.parseDouble(args[4]));
            final int HYDROGEN_SHIFTS = Integer.parseInt(args[5]);

            final HashMap<Integer, ArrayList<Double>> bbPos2Mass = parseBBCSV(bbFile);
            final int numBBs = bbPos2Mass.keySet().size();

            final BitSet[] bbCombinations = new BitSet[(int) Math.pow(2,numBBs) + numBBs - 1];
            int k = 0;
            BitSet bbCombination = new BitSet(numBBs);
            while (bbCombination.cardinality() < numBBs) {
                bbCombinations[k++] = bbCombination;
                bbCombination = incrementBitSet(bbCombination);
            }
            bbCombination = new BitSet(numBBs + 1);
            bbCombination.set(numBBs);
            for (int i = 0; i < numBBs; i++) {
                BitSet clonedBBCombination = (BitSet) bbCombination.clone();
                clonedBBCombination.set(i);
                bbCombinations[k++] = clonedBBCombination;
            }

            try (BufferedWriter fileWriter = Files.newBufferedWriter(outputFile.toPath(), StandardCharsets.UTF_8)) {
                // Header:
                StringBuilder strBuilder = new StringBuilder("fileName");
                for (BitSet bitSet : bbCombinations) {
                    strBuilder.append(",").append(getBBCombinationString(bitSet, numBBs));
                }
                fileWriter.write(strBuilder.toString());
                fileWriter.newLine();

                // Data rows:
                SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
                final String[] fileNames = Objects.requireNonNull(spectraDir.list());
                for (final String fileName : fileNames) {
                    //Parse all information:
                    final Ms2Experiment ms2Experiment = MsIO.readExperimentFromFile(new File(spectraDir, fileName)).next();
                    final ProcessedInput processedInput = new Sirius().preprocessForMs2Analysis(ms2Experiment);
                    final int[] bbIndices = parseName(ms2Experiment.getName(), numBBs);
                    final Spectrum<Peak> msrdSpectrum = getMsrdMs2Spectrum(processedInput);
                    final MolecularGraph molecule = getMolecularGraph(processedInput, smilesParser);
                    final PrecursorIonType ionization = ms2Experiment.getPrecursorIonType();

                    // Predict the spectrum / initialise an object of BBBarcodeSpectrumPredictor
                    final double[] bbMasses = new double[numBBs];
                    for(int i = 0; i < numBBs; i++)
                        bbMasses[i] = bbPos2Mass.get(i).get(bbIndices[i]);
                    final BBBarcodeSpectrumPredictor predictor = new BBBarcodeSpectrumPredictor(molecule, bbMasses, scaffoldMass, HYDROGEN_SHIFTS, ionization);
                    predictor.predictSpectrum();

                    // For each building block combination,
                    // test if at least one corresponding peak is contained in msrdSpectrum
                    strBuilder = new StringBuilder(fileName);
                    for(final BitSet bitSet : bbCombinations){
                        final List<Peak> peaks = predictor.getBitset2Peaks().get(bitSet);
                        boolean containsBBPeak = false;
                        for(Peak peak : peaks){
                            if(containsBBPeak(msrdSpectrum, peak, deviation)){
                                containsBBPeak = true;
                                break;
                            }
                        }
                        strBuilder.append(",").append(containsBBPeak ? 1 : 0);
                    }
                    fileWriter.write(strBuilder.toString());
                    fileWriter.newLine();
                }
            } catch (IOException e) {
                throw new RuntimeException(e);
            }
        } catch (UnknownElementException | InvalidSmilesException e) {
            throw new RuntimeException(e);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    private static String getBBCombinationString(BitSet bbCombination, int numBBs){
        if(bbCombination.length() <= numBBs){
            int bbIdx = bbCombination.nextSetBit(0);
            if(bbIdx < 0) return "S";

            final StringBuilder strBuilder = new StringBuilder("S[").append(bbIdx);
            bbIdx = bbCombination.nextSetBit(bbIdx+1);
            while(bbIdx >= 0){
                strBuilder.append(';').append(bbIdx);
                bbIdx = bbCombination.nextSetBit(bbIdx+1);
            }
            strBuilder.append("]");
            return strBuilder.toString();
        }else{
            return Integer.toString(bbCombination.nextSetBit(0));
        }
    }

    private static boolean containsBBPeak(Spectrum<Peak> msrdSpectrum, Peak peak, Deviation deviation){
        return Spectrums.binarySearch(msrdSpectrum, peak.getMass(), deviation) >= 0;
    }

    private static MolecularGraph getMolecularGraph(ProcessedInput processedInput, SmilesParser smiParser) throws InvalidSmilesException {
        final String smiles = processedInput.getAnnotation(Smiles.class).orElseThrow().toString();
        return new MolecularGraph(smiParser.parseSmiles(smiles));
    }

    public static int[] parseName(String name, int numBBs){
        final String[] strIndices = name.split("-");
        if(strIndices.length != numBBs) throw new RuntimeException("Name is invalid!");
        final int[] indices = new int[strIndices.length];
        for(int i = 0; i < strIndices.length; i++){
            indices[i] = Integer.parseInt(strIndices[i]) - 1;
        }
        return indices;
    }

    public static Spectrum<Peak> getMsrdMs2Spectrum(ProcessedInput processedInput) throws IOException{
        List<ProcessedPeak> mergedPeaks = processedInput.getMergedPeaks();
        SimpleMutableSpectrum spec = new SimpleMutableSpectrum(mergedPeaks.size());
        for(Peak mergedPeak : mergedPeaks) spec.addPeak(mergedPeak);
        return spec;
    }

    public static HashMap<Integer, ArrayList<Double>> parseBBCSV(File file) throws IOException, UnknownElementException {
        try(BufferedReader fileReader = Files.newBufferedReader(file.toPath())){
            String currentLine = fileReader.readLine();

            final String[] columnNames = currentLine.split(",");
            final HashMap<String, Integer> name2Idx = new HashMap<>();
            for(int i = 0; i < columnNames.length; i++) name2Idx.put(columnNames[i], i);

            final HashMap<Integer, ArrayList<Double>> bbPos2Masses = new HashMap<>();
            currentLine = fileReader.readLine();
            while(currentLine != null){
                final String[] dataArray = currentLine.split(",");
                final int bbPos = Integer.parseInt(dataArray[name2Idx.get("bb_pos")]);
                final double bbMass = getBBMass(dataArray[name2Idx.get("formula")], dataArray[name2Idx.get("reaction_loss")]);
                bbPos2Masses.computeIfAbsent(bbPos, k -> new ArrayList<>()).add(bbMass);
                currentLine = fileReader.readLine();
            }

            return bbPos2Masses;
        }
    }

    public static double getBBMass(String mf, String rl) throws UnknownElementException {
        MolecularFormula molFormula = MolecularFormula.parse(mf);
        MolecularFormula reactionLoss = MolecularFormula.parse(rl);
        return molFormula.getMass() - reactionLoss.getMass();
    }
}
