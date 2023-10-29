package de.unijena.bioinf.cmlSpectrumPrediction;

import com.fasterxml.jackson.core.JsonFactory;
import com.fasterxml.jackson.core.JsonGenerator;
import com.fasterxml.jackson.core.JsonParser;
import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.node.ArrayNode;
import com.fasterxml.jackson.databind.node.ObjectNode;
import de.unijena.bioinf.ChemistryBase.chem.MolecularFormula;
import de.unijena.bioinf.ChemistryBase.chem.PrecursorIonType;
import de.unijena.bioinf.ChemistryBase.data.JSONDocumentType;
import de.unijena.bioinf.ChemistryBase.jobs.SiriusJobs;
import de.unijena.bioinf.ChemistryBase.ms.*;
import de.unijena.bioinf.ChemistryBase.ms.ft.FTree;
import de.unijena.bioinf.ChemistryBase.ms.ft.Fragment;
import de.unijena.bioinf.ChemistryBase.ms.ft.model.Whiteset;
import de.unijena.bioinf.ChemistryBase.ms.utils.OrderedSpectrum;
import de.unijena.bioinf.ChemistryBase.ms.utils.SimpleMutableSpectrum;
import de.unijena.bioinf.ChemistryBase.ms.utils.SimpleSpectrum;
import de.unijena.bioinf.FragmentationTreeConstruction.computation.FasterTreeComputationInstance;
import de.unijena.bioinf.FragmentationTreeConstruction.computation.FragmentationPatternAnalysis;
import de.unijena.bioinf.babelms.MsIO;
import de.unijena.bioinf.cmlFragmentation.*;
import de.unijena.bioinf.fragmenter.CombinatorialFragment;
import de.unijena.bioinf.fragmenter.CriticalPathSubtreeCalculator;
import de.unijena.bioinf.fragmenter.EMFragmenterScoring2;
import de.unijena.bioinf.fragmenter.MolecularGraph;
import de.unijena.bioinf.jjobs.JobManager;
import de.unijena.bioinf.sirius.ProcessedInput;
import de.unijena.bioinf.sirius.ProcessedPeak;
import de.unijena.bioinf.sirius.Sirius;
import de.unijena.bionf.spectral_alignment.RecallSpectralAlignment;
import de.unijena.bionf.spectral_alignment.WeightedRecallSpectralAlignment;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

public class Main {

    public static class SimpleFragmentationRule implements FragmentationRules {

        private List<String> allowedElements;

        public SimpleFragmentationRule(String[] allowedElements){
            this.allowedElements = Arrays.asList(allowedElements);
        }

        @Override
        public boolean match(IBond bond) {
            if(!bond.getOrder().equals(IBond.Order.SINGLE)) return false;
            String atom1Symbol = bond.getAtom(0).getSymbol();
            String atom2Symbol = bond.getAtom(1).getSymbol();
            return (this.allowedElements.contains(atom1Symbol) || this.allowedElements.contains(atom2Symbol));
        }
    }

    public enum FragmentationPredictorType {
        ITERATIVE, RULE_BASED;
    }

    private static SimpleSpectrum parseMsrdSpectrum(ProcessedInput processedMs2Experiment) throws IOException {
        List<ProcessedPeak> mergedPeaks = processedMs2Experiment.getMergedPeaks();
        SimpleMutableSpectrum s = new SimpleMutableSpectrum(mergedPeaks.size());
        for(ProcessedPeak peak : mergedPeaks) s.addPeak(peak);
        s.removePeakAt(s.size()-1);
        return new SimpleSpectrum(s);
    }

    private static HashMap<Peak, CombinatorialFragment> getEpimetheusMapping(MolecularGraph molecule, FTree fTree, Spectrum<Peak> spectrum){
        final EMFragmenterScoring2 scoring = new EMFragmenterScoring2(molecule, fTree);

        final HashSet<MolecularFormula> mfSet = new HashSet<>();
        for (Fragment ft : fTree.getFragmentsWithoutRoot()) {
            mfSet.add(ft.getFormula());
            mfSet.add(ft.getFormula().add(MolecularFormula.getHydrogen()));
            mfSet.add(ft.getFormula().add(MolecularFormula.getHydrogen().multiply(2)));
            if (ft.getFormula().numberOfHydrogens()>0) mfSet.add(ft.getFormula().subtract(MolecularFormula.getHydrogen()));
            if (ft.getFormula().numberOfHydrogens()>1) mfSet.add(ft.getFormula().subtract(MolecularFormula.getHydrogen().multiply(2)));
        }

        final CriticalPathSubtreeCalculator subtreeCalc = new CriticalPathSubtreeCalculator(fTree, molecule, scoring, true);
        subtreeCalc.setMaxNumberOfNodes(100000);
        subtreeCalc.initialize((node, nnodes, nedges) -> {
            if(mfSet.contains(node.getFragment().getFormula())) return true;
            return (node.getTotalScore() > -5f);
        });
        subtreeCalc.computeSubtree();

        final HashMap<Fragment, ArrayList<CombinatorialFragment>> mapping = subtreeCalc.computeMapping();
        final HashMap<Peak, CombinatorialFragment> peak2Fragment = new HashMap<>();
        for(Fragment ft : mapping.keySet()){
            if(ft.isRoot()) continue;
            Peak peak = spectrum.getPeakAt(ft.getPeakId());
            peak2Fragment.put(peak, mapping.get(ft).get(0));
        }
        return peak2Fragment;
    }

    private static void writeSpectrumComparison2Json(JsonGenerator jsonGenerator, String name, Deviation deviation, OrderedSpectrum<Peak> msrdSpectrum, OrderedSpectrum<Peak> predSpectrum, Map<Peak, CombinatorialFragment> msrdPeak2Fragment, Map<Peak, CombinatorialFragment> predPeak2Fragment) throws IOException {
        RecallSpectralAlignment recall = new RecallSpectralAlignment(deviation);
        WeightedRecallSpectralAlignment weightedRecall = new WeightedRecallSpectralAlignment(deviation);
        List<Peak> matchedMsrdPeaks = recall.getMatchedMsrdPeaks(msrdSpectrum, predSpectrum);

        jsonGenerator.writeStartObject();
        jsonGenerator.writeStringField("name", name);
        jsonGenerator.writeNumberField("recall", recall.score(msrdSpectrum, predSpectrum).similarity);
        jsonGenerator.writeNumberField("weighted_recall", weightedRecall.score(msrdSpectrum, predSpectrum).similarity);
        writeSpectrum2Json(jsonGenerator, "msrd_spectrum", msrdSpectrum, msrdPeak2Fragment, matchedMsrdPeaks);
        writeSpectrum2Json(jsonGenerator, "pred_spectrum", predSpectrum, predPeak2Fragment, Collections.emptyList());
        jsonGenerator.writeEndObject();
    }

    private static void writeSpectrum2Json(JsonGenerator jsonGenerator, String fieldName, OrderedSpectrum<Peak> spectrum, Map<Peak, CombinatorialFragment> peak2Fragment, List<Peak> matchedMsrdPeaks) throws IOException {
        jsonGenerator.writeFieldName(fieldName);

        jsonGenerator.writeStartArray();
        for(Peak peak : spectrum){
            BitSet bitSet = Optional.ofNullable(peak2Fragment.get(peak)).map(CombinatorialFragment::getBitSet).orElse(new BitSet());
            int[] atomIndices = new int[bitSet.cardinality()];
            int k = 0;
            for(int idx = bitSet.nextSetBit(0); idx >= 0; idx = bitSet.nextSetBit(idx+1))
                atomIndices[k++] = idx;

            writePeak2Json(jsonGenerator, peak, atomIndices, matchedMsrdPeaks.contains(peak));
        }
        jsonGenerator.writeEndArray();
    }

    private static void writePeak2Json(JsonGenerator jsonGenerator, Peak peak, int[] atomIndices, boolean isMatchedMsrdPeak) throws IOException {
        jsonGenerator.writeStartObject();
        jsonGenerator.writeNumberField("mz", peak.getMass());
        jsonGenerator.writeNumberField("intensity", peak.getIntensity());

        jsonGenerator.writeFieldName("atom_indices");
        jsonGenerator.writeArray(atomIndices, 0, atomIndices.length);

        jsonGenerator.writeBooleanField("isMatched", isMatchedMsrdPeak);
        jsonGenerator.writeEndObject();
    }




    public static void main(String[] args) {
        try {
            // GENERAL INITIALISATION:
            final String smiles = "CCC(CC)N1C2=C(C=C(C=C2)C(=O)NC(C(C)C)C(=O)N)N=C1CC(C)C";
            final File msFile = new File("C:\\Users\\Nutzer\\Documents\\Bioinformatik_PhD\\AS-MS-Project\\LCMS_Benzimidazole\\BAMS-14-3\\ProjectSpaces\\filtered_by_hand_PS_5.7.2\\1915_230220_BAMS-14-3_01_1915\\spectrum.ms");
            final int NUM_FRAGMENTS = 50;
            final int NUM_H_SHIFTS = 2;
            final PrecursorIonType ionization = PrecursorIonType.fromString("[M+H]+");
            final Deviation deviation = new Deviation(5);

            final SmilesParser smiParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
            final MolecularGraph molecule = new MolecularGraph(smiParser.parseSmiles(smiles));
            final MolecularFormula mf = molecule.getFormula();

            // PROCESS MS-INPUT AND COMPUTE PEAK-FRAGMENT MAPPING WITH EPIMETHEUS:
            final MutableMs2Experiment ms2Experiment = new MutableMs2Experiment(MsIO.readExperimentFromFile(msFile).next());
            ms2Experiment.setMolecularFormula(mf);
            Set<MolecularFormula> wh = Collections.singleton(mf);
            ms2Experiment.setAnnotation(Whiteset.class, Whiteset.ofMeasuredOrNeutral(wh));

            final ProcessedInput processedMs2Experiment = new Sirius().preprocessForMs2Analysis(ms2Experiment);
            final FTree fTree = new Sirius().compute(processedMs2Experiment.getExperimentInformation(), mf).getTree();

            final SimpleSpectrum msrdSpectrum = parseMsrdSpectrum(processedMs2Experiment);
            final HashMap<Peak, CombinatorialFragment> epimetheusPeak2Fragment = getEpimetheusMapping(molecule, fTree, msrdSpectrum);

            // PREDICTION WITH ICEBERG:
            final ICEBERGSpectrumPredictor icebergPredictor = new ICEBERGSpectrumPredictor(smiles, ionization, NUM_FRAGMENTS);
            final SimpleSpectrum icebergSpectrum = new SimpleSpectrum(icebergPredictor.predictSpectrum());

            // PREDICT THE FRAGMENTATION PROCESS & A SPECTRUM WITH MY METHOD:
            final FragStepDependentScoring scoring = new FragStepDependentScoring(molecule);
            final FragmentationPredictorType type = FragmentationPredictorType.ITERATIVE;
            AbstractFragmentationPredictor fragmentationPredictor;
            switch (type) {
                case RULE_BASED -> {
                    String[] allowedElements = new String[]{"N", "O", "P", "S"};
                    SimpleFragmentationRule fragRule = new SimpleFragmentationRule(allowedElements);
                    fragmentationPredictor = new RuleBasedFragmentation(molecule, scoring, NUM_FRAGMENTS, fragRule, (node, nnodes, nedges) -> true);
                }
                default -> fragmentationPredictor = new PrioritizedIterativeFragmentation(molecule, scoring, NUM_FRAGMENTS);
            }
            fragmentationPredictor.predictFragmentation();

            // Build a barcode spectrum:
            BarcodeSpectrumPredictor spectrumPredictor = new BarcodeSpectrumPredictor(fragmentationPredictor, ionization.isPositive(), NUM_H_SHIFTS);
            SimpleSpectrum predictedSpectrum = new SimpleSpectrum(spectrumPredictor.predictSpectrum());


            // WRITE RESULTS INTO A JSON-FILE:
            String outputFilePath = "C:\\Users\\Nutzer\\Documents\\Repositories\\sirius-libs\\affinity_selection_ms\\src\\main\\resources\\visualization.json";
            try(FileWriter writer = new FileWriter(outputFilePath)) {
                JsonFactory factory = new JsonFactory();
                factory.enable(JsonParser.Feature.ALLOW_NON_NUMERIC_NUMBERS);
                JsonGenerator jsonGenerator = factory.createGenerator(writer);
                jsonGenerator.useDefaultPrettyPrinter();

                jsonGenerator.writeStartObject();
                jsonGenerator.writeStringField("smiles", smiles);
                jsonGenerator.writeNumberField("precursor_mz", ionization.neutralMassToPrecursorMass(molecule.getFormula().getMass()));
                jsonGenerator.writeStringField("ionization", ionization.toString());

                jsonGenerator.writeFieldName("spectra");
                jsonGenerator.writeStartArray();
                writeSpectrumComparison2Json(jsonGenerator, "my_prediction", deviation, msrdSpectrum, predictedSpectrum, epimetheusPeak2Fragment, spectrumPredictor.getPeak2FragmentMapping());
                writeSpectrumComparison2Json(jsonGenerator, "iceberg_prediction", deviation, msrdSpectrum, icebergSpectrum, epimetheusPeak2Fragment, icebergPredictor.getPeak2FragmentMapping());
                jsonGenerator.writeEndArray();
                jsonGenerator.writeEndObject();

                jsonGenerator.close();
            }

        } catch (InvalidSmilesException e) {
            throw new RuntimeException(e);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }
}
