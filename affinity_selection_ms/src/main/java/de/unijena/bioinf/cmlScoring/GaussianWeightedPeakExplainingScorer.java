package de.unijena.bioinf.cmlScoring;

import de.unijena.bioinf.ChemistryBase.algorithm.scoring.SScored;
import de.unijena.bioinf.ChemistryBase.algorithm.scoring.Score;
import de.unijena.bioinf.ChemistryBase.chem.MolecularFormula;
import de.unijena.bioinf.ChemistryBase.chem.PrecursorIonType;
import de.unijena.bioinf.ChemistryBase.chem.utils.UnknownElementException;
import de.unijena.bioinf.ChemistryBase.ms.*;
import de.unijena.bioinf.ChemistryBase.ms.utils.OrderedSpectrum;
import de.unijena.bioinf.ChemistryBase.ms.utils.SimpleSpectrum;
import de.unijena.bioinf.ChemistryBase.ms.utils.Spectrums;
import de.unijena.bioinf.babelms.MsIO;
import de.unijena.bioinf.cmlSpectrumPrediction.BBFragmentSpectrumPredictor;
import de.unijena.bioinf.datastructures.BuildingBlock;
import de.unijena.bioinf.datastructures.CMLMolecule;
import de.unijena.bioinf.datastructures.CombinatorialMoleculeLibrary;
import de.unijena.bioinf.datastructures.Scaffold;
import de.unijena.bioinf.io.BuildingBlockReader;
import de.unijena.bioinf.sirius.ProcessedInput;
import de.unijena.bioinf.sirius.ProcessedPeak;
import de.unijena.bioinf.sirius.Sirius;
import de.unijena.bionf.spectral_alignment.GaussianSpectralMatching;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;

public class GaussianWeightedPeakExplainingScorer implements CMLScorer<Score.DoubleScore> {

    private final Deviation massDeviation;
    private final int allowedHydrogenShifts;
    private final double hydrogenShiftProbability;
    private final HashMap<String, Double> fragmentType2Frequency;

    public GaussianWeightedPeakExplainingScorer(double ppm, int allowedHydrogenShifts, HashMap<String, Double> fragmentType2Frequency, double hydrogenShiftProbability) {
        this.massDeviation = new Deviation(ppm);
        this.allowedHydrogenShifts = allowedHydrogenShifts;
        this.hydrogenShiftProbability = hydrogenShiftProbability;
        this.fragmentType2Frequency = fragmentType2Frequency;
    }


    @Override
    public SScored<CMLMolecule, Score.DoubleScore> score(Ms2Experiment exp, CMLMolecule molecule) {
        final PrecursorIonType precursorIonType = exp.getPrecursorIonType().isIonizationUnknown() ? PrecursorIonType.fromString("[M+H]+") : exp.getPrecursorIonType();
        final OrderedSpectrum<Peak> msrdSpectrum = this.getMeasuredSpectrum(exp);

        final BBFragmentSpectrumPredictor spectrumPredictor = new BBFragmentSpectrumPredictor(molecule, this.fragmentType2Frequency, precursorIonType, this.allowedHydrogenShifts, this.hydrogenShiftProbability);
        final OrderedSpectrum<Peak> predSpectrum = new SimpleSpectrum(spectrumPredictor.predictSpectrum());

        final GaussianSpectralMatching spectralMatching = new GaussianSpectralMatching(this.massDeviation);
        double score = spectralMatching.score(msrdSpectrum, predSpectrum).similarity;
        score = this.normalizeGaussianDotProduct(score, msrdSpectrum, predSpectrum, spectralMatching);
        return new SScored<>(molecule, new Score.DoubleScore(score));
    }

    private double normalizeGaussianDotProduct(double score, OrderedSpectrum<Peak> msrdSpectrum, OrderedSpectrum<Peak> predSpectrum, GaussianSpectralMatching spectralMatching) {
        double normMsrdSpectrum = Math.sqrt(spectralMatching.score(msrdSpectrum,msrdSpectrum).similarity);
        double normPredSpectrum = Math.sqrt(spectralMatching.score(predSpectrum,predSpectrum).similarity);
        return score / (normMsrdSpectrum * normPredSpectrum);
    }

    private OrderedSpectrum<Peak> getMeasuredSpectrum(Ms2Experiment exp) {
        ProcessedInput processedInput = new Sirius().preprocessForMs2Analysis(exp);
        List<ProcessedPeak> mergedPeaks = processedInput.getMergedPeaks();
        mergedPeaks.removeLast(); // the peaks in mergedPeaks are sorted according to their mass in ascending order and the last peak is the precursor peak

        MutableMs2Spectrum msrdSpectrum = new MutableMs2Spectrum();
        msrdSpectrum.setPrecursorMz(exp.getIonMass());
        mergedPeaks.forEach(msrdSpectrum::addPeak);

        Spectrums.normalizeToSum(msrdSpectrum, 1d);

        return new SimpleSpectrum(msrdSpectrum);
    }

    public static void main(String[] arx) {
        try {
            // LIBRARY INIT:
            BuildingBlock[][] bbs = BuildingBlockReader.readBuildingBlocks(new File("/home/nils/Dokumente/Bioinformatik_PhD/AS-MS-Project/Data/benzimidazole/ENL173/bbs_and_compounds/ENL173_CustomDB_BBs.csv"));
            Scaffold scaffold = new Scaffold(MolecularFormula.parse("C8H3N2O"), null);
            CombinatorialMoleculeLibrary cmlLibrary = new CombinatorialMoleculeLibrary(bbs, scaffold);

            List<CMLMolecule> mols = cmlLibrary.generateMolecules();
            HashMap<String, CMLMolecule> name2Molecule = new HashMap<>(mols.size());
            mols.forEach(mol -> name2Molecule.put(mol.getName(), mol));

            // GENERAL INIT:
            HashMap<String, Double> fragmentType2Frequency = new HashMap<>();
            fragmentType2Frequency.put("S[0]", 0d);
            fragmentType2Frequency.put("S[1]", 0.01293103448275862);
            fragmentType2Frequency.put("S[0;1]", 0.017241379310344827);
            fragmentType2Frequency.put("S[2]", 0.7327586206896551);
            fragmentType2Frequency.put("S[0;2]", 0.5);
            fragmentType2Frequency.put("S[1;2]", 0.875);
            fragmentType2Frequency.put("0", 0.12931034482758622);
            fragmentType2Frequency.put("1", 0.34051724137931033);
            fragmentType2Frequency.put("2", 0d);

            // MEASURED SPECTRUM:
            File msFile = new File("/home/nils/Dokumente/Bioinformatik_PhD/AS-MS-Project/Ergebnisse/CML_Scoring_Evaluation/BBFragment_SimplePeakExplaining_Scoring/ENL161/spectra/ENL161_spectrum_100.ms");
            Ms2Experiment exp = MsIO.readExperimentFromFile(msFile).next();

            // MOLECULE SCORING:
            GaussianWeightedPeakExplainingScorer scorer = new GaussianWeightedPeakExplainingScorer(5, 2, fragmentType2Frequency, 0.5);

            CMLMolecule molecule = name2Molecule.get("23-15-37");
            System.out.println(scorer.score(exp, molecule).getScore());

        } catch (UnknownElementException | IOException e) {
            throw new RuntimeException(e);
        }
    }
}
