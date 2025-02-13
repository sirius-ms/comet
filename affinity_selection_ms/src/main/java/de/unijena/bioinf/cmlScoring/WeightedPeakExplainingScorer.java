package de.unijena.bioinf.cmlScoring;


import de.unijena.bioinf.ChemistryBase.algorithm.scoring.SScored;
import de.unijena.bioinf.ChemistryBase.algorithm.scoring.Score;
import de.unijena.bioinf.ChemistryBase.chem.PrecursorIonType;
import de.unijena.bioinf.ChemistryBase.ms.*;
import de.unijena.bioinf.ChemistryBase.ms.utils.OrderedSpectrum;
import de.unijena.bioinf.ChemistryBase.ms.utils.SimpleSpectrum;
import de.unijena.bioinf.ChemistryBase.ms.utils.Spectrums;
import de.unijena.bioinf.cmlSpectrumPrediction.BBFragmentSpectrumPredictor;
import de.unijena.bioinf.datastructures.CMLMolecule;
import de.unijena.bioinf.sirius.ProcessedInput;
import de.unijena.bioinf.sirius.ProcessedPeak;
import de.unijena.bioinf.sirius.Sirius;
import de.unijena.bionf.spectral_alignment.ModifiedCosine;
import de.unijena.bionf.spectral_alignment.SpectralSimilarity;

import java.util.HashMap;
import java.util.List;

public class WeightedPeakExplainingScorer implements CMLScorer<Score.DoubleScore> {

    private final Deviation massDeviation;
    private final int numHydrogenShifts;
    private final double hydrogenShiftProbability;
    private final HashMap<String, Double> fragmentType2Frequency;

    public WeightedPeakExplainingScorer(double ppm, int numHydrogenShifts, double hydrogenShiftProbability, HashMap<String, Double> fragmentType2Frequency) {
        this.massDeviation = new Deviation(ppm);
        this.numHydrogenShifts = numHydrogenShifts;
        this.hydrogenShiftProbability = hydrogenShiftProbability;
        this.fragmentType2Frequency = fragmentType2Frequency;
    }

    @Override
    public SScored<CMLMolecule, Score.DoubleScore> score(Ms2Experiment exp, CMLMolecule molecule) {
        final PrecursorIonType precursorIonType = exp.getPrecursorIonType().isIonizationUnknown() ? PrecursorIonType.fromString("[M+H]+") : exp.getPrecursorIonType();
        final OrderedSpectrum<Peak> msrdSpectrum = this.getMeasuredSpectrum(exp); // is normalized to sum of peak intensities equals 1

        final BBFragmentSpectrumPredictor spectrumPredictor = new BBFragmentSpectrumPredictor(molecule, this.fragmentType2Frequency, precursorIonType, this.numHydrogenShifts, this.hydrogenShiftProbability);
        final OrderedSpectrum<Peak> predSpectrum = new SimpleSpectrum(spectrumPredictor.predictSpectrum());

        final ModifiedCosine modifiedCosine = new ModifiedCosine(this.massDeviation);
        final SpectralSimilarity spectralSimilarity = modifiedCosine.score(msrdSpectrum, predSpectrum, exp.getIonMass(), spectrumPredictor.getSpectrum().getPrecursorMz());
        return new SScored<>(molecule, new Score.DoubleScore(spectralSimilarity.similarity));
    }

    private OrderedSpectrum<Peak> getMeasuredSpectrum(Ms2Experiment exp){
        ProcessedInput processedInput = new Sirius().preprocessForMs2Analysis(exp);
        List<ProcessedPeak> mergedPeaks = processedInput.getMergedPeaks();
        mergedPeaks.removeLast();

        MutableMs2Spectrum msrdSpectrum = new MutableMs2Spectrum();
        msrdSpectrum.setPrecursorMz(exp.getIonMass());
        mergedPeaks.forEach(msrdSpectrum::addPeak);

        Spectrums.normalizeToSum(msrdSpectrum, 1d);

        return new SimpleSpectrum(msrdSpectrum);
    }
}
