package de.unijena.bioinf.cmlSpectrumPrediction;

import de.unijena.bioinf.ChemistryBase.chem.MolecularFormula;
import de.unijena.bioinf.ChemistryBase.chem.PrecursorIonType;
import de.unijena.bioinf.ChemistryBase.ms.*;
import de.unijena.bioinf.ChemistryBase.ms.utils.SimpleMutableSpectrum;
import de.unijena.bioinf.ChemistryBase.ms.utils.Spectrums;
import de.unijena.bioinf.datastructures.*;
import lombok.Getter;

import java.util.HashMap;
import java.util.List;

@Getter
public class BBFragmentSpectrumPredictor implements SpectrumPredictor<Peak, Ms2Spectrum<Peak>>{

    private static final double HYDROGEN_MASS = MolecularFormula.getHydrogen().getMass();

    private final CMLMolecule molecule;
    private final PrecursorIonType precursorIonType;
    private final HashMap<String, Double> fragmentType2Frequency;
    private final int numHydrogenShifts;
    private final double hydrogenShiftProbability;
    private Ms2Spectrum<Peak> spectrum;

    public BBFragmentSpectrumPredictor(CMLMolecule molecule, HashMap<String, Double> fragmentType2Frequency, PrecursorIonType precursorIonType, int numHydrogenShifts, double hydrogenShiftProbability) {
        this.molecule = molecule;
        this.precursorIonType = precursorIonType;
        this.fragmentType2Frequency = fragmentType2Frequency;
        this.numHydrogenShifts = numHydrogenShifts;
        this.hydrogenShiftProbability = hydrogenShiftProbability;
    }

    @Override
    public Ms2Spectrum<Peak> predictSpectrum() {
        final List<BBFragment> bbFragments = this.molecule.createAllBBFragments();
        final double precursorMz = this.precursorIonType.neutralMassToPrecursorMass(this.molecule.getMass());
        final SimpleMutableSpectrum spec = new SimpleMutableSpectrum(bbFragments.size() * (2*this.numHydrogenShifts + 1));

        for(final BBFragment fragment : bbFragments) {
            final double fragmentMz = this.precursorIonType.neutralMassToPrecursorMass(fragment.getMass());
            final double fragmentIntensity = this.fragmentType2Frequency.get(fragment.getFragmentTypeString());

            for(int h = -this.numHydrogenShifts; h <= this.numHydrogenShifts; h++) {
                double shiftedFragmentMz = fragmentMz + h * HYDROGEN_MASS;
                double shiftedFragmentIntensity = fragmentIntensity * Math.pow(this.hydrogenShiftProbability, Math.abs(h));
                SimplePeak peak = new SimplePeak(shiftedFragmentMz, shiftedFragmentIntensity);
                spec.addPeak(peak);
            }
        }

        Spectrums.sortSpectrumByMass(spec);
        Spectrums.normalizeToSum(spec, 1d);

        this.spectrum = new MutableMs2Spectrum(spec, precursorMz, null, 2);
        return this.spectrum;
    }

    /*
    public static void main(String[] args){
        try {
            // LIBRARY INIT:
            BuildingBlock[][] bbs = BuildingBlockReader.readBuildingBlocks(new File("/home/nils/Dokumente/Bioinformatik_PhD/AS-MS-Project/Data/benzimidazole/ENL161/bbs_and_compounds/ENL161_CustomDB_BBs.csv"));
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

            PrecursorIonType ionType = PrecursorIonType.fromString("[M+H]+");
            int numHydrogenShifts = 2;

            // MOLECULE:
            CMLMolecule molecule = name2Molecule.get("2-3-7");

            // SPECTRUM PREDICTOR:
            BBFragmentSpectrumPredictor spectrumPredictor = new BBFragmentSpectrumPredictor(molecule, fragmentType2Frequency, ionType, numHydrogenShifts);
            Ms2Spectrum<Peak> predictedSpectrum = spectrumPredictor.predictSpectrum();

            predictedSpectrum.forEach(p -> System.out.print("(" + p.getMass() + "," + p.getIntensity() + ")\t"));
            System.out.println();

            // MEASURED SPECTRUM:
            File msFile = new File("/home/nils/Dokumente/Bioinformatik_PhD/AS-MS-Project/Data/benzimidazole/ENL161/annotation/annotated_spectra/post_filtered_spectra/manual/ENL161_spectrum_170.ms");
            Ms2Experiment exp = MsIO.readExperimentFromFile(msFile).next();
            ProcessedInput processedInput = new Sirius().preprocessForMs2Analysis(exp);
            List<ProcessedPeak> mergedPeaks = processedInput.getMergedPeaks();
            mergedPeaks.removeLast();

            // normalize the fragment peaks of the measured spectrum:
            double sumIntensity = 0d;
            for(final ProcessedPeak p : mergedPeaks) sumIntensity += p.getIntensity();

            for(final ProcessedPeak p : mergedPeaks) System.out.print("(" + p.getMass() + "," + (p.getIntensity() / sumIntensity) + ")\t");

        } catch (IOException | UnknownElementException e){
            e.printStackTrace();
        }
    }
     */
}
