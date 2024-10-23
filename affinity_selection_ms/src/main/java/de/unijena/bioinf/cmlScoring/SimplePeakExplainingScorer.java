package de.unijena.bioinf.cmlScoring;

import de.unijena.bioinf.ChemistryBase.algorithm.scoring.SScored;
import de.unijena.bioinf.ChemistryBase.algorithm.scoring.Score.DoubleScore;
import de.unijena.bioinf.ChemistryBase.chem.MolecularFormula;
import de.unijena.bioinf.ChemistryBase.chem.PrecursorIonType;
import de.unijena.bioinf.ChemistryBase.chem.utils.UnknownElementException;
import de.unijena.bioinf.ChemistryBase.ms.Deviation;
import de.unijena.bioinf.ChemistryBase.ms.Ms2Experiment;
import de.unijena.bioinf.babelms.MsIO;
import de.unijena.bioinf.datastructures.*;
import de.unijena.bioinf.io.BuildingBlockReader;
import de.unijena.bioinf.sirius.ProcessedInput;
import de.unijena.bioinf.sirius.ProcessedPeak;
import de.unijena.bioinf.sirius.Sirius;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;

public class SimplePeakExplainingScorer implements CMLScorer{

    private final Deviation massDeviation;
    private final int allowedHydrogenShifts;

    public SimplePeakExplainingScorer(double ppm, int allowedHydrogenShifts) {
        this.massDeviation = new Deviation(ppm);
        this.allowedHydrogenShifts = allowedHydrogenShifts;
    }

    @Override
    public SScored<CMLMolecule, DoubleScore> score(Ms2Experiment exp, CMLMolecule molecule) {
        // PREPROCESSING:
        // 1. Preprocessing of the Ms2Experiment (i.e. merging the fragment peaks of several measurements etc.):
        ProcessedInput processedInput = new Sirius().preprocessForMs2Analysis(exp);
        List<ProcessedPeak> fragmentPeaks = processedInput.getMergedPeaks();
        fragmentPeaks.removeLast(); // remove precursor/parent peak

        // 2. Create all building block fragments of 'molecule':
        List<BBFragment> bbFragments = molecule.createAllBBFragments();

        // 3. Additional objects/variables:
        PrecursorIonType ionType = exp.getPrecursorIonType().isIonizationUnknown() ? PrecursorIonType.fromString("[M+H]+") : exp.getPrecursorIonType();
        double hydrogenMass = MolecularFormula.getHydrogen().getMass();

        // MAPPING OF THE BUILDING BLOCK FRAGMENTS ONTO THE FRAGMENT PEAKS:
        double explainedIntensity = 0d;
        for(ProcessedPeak peak : fragmentPeaks) {
            double peakMz = peak.getMass();
            boolean isMatched = false;

            int idx = 0;
            while(!isMatched && idx < bbFragments.size()) {
                BBFragment fragment = bbFragments.get(idx);
                double fragmentMz = ionType.neutralMassToPrecursorMass(fragment.getMass());
                for(int h = -allowedHydrogenShifts; h <= allowedHydrogenShifts; h++) {
                    double shiftedFragmentMz = fragmentMz + h * hydrogenMass;

                    double allowedDeviation = massDeviation.absoluteFor(Math.min(peakMz, shiftedFragmentMz));
                    double observedDeviation = Math.abs(peakMz - shiftedFragmentMz);
                    if(observedDeviation <= allowedDeviation) {
                        explainedIntensity += peak.getRelativeIntensity();
                        isMatched = true;
                        break;
                    }
                }
                idx++;
            }
        }

        double summedTotalIntensity = 0d;
        for(ProcessedPeak peak : fragmentPeaks) summedTotalIntensity += peak.getRelativeIntensity();
        DoubleScore score = new DoubleScore(explainedIntensity / summedTotalIntensity);

        return new SScored<>(molecule, score);
    }

    public static void main(String[] args){
        try{
            File msFile = new File("/home/nils/Dokumente/Bioinformatik_PhD/AS-MS-Project/Data/benzimidazole/ENL161/annotation/annotated_spectra/post_filtered_spectra/manual/ENL161_spectrum_42.ms");
            Ms2Experiment exp = MsIO.readExperimentFromFile(msFile).next();

            File bbsFile = new File("/home/nils/Dokumente/Bioinformatik_PhD/AS-MS-Project/Data/benzimidazole/ENL161/bbs_and_compounds/ENL161_CustomDB_BBs.csv");
            Scaffold scaffold = new Scaffold(MolecularFormula.parse("C8H3N2O"), "");
            BuildingBlock[][] bbs = BuildingBlockReader.readBuildingBlocks(bbsFile);

            CombinatorialMoleculeLibrary cmlLibrary = new CombinatorialMoleculeLibrary(bbs, scaffold);
            List<CMLMolecule> mols = cmlLibrary.generateMolecules();
            HashMap<String,CMLMolecule> name2Mol = new HashMap<>();
            mols.forEach(m -> name2Mol.put(m.getName(), m));

            CMLMolecule molecule = name2Mol.get("2-10-7");

            SimplePeakExplainingScorer scorer = new SimplePeakExplainingScorer(5, 2);
            System.out.println(scorer.score(exp, molecule)); // expected 0.8551876379690948

        }catch(IOException | UnknownElementException e){
            e.printStackTrace();
        }
    }
}
