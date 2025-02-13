package de.unijena.bioinf.cmlScoring;

import de.unijena.bioinf.ChemistryBase.algorithm.scoring.SScored;
import de.unijena.bioinf.ChemistryBase.algorithm.scoring.Score.DoubleScore;
import de.unijena.bioinf.ChemistryBase.chem.MolecularFormula;
import de.unijena.bioinf.ChemistryBase.chem.PrecursorIonType;
import de.unijena.bioinf.ChemistryBase.ms.Deviation;
import de.unijena.bioinf.ChemistryBase.ms.Ms2Experiment;
import de.unijena.bioinf.datastructures.*;
import de.unijena.bioinf.sirius.ProcessedInput;
import de.unijena.bioinf.sirius.ProcessedPeak;
import de.unijena.bioinf.sirius.Sirius;

import java.util.List;

public class SimplePeakExplainingScorer implements CMLScorer<DoubleScore>{

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
        // todo: the ionization should also match the possible ionization of the candidate molecule (for the former evaluation this wasn't a problem, but in general it can be a problem)
        // todo: if the ionization is unknown, the candidate molecule can maybe explain the measured precursor-peak with e.g. [M+Na]+, or [M+K]+ and so on...(it doesn't have to be [M+H]+)
        // todo: therefore, this method should get as additional parameter the precursor ion type
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

/*
    public static void main(String[] args){
        try{
            File msFile = new File("/home/nils/Dokumente/Bioinformatik_PhD/AS-MS-Project/Data/benzimidazole/ENL161/annotation/annotated_spectra/post_filtered_spectra/manual/ENL161_spectrum_170.ms");
            Ms2Experiment exp = MsIO.readExperimentFromFile(msFile).next();

            File bbsFile = new File("/home/nils/Dokumente/Bioinformatik_PhD/AS-MS-Project/Data/benzimidazole/ENL161/bbs_and_compounds/ENL161_CustomDB_BBs.csv");
            Scaffold scaffold = new Scaffold(MolecularFormula.parse("C8H3N2O"), "");
            BuildingBlock[][] bbs = BuildingBlockReader.readBuildingBlocks(bbsFile);

            CombinatorialMoleculeLibrary cmlLibrary = new CombinatorialMoleculeLibrary(bbs, scaffold);
            List<CMLMolecule> mols = cmlLibrary.generateMolecules();
            HashMap<String,CMLMolecule> name2Mol = new HashMap<>();
            mols.forEach(m -> name2Mol.put(m.getName(), m));

            CMLMolecule molecule = name2Mol.get("2-3-7");

            SimplePeakExplainingScorer scorer = new SimplePeakExplainingScorer(5, 2);
            System.out.println("Scorer:\t" + scorer.score(exp, molecule).getScore()); // expected 0.8551876379690948


            MolecularGraph mol = new MolecularGraph(new SmilesParser(SilentChemObjectBuilder.getInstance()).parseSmiles("CC(C(=O)N)NC(=O)C1=CC2=C(C=C1)N(C(=N2)CCCCO)CCCCCC(=O)O"));
            double[] bbMasses = new double[bbs.length];
            for(int i = 0; i < bbs.length; i++){
                bbMasses[i] = bbs[i][molecule.getBbsIndices()[i]].getMass();
            }
            PrecursorIonType ionType = exp.getPrecursorIonType().isIonizationUnknown() ? PrecursorIonType.fromString("[M+H]+") : exp.getPrecursorIonType();
            BBBarcodeSpectrumPredictor specPred = new BBBarcodeSpectrumPredictor(mol, bbMasses, scaffold.getMass(), 2, ionType);
            SimpleSpectrum predSpectrum = new SimpleSpectrum(specPred.predictSpectrum());

            ProcessedInput processedInput = new Sirius().preprocessForMs2Analysis(exp);
            SimpleMutableSpectrum msrdSpectrum = new SimpleMutableSpectrum();
            List<ProcessedPeak> mergedPeaks = processedInput.getMergedPeaks();
            mergedPeaks.removeLast();

            double sumIntensity = 0d;
            for(ProcessedPeak peak : mergedPeaks)sumIntensity += peak.getRelativeIntensity();
            System.out.println(sumIntensity);

            for(ProcessedPeak peak : mergedPeaks){
                msrdSpectrum.addPeak(peak.getMass(), peak.getRelativeIntensity() / sumIntensity);
            }

            ModifiedCosine spectralAlign = new ModifiedCosine(new Deviation(5));
            System.out.println("Spec_Pred:\t " + spectralAlign.score(new SimpleSpectrum(msrdSpectrum), predSpectrum, ionType.neutralMassToPrecursorMass(mol.getFormula().getMass()), ionType.neutralMassToPrecursorMass(molecule.getMass())));


        }catch(IOException | UnknownElementException e){
            e.printStackTrace();
        } catch (InvalidSmilesException e) {
            throw new RuntimeException(e);
        }
    }
 */
}
