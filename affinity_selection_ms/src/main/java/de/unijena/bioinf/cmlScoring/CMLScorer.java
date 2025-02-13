package de.unijena.bioinf.cmlScoring;

import de.unijena.bioinf.ChemistryBase.algorithm.scoring.SScored;
import de.unijena.bioinf.ChemistryBase.algorithm.scoring.Score;
import de.unijena.bioinf.ChemistryBase.ms.Ms2Experiment;
import de.unijena.bioinf.datastructures.CMLMolecule;


public interface CMLScorer<S extends Score> {

    // todo this method should possible get the precursorIonType with which the molecule can explain the precursor peak of the Ms2Experiment
    // e.g. the molecule is only a candidate because it can explain the precursor peak with the ionization of [M+Na]+ and this wouldn't be detected if
    // 'exp' doesn't have any information about the ionization
    SScored<CMLMolecule, S> score(Ms2Experiment exp, CMLMolecule molecule);

}
