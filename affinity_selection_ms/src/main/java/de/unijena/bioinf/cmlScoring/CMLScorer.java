package de.unijena.bioinf.cmlScoring;

import de.unijena.bioinf.ChemistryBase.algorithm.scoring.SScored;
import de.unijena.bioinf.ChemistryBase.algorithm.scoring.Score;
import de.unijena.bioinf.ChemistryBase.ms.Ms2Experiment;
import de.unijena.bioinf.datastructures.CMLMolecule;


public interface CMLScorer<S extends Score> {

    SScored<CMLMolecule, S> score(Ms2Experiment exp, CMLMolecule molecule);

}
