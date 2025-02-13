package de.unijena.bioinf.datastructures;

import de.unijena.bioinf.ChemistryBase.chem.PrecursorIonType;
import de.unijena.bioinf.ChemistryBase.ms.Ms2Experiment;

public record CMLCandidate(Ms2Experiment experiment, CMLMolecule candidate, PrecursorIonType adduct) {

    public String getName(){
        return "(" + candidate.getName() + ";" + adduct.toString() + ")";
    }

}
