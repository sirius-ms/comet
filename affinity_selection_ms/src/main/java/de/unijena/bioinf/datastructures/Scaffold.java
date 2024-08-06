package de.unijena.bioinf.datastructures;

import de.unijena.bioinf.ChemistryBase.chem.MolecularFormula;

public class Scaffold {

    /**
     * The SMILES string of the scaffold structure before synthesis.
     */
    private final String smiles;

    /**
     * The molecular formula of the scaffold structure in the final molecule structure.
     */
    private final MolecularFormula molecularFormula;

    /**
     * The monoisotopic mass of the scaffold structure in the final molecule structure.
     */
    private final double mass;

    public Scaffold(MolecularFormula molecularFormula, String smiles){
        this.molecularFormula = molecularFormula;
        this.mass = molecularFormula.getMass();
        this.smiles = smiles;
    }


    public String getSmiles() {
        return smiles;
    }

    public MolecularFormula getMolecularFormula() {
        return molecularFormula;
    }

    public double getMass() {
        return mass;
    }
}
