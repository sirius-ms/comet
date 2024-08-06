package de.unijena.bioinf.datastructures;

import de.unijena.bioinf.ChemistryBase.chem.MolecularFormula;

public class BuildingBlock {

    /**
     * The monoisotopic mass of this building block in the final molecule (after synthesis).
     */
    private final double mass;

    /**
     * The molecular formula of this building block in the final molecule (after synthesis).
     */
    private final MolecularFormula molecularFormula;

    /**
     * The SMILES string of this building block before synthesis.
     */
    private final String smiles;

    public BuildingBlock(MolecularFormula molecularFormula, String smiles){
        this.molecularFormula = molecularFormula;
        this.mass = molecularFormula.getMass();
        this.smiles = smiles;
    }

    public double getMass(){
        return this.mass;
    }

    public MolecularFormula getMolecularFormula(){
        return this.molecularFormula;
    }

    public String getSmiles(){
        return this.smiles;
    }



}
