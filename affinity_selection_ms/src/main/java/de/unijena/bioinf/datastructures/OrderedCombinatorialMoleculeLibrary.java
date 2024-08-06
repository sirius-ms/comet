package de.unijena.bioinf.datastructures;

import de.unijena.bioinf.ChemistryBase.chem.MolecularFormula;

import java.util.ArrayList;
import java.util.Comparator;

public class OrderedCombinatorialMoleculeLibrary extends CombinatorialMoleculeLibrary {

    private ArrayList<CMLMolecule> orderedLibrary;

    public OrderedCombinatorialMoleculeLibrary(BuildingBlock[][] buildingBlocks, Scaffold scaffold){
        super(buildingBlocks, scaffold);
        this.orderedLibrary = this.generateMolecules();
    }

    public OrderedCombinatorialMoleculeLibrary(String name, BuildingBlock[][] buildingBlocks, Scaffold scaffold){
        super(name, buildingBlocks, scaffold);
        this.orderedLibrary = this.generateMolecules();
    }

    public static OrderedCombinatorialMoleculeLibrary emptyLibrary(String name){
        return new OrderedCombinatorialMoleculeLibrary(name, new BuildingBlock[0][0], new Scaffold(MolecularFormula.emptyFormula(), null));
    }

    public static OrderedCombinatorialMoleculeLibrary emptyLibrary(){
        return emptyLibrary(null);
    }

    @Override
    public ArrayList<CMLMolecule> generateMolecules(){
        this.orderedLibrary = super.generateMolecules();
        this.orderedLibrary.sort(Comparator.comparingDouble(CMLMolecule::getMass));
        return this.orderedLibrary;
    }

    public CMLMolecule get(int index){
        return this.orderedLibrary.get(index);
    }

}
