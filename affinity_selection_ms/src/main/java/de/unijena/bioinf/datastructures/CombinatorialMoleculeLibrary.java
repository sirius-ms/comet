package de.unijena.bioinf.datastructures;

import de.unijena.bioinf.ChemistryBase.chem.MolecularFormula;
import lombok.Getter;

import java.util.ArrayList;

@Getter
public class CombinatorialMoleculeLibrary {

    private final int numBBs;

    private final int numOfMolecules;

    private final BuildingBlock[][] buildingBlocks;

    private final Scaffold scaffold;

    private String name;

    public CombinatorialMoleculeLibrary(BuildingBlock[][] buildingBlocks, Scaffold scaffold) {
        this.numBBs = buildingBlocks.length;
        this.buildingBlocks = buildingBlocks;
        this.scaffold = scaffold;
        this.numOfMolecules = this.numberOfMolecules();
    }

    public CombinatorialMoleculeLibrary(String name, BuildingBlock[][] buildingBlocks, Scaffold scaffold) {
        this(buildingBlocks, scaffold);
        this.name = name;
    }

    public static CombinatorialMoleculeLibrary emptyLibrary(String name){
        return new CombinatorialMoleculeLibrary(name, new BuildingBlock[0][0], new Scaffold(MolecularFormula.emptyFormula(), null));
    }

    public static CombinatorialMoleculeLibrary emptyLibrary(){
        return emptyLibrary(null);
    }

    private int numberOfMolecules(){
        if(this.numBBs == 0) return 0;
        int num = 1;
        for(int pos = 0; pos < this.numBBs; pos++) num = num * this.buildingBlocks[pos].length;
        return num;
    }

    public ArrayList<CMLMolecule> generateMolecules(){
        if(this.numOfMolecules == 0) return new ArrayList<>();

        final ArrayList<CMLMolecule> molecules = new ArrayList<>(this.numOfMolecules);
        int[] currentCombination = new int[this.numBBs];

        while(true){
            // Generate a molecule contained in this library based on 'currentCombination':
            final CMLMolecule molecule = new CMLMolecule(currentCombination.clone(), this);
            molecules.add(molecule);

            // Determine the next possible combination of building blocks.
            // If there is none, stop the whole generation process.
            int k = this.numBBs-1;
            while(k >= 0 && currentCombination[k] == this.buildingBlocks[k].length - 1){
                currentCombination[k] = 0;
                k--;
            }
            if(k >= 0){
                currentCombination[k]++;
            }else{
                break;
            }
        }

        return molecules;
    }
}
