package de.unijena.bioinf.datastructures;


import lombok.Getter;
import org.jetbrains.annotations.NotNull;

import java.util.*;

@Getter
public class CMLMolecule {

    /**
     * The identifier or name of this compound.<br>
     * It encodes which building blocks are contained in this molecule.
     * E.g. "1-6-2" for a compound which contains the building block 1 at pos. 1, the bb. 6 at pos. 2 and the bb. 2 at pos. 3.
     */
    private final String name;

    /**
     * The SMILES string of this molecule.
     */
    private String smiles;

    /**
     * The neutral monoisotopic mass of this molecule.
     */
    private final double mass;

    /**
     * The indices of the building blocks contained in this molecule.
     */
    private final int[] bbsIndices;

    /**
     * A pointer back to the library this molecule belongs to.
     */
    private final CombinatorialMoleculeLibrary library;


    public CMLMolecule(String smiles, int[] bbsIndices, CombinatorialMoleculeLibrary library){
        this.smiles = smiles;
        this.bbsIndices = bbsIndices;
        this.library = library;
        StringBuilder compoundName = new StringBuilder(Integer.toString(bbsIndices[0] + 1));
        for(int i = 1; i < bbsIndices.length; i++) compoundName.append("-").append(bbsIndices[i] + 1);
        this.name = compoundName.toString();

        BuildingBlock[][] allBuildingBlocks = library.getBuildingBlocks();
        double mass = library.getScaffold().getMass();
        for(int pos = 0; pos < bbsIndices.length; pos++){
            final int bbIdx = bbsIndices[pos];
            mass += allBuildingBlocks[pos][bbIdx].getMass();
        }
        this.mass = mass;
    }

    public CMLMolecule(int[] bbsIndices, CombinatorialMoleculeLibrary library){
        this(null, bbsIndices, library);
    }

    public BuildingBlock[] getBuildingBlocks(){
        final BuildingBlock[] buildingBlocks = new BuildingBlock[this.bbsIndices.length];
        final BuildingBlock[][] allBuildingBlocks = this.library.getBuildingBlocks();
        for(int pos = 0; pos < this.bbsIndices.length; pos++){
            buildingBlocks[pos] = allBuildingBlocks[pos][this.bbsIndices[pos]];
        }
        return buildingBlocks;
    }

    public List<BBFragment> createAllBBFragments(){
        final List<BBFragment> fragments = new ArrayList<>(2^this.bbsIndices.length - 2);

        // 1. Generate all possible fragments containing the scaffold - without the molecule or the scaffold itself:
        BitSet presentBBs = new BitSet(this.bbsIndices.length);
        presentBBs.set(0, true); // we don't want to create a fragment which is just the scaffold
        while(presentBBs.cardinality() < this.bbsIndices.length){
            final BBFragment frag = new BBFragment(this, presentBBs, true);
            fragments.add(frag);
            presentBBs = this.nextCombination(presentBBs);
        }

        // 2. Generate all fragments consisting of only one building block:
        for(int pos = 0; pos < this.bbsIndices.length; pos++){
            final BitSet singleBBBitSet = new BitSet(this.bbsIndices.length);
            singleBBBitSet.set(pos, true);
            final BBFragment frag = new BBFragment(this, singleBBBitSet, false);
            fragments.add(frag);
        }

        return fragments;
    }

    public List<BBFragment> createSpecificBBFragments(Collection<String> fragmentTypes){
        return this.createAllBBFragments().stream().filter(f -> fragmentTypes.contains(f.getFragmentTypeString())).toList();
    }

    private @NotNull BitSet nextCombination(@NotNull BitSet presentBBs){
        BitSet newPresentBBs = (BitSet) presentBBs.clone();
        int idx = 0;
        while(newPresentBBs.get(idx)){
            newPresentBBs.set(idx, false);
            idx++;
        }
        newPresentBBs.set(idx, true);
        return newPresentBBs;
    }

}
