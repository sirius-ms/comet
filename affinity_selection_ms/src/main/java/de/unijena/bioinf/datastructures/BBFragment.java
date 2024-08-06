package de.unijena.bioinf.datastructures;

import lombok.Getter;

import java.util.BitSet;

@Getter
public class BBFragment {

    private final CMLMolecule molecule;
    private final BitSet presentBBs;
    private final boolean presentScaffold;
    private final double mass;

    public BBFragment(CMLMolecule molecule, BitSet presentBBs, boolean presentScaffold){
        this.molecule = molecule;
        this.presentBBs = presentBBs;
        this.presentScaffold = presentScaffold;

        final BuildingBlock[] buildingBlocks = molecule.getBuildingBlocks();
        double mass = presentScaffold ? molecule.getLibrary().getScaffold().getMass() : 0d;
        for(int i = presentBBs.nextSetBit(0); i >= 0; i = presentBBs.nextSetBit(i+1)){
            mass += buildingBlocks[i].getMass();
        }
        this.mass = mass;
    }

    public String getFragmentTypeString(){
        if(this.presentBBs.isEmpty()){
            return this.presentScaffold ? "S" : "";
        }
        if(this.presentBBs.cardinality() == 1){
            final int bbIdx = this.presentBBs.nextSetBit(0);
            return this.presentScaffold ? "S[" + bbIdx + "]" : Integer.toString(bbIdx);
        }

        final StringBuilder strBuilder = this.presentScaffold ? new StringBuilder("S[") : new StringBuilder("["); // the scaffold should always be present if this fragments contains several building blocks!!!
        int k = this.presentBBs.nextSetBit(0);
        strBuilder.append(k);
        for(int i = this.presentBBs.nextSetBit(k+1); i >= 0; i = presentBBs.nextSetBit(i+1)){
            strBuilder.append(";").append(i);
        }
        strBuilder.append("]");
        return strBuilder.toString();
    }

    @Override
    public String toString(){
        return "(" + this.molecule.getName() + ", " + this.getFragmentTypeString() + ")";
    }



}
