package de.unijena.bioinf.cmlFragmentation;

import de.unijena.bioinf.fragmenter.*;

import java.util.ArrayList;
import java.util.List;

public abstract class AbstractFragmentationPredictor implements FragmentationPredictor{

    protected final MolecularGraph molecule;
    protected final ArrayList<CombinatorialFragment> fragments; // fragments of 'molecule' without 'molecule'
    protected CombinatorialGraph graph;

    public AbstractFragmentationPredictor(MolecularGraph molecule){
        this.molecule = molecule;
        this.fragments = new ArrayList<>();
    }

    @Override
    public List<CombinatorialFragment> getFragments(){
        return this.fragments;
    }

    public List<CombinatorialFragment> getFragmentsWithPrecursor(){
        ArrayList<CombinatorialFragment> allFragments = new ArrayList<>(this.fragments);
        allFragments.add(this.molecule.asFragment());
        return allFragments;
    }

    @Override
    public CombinatorialGraph getFragmentationGraph(){
        return this.graph;
    }

    @Override
    public MolecularGraph getMolecule(){
        return this.molecule;
    }

}
