package de.unijena.bioinf.cmlSpectrumPrediction;

import de.unijena.bioinf.ChemistryBase.chem.MolecularFormula;
import de.unijena.bioinf.ChemistryBase.chem.PrecursorIonType;
import de.unijena.bioinf.ChemistryBase.ms.*;
import de.unijena.bioinf.ChemistryBase.ms.utils.SimpleMutableSpectrum;
import de.unijena.bioinf.fragmenter.CombinatorialFragment;
import de.unijena.bioinf.fragmenter.MolecularGraph;
import de.unijena.bioinf.projectspace.*;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;
import java.util.Iterator;

public class BBBarcodeSpectrumPredictor extends AbstractMs2SpectrumPredictor<Peak>{

    private final MolecularGraph molecule;
    private final double[] bbMasses;
    private final double scaffoldMass;
    private final int hydrogenShifts;
    private HashMap<BitSet, ArrayList<Peak>> bitset2Peaks;

    public BBBarcodeSpectrumPredictor(MolecularGraph molecule, double[] bbMasses, double scaffoldMass, int hydrogenShifts, PrecursorIonType precursorIonType){
        super(null, precursorIonType);
        this.molecule = molecule;
        this.bbMasses = bbMasses;
        this.scaffoldMass = scaffoldMass;
        this.hydrogenShifts = hydrogenShifts;
    }

    @Override
    public Ms2Spectrum<Peak> predictSpectrum() {
        this.bitset2Peaks = new HashMap<>();
        final double precursorMz = this.precursorIonType.neutralMassToPrecursorMass(this.molecule.getFormula().getMass());
        final SimpleMutableSpectrum spec = new SimpleMutableSpectrum((int) (Math.pow(2, this.bbMasses.length)-1+this.bbMasses.length) * (2*this.hydrogenShifts+1)); // no precursor peak

        BitSet bbCombination = new BitSet(this.bbMasses.length+1);
        while(bbCombination.cardinality() < this.bbMasses.length){
            double fragmentMass = this.scaffoldMass;
            for(int bbIdx = bbCombination.nextSetBit(0); bbIdx >= 0; bbIdx = bbCombination.nextSetBit(bbIdx+1)){
                fragmentMass += this.bbMasses[bbIdx];
            }
            this.storeFragmentPeaks(fragmentMass, spec, bbCombination);
            bbCombination = this.incrementBitSet(bbCombination);
        }

        bbCombination.clear();
        bbCombination.set(this.bbMasses.length);
        for(int i = 0; i < this.bbMasses.length; i++){
            double bbMass = this.bbMasses[i];
            BitSet clonedBBCombination = (BitSet) bbCombination.clone();
            clonedBBCombination.set(i);
            this.storeFragmentPeaks(bbMass, spec, clonedBBCombination);
        }

        this.spectrum = new MutableMs2Spectrum(spec, precursorMz, null, 2);
        return this.spectrum;
    }

    private void storeFragmentPeaks(double fragmentNeutralMass, SimpleMutableSpectrum spec, BitSet bbCombination){
        final double hydrogenMass = MolecularFormula.getHydrogen().getMass();
        final double fragmentMz = this.precursorIonType.neutralMassToPrecursorMass(fragmentNeutralMass);

        for(int h = -this.hydrogenShifts; h <= this.hydrogenShifts; h++) {
            final double shiftedFragmentMz = fragmentMz + h * hydrogenMass;
            final SimplePeak peak = new SimplePeak(shiftedFragmentMz, 1d);
            spec.addPeak(peak);
            this.peak2fragment.put(peak, new CombinatorialFragment(this.molecule, new BitSet(), new BitSet()));
            this.bitset2Peaks.computeIfAbsent(bbCombination, k -> new ArrayList<>()).add(peak);
        }
    }

    private BitSet incrementBitSet(BitSet bitSet){
        BitSet newBitSet = (BitSet) bitSet.clone();
        int i = 0;
        while(newBitSet.get(i)){
            newBitSet.set(i, false);
            i++;
        }
        newBitSet.set(i, true);
        return newBitSet;
    }

    @Override
    public MolecularGraph getPrecursorMolecule(){
        return this.molecule;
    }

    public HashMap<BitSet, ArrayList<Peak>> getBitset2Peaks(){
        return this.bitset2Peaks;
    }
}
