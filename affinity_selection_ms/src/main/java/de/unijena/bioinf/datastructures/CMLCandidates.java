package de.unijena.bioinf.datastructures;

import de.unijena.bioinf.ChemistryBase.chem.PrecursorIonType;
import de.unijena.bioinf.ChemistryBase.ms.Deviation;
import de.unijena.bioinf.ChemistryBase.ms.Ms2Experiment;
import de.unijena.bioinf.ms.annotations.Ms2ExperimentAnnotation;
import lombok.Getter;

import java.util.*;

@Getter
public class CMLCandidates implements Ms2ExperimentAnnotation {

    private final OrderedCombinatorialMoleculeLibrary orderedLibrary;
    private final Ms2Experiment experiment;
    private final double ms1Deviation;
    private final ArrayList<CMLMolecule> candidates;

    public CMLCandidates(OrderedCombinatorialMoleculeLibrary orderedLibrary, Ms2Experiment exp, double ms1Dev){
        this.orderedLibrary = orderedLibrary;
        this.experiment = exp;
        this.ms1Deviation = ms1Dev;
        this.candidates = this.findCandidates(exp, ms1Dev);
    }

    public ArrayList<CMLMolecule> findCandidates(Ms2Experiment exp, double ms1Dev) {
        final PrecursorIonType ionType = exp.getPrecursorIonType().isIonizationUnknown() ? PrecursorIonType.fromString("[M+H]+") : exp.getPrecursorIonType();
        final double neutralMass = ionType.precursorMassToNeutralMass(exp.getIonMass());
        final Deviation deviation = new Deviation(ms1Dev);

        final ArrayList<CMLMolecule> candidates = new ArrayList<>();

        // Binary search for finding potential candidates:
        int start = 0; // inclusive
        int end = orderedLibrary.getNumOfMolecules(); // exclusive

        while(start < end){
            final int mid = (start + end) >>> 1;
            final CMLMolecule compound = orderedLibrary.get(mid);
            final double compoundMass = compound.getMass();

            if(deviation.inErrorWindow(neutralMass, compoundMass)){
                candidates.add(compound);

                // candidate was found
                // now, look around this compound for more potential candidates
                if(mid > start){
                    for(int currentIdx = mid-1; currentIdx >= start; currentIdx--){
                        final CMLMolecule currentCompound = orderedLibrary.get(currentIdx);
                        if(deviation.inErrorWindow(neutralMass, currentCompound.getMass())){
                            candidates.add(currentCompound);
                        }else{
                            break;
                        }
                    }
                }
                if(mid < end-1){
                    for(int currentIdx = mid+1; currentIdx < end; currentIdx++){
                        final CMLMolecule currentCompound = orderedLibrary.get(currentIdx);
                        if(deviation.inErrorWindow(neutralMass, currentCompound.getMass())){
                            candidates.add(currentCompound);
                        }else{
                            break;
                        }
                    }
                }
                break;
            }else if(neutralMass < compoundMass){
                end = mid;
            }else{
                start = mid+1;
            }
        }
        return candidates;
    }

    public List<String> getCandidateNames(){
        return this.candidates.stream().map(CMLMolecule::getName).toList();
    }

    @Override
    public String toString() {
        return "[" + this.orderedLibrary.getName() + "; {"+  String.join(",", this.getCandidateNames()) + "}]";
    }
}
