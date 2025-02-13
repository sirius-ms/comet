package de.unijena.bioinf.datastructures;

import de.unijena.bioinf.ChemistryBase.chem.PrecursorIonType;
import de.unijena.bioinf.ChemistryBase.ms.Deviation;
import de.unijena.bioinf.ChemistryBase.ms.Ms2Experiment;
import de.unijena.bioinf.ms.annotations.Ms2ExperimentAnnotation;
import lombok.Getter;
import org.jetbrains.annotations.NotNull;

import java.util.*;

@Getter
public class CMLCandidates implements Ms2ExperimentAnnotation, Iterable<CMLCandidate> {

    private final OrderedCombinatorialMoleculeLibrary orderedLibrary;
    private final Ms2Experiment experiment;
    private final Deviation ms1Deviation;
    private final List<PrecursorIonType> fallbackAdducts;
    private final ArrayList<CMLCandidate> candidates;

    public CMLCandidates(OrderedCombinatorialMoleculeLibrary orderedLibrary, Ms2Experiment exp, Deviation ms1Deviation, List<PrecursorIonType> fallbackAdducts) {
        this.orderedLibrary = orderedLibrary;
        this.experiment = exp;
        this.ms1Deviation = ms1Deviation;
        this.fallbackAdducts = fallbackAdducts;
        this.candidates = this.findCandidates();
    }

    /**
     * Returns a list of candidate molecules in the given combinatorial molecule library based on the measured precursor m/z.<br>
     * If the ionization of the precursor molecule is unknown, candidate molecules are retrieved for each provided fallback adduct
     * by performing several binary searches; otherwise, the detected adduct is used for only one binary search.
     *
     * @return an {@link ArrayList} of {@link CMLCandidate} molecules
     */
    private ArrayList<CMLCandidate> findCandidates() {
        final ArrayList<CMLCandidate> candidates = new ArrayList<>();
        final List<PrecursorIonType> consideredAdducts = this.experiment.getPrecursorIonType().isIonizationUnknown() ?
                this.fallbackAdducts :
                Collections.singletonList(this.experiment.getPrecursorIonType());

        for(final PrecursorIonType ionType : consideredAdducts) {
            final double neutralMass = ionType.precursorMassToNeutralMass(this.experiment.getIonMass());

            // Binary search for finding potential candidates:
            int start = 0; // inclusive
            int end = orderedLibrary.getNumOfMolecules(); // exclusive

            while (start < end) {
                final int mid = (start + end) >>> 1;
                final CMLMolecule compound = orderedLibrary.get(mid);
                final double compoundMass = compound.getMass();

                if (this.ms1Deviation.inErrorWindow(neutralMass, compoundMass)) {
                    candidates.add(new CMLCandidate(this.experiment, compound, ionType));

                    // candidate was found
                    // now, look around this compound for more potential candidates
                    if (mid > start) {
                        for (int currentIdx = mid - 1; currentIdx >= start; currentIdx--) {
                            final CMLMolecule currentCompound = orderedLibrary.get(currentIdx);
                            if (this.ms1Deviation.inErrorWindow(neutralMass, currentCompound.getMass())) {
                                candidates.add(new CMLCandidate(this.experiment, currentCompound, ionType));
                            } else {
                                break;
                            }
                        }
                    }
                    if (mid < end - 1) {
                        for (int currentIdx = mid + 1; currentIdx < end; currentIdx++) {
                            final CMLMolecule currentCompound = orderedLibrary.get(currentIdx);
                            if (this.ms1Deviation.inErrorWindow(neutralMass, currentCompound.getMass())) {
                                candidates.add(new CMLCandidate(this.experiment, currentCompound, ionType));
                            } else {
                                break;
                            }
                        }
                    }
                    break;
                } else if (neutralMass < compoundMass) {
                    end = mid;
                } else {
                    start = mid + 1;
                }
            }
        }
        return candidates;
    }

    public List<String> getCandidateNames(){
        return this.candidates.stream().map(CMLCandidate::getName).toList();
    }

    @Override
    public String toString() {
        return "[" + this.orderedLibrary.getName() + "; {"+  String.join(",", this.getCandidateNames()) + "}]";
    }

    @Override
    public @NotNull Iterator<CMLCandidate> iterator() {
        return this.candidates.iterator();
    }
}
