package de.unijena.bioinf.sirius.plugins;

import de.unijena.bioinf.ChemistryBase.chem.*;
import de.unijena.bioinf.ChemistryBase.chem.utils.ValenceFilter;
import de.unijena.bioinf.ChemistryBase.ms.Deviation;
import de.unijena.bioinf.ChemistryBase.ms.MS2MassDeviation;
import de.unijena.bioinf.ChemistryBase.ms.PossibleAdducts;
import de.unijena.bioinf.ChemistryBase.ms.ft.model.Decomposition;
import de.unijena.bioinf.ChemistryBase.ms.ft.model.FormulaSearchSettings;
import de.unijena.bioinf.ChemistryBase.ms.ft.model.Whiteset;
import de.unijena.bioinf.FragmentationTreeConstruction.computation.SiriusPlugin;
import de.unijena.bioinf.FragmentationTreeConstruction.computation.scoring.DBPairedScorer;
import de.unijena.bioinf.sirius.ProcessedInput;
import de.unijena.bioinf.sirius.ProcessedPeak;
import it.unimi.dsi.fastutil.objects.Object2DoubleOpenHashMap;

import java.io.IOException;
import java.io.ObjectInputStream;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.zip.GZIPInputStream;

public class BottomUpSearch extends SiriusPlugin {

    public final static MolecularFormulaMap MOLECULAR_FORMULA_MAP;
    static {
        MolecularFormulaMap map = null;
        try (final ObjectInputStream OIN =  new ObjectInputStream(new GZIPInputStream(BottomUpSearch.class.getResourceAsStream("/bioformulas.bin.gz")))) {
            map = (MolecularFormulaMap) OIN.readObject();
        } catch (IOException | ClassNotFoundException e) {
            throw new RuntimeException(e);
        }
        MOLECULAR_FORMULA_MAP = map;
    }

    @Override
    public void initializePlugin(PluginInitializer initializer) {
        initializer.addFragmentScorer(new DBPairedScorer());
    }

    private static List<Decomposition> generateDecompositions(ProcessedInput input, boolean addToWhiteset) {
        if (addToWhiteset && input.getAnnotationOrThrow(Whiteset.class).isFinalized()) return null;
        final ValenceFilter filter = new ValenceFilter();

        final Object2DoubleOpenHashMap<Decomposition> weighting = new Object2DoubleOpenHashMap<>();
        final Deviation dev = input.getAnnotation(MS2MassDeviation.class).map(x->x.allowedMassDeviation).orElse(new Deviation(5));
        final Set<IonMode> ionModes = input.getAnnotationOrThrow(PossibleAdducts.class).getIonModes();
        if (!input.getExperimentInformation().getPrecursorIonType().isIonizationUnknown()) { //todo ElementFilter this should be already decided. This is the wrong location to check again. Or has this to do with recalibration?
            ionModes.clear();
            ionModes.add((IonMode)input.getExperimentInformation().getPrecursorIonType().getIonization());
        }
        for (ProcessedPeak peak :input.getMergedPeaks()) {
            final double rootLossMass = input.getExperimentInformation().getIonMass() - peak.getMass();
            final MolecularFormula[] losses = MOLECULAR_FORMULA_MAP.searchMass(rootLossMass, dev);
            for (Ionization ionMode : ionModes) {
                final double mz = ionMode.subtractFromMass(peak.getMass());
                if (mz < 1) continue;
                for (MolecularFormula fragment : MOLECULAR_FORMULA_MAP.searchMass(mz, dev)) {
                    for (MolecularFormula loss : losses) {
                        final MolecularFormula together = fragment.add(loss);
                        if (dev.inErrorWindow(input.getExperimentInformation().getIonMass(), ionMode.addToMass(together.getMass())) && !together.maybeCharged() && filter.isValid(together, ionMode)) {
                            weighting.addTo(new Decomposition(together, ionMode, 0d), Math.max(peak.getRelativeIntensity(),  1e-3));
                        }
                    }
                }
            }
        }
        final HashSet<Decomposition> decompositions = new HashSet<>();
        weighting.forEach((formula, weight)->{
            if (weight >= 0.05) {
                decompositions.add(formula);
            }
        });
        if (addToWhiteset) {
            Set<MolecularFormula> formulas = decompositions.stream().map(Decomposition::getCandidate).collect(Collectors.toSet());
            Whiteset ws = input.getAnnotationOrThrow(Whiteset.class);

            FormulaSearchSettings formulaSearchSettings = input.getAnnotation(FormulaSearchSettings.class, FormulaSearchSettings::bottomUpOnly);
            if (formulaSearchSettings.applyFormulaConstraintsToBottomUp) {
                final FormulaConstraints formulaConstraints = input.getAnnotationOrThrow(FormulaConstraints.class);
                PossibleAdducts possibleAdducts = input.getAnnotationOrThrow(PossibleAdducts.class);
                formulas = Whiteset.filterMeasuredFormulas(formulas, formulaConstraints, possibleAdducts.getAdducts());
            }

            Whiteset bottomUpWs = Whiteset.ofMeasuredFormulas(formulas, BottomUpSearch.class);
            ws = ws.setRequiresBottomUp(false).add(bottomUpWs);
            input.setAnnotation(Whiteset.class, ws);
        }
        return new ArrayList<>(decompositions);
    }

    public static boolean generateDecompositionsAndSaveToWhiteset(ProcessedInput input) {
        return generateDecompositions(input, true) == null ? false : true;
    }

    // TODO: refactor
    @Override
    protected void beforeDecomposing(ProcessedInput input) {
        super.beforeDecomposing(input);

        //run if enabled
        //todo ElementFilter: the following comment relates to the "!input.getExperimentInformation().getPrecursorIonType().isIonizationUnknown()" check -> maybe not necessary anymore with new whiteset workflow?
        // urks, we have to rewrite this horrible recalibration routine. if we recalibrate a tree, we have to fix its molecular formula.
        // same is when the user specifies a molecular formula in the file
        if (!input.getAnnotationOrThrow(Whiteset.class).stillRequiresBottomUpBeAdded()) return; //has already been performed at the isotope pattern analysis stage.

        generateDecompositionsAndSaveToWhiteset(input);

    }


}
