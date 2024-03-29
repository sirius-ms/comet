/*
 *  This file is part of the SIRIUS Software for analyzing MS and MS/MS data
 *
 *  Copyright (C) 2013-2020 Kai Dührkop, Markus Fleischauer, Marcus Ludwig, Martin A. Hoffman, Fleming Kretschmer, Marvin Meusel and Sebastian Böcker,
 *  Chair of Bioinformatics, Friedrich-Schiller University.
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Affero General Public License
 *  as published by the Free Software Foundation; either
 *  version 3 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License along with SIRIUS.  If not, see <https://www.gnu.org/licenses/agpl-3.0.txt>
 */

package de.unijena.bioinf.ms.frontend.subtools.fingerblast;

import de.unijena.bioinf.ChemistryBase.algorithm.scoring.FormulaScore;
import de.unijena.bioinf.ChemistryBase.algorithm.scoring.SScored;
import de.unijena.bioinf.ChemistryBase.fp.ProbabilityFingerprint;
import de.unijena.bioinf.ChemistryBase.fp.Tanimoto;
import de.unijena.bioinf.ChemistryBase.jobs.SiriusJobs;
import de.unijena.bioinf.ChemistryBase.ms.ft.FTree;
import de.unijena.bioinf.canopus.CanopusResult;
import de.unijena.bioinf.chemdb.FingerprintCandidate;
import de.unijena.bioinf.fingerid.CSIPredictor;
import de.unijena.bioinf.fingerid.FingerIdResult;
import de.unijena.bioinf.fingerid.FingerblastJJob;
import de.unijena.bioinf.fingerid.FingerprintResult;
import de.unijena.bioinf.fingerid.blast.FBCandidates;
import de.unijena.bioinf.fingerid.blast.FingerblastResult;
import de.unijena.bioinf.fingerid.predictor_types.PredictorTypeAnnotation;
import de.unijena.bioinf.jjobs.BasicJJob;
import de.unijena.bioinf.jjobs.JJob;
import de.unijena.bioinf.jjobs.JobSubmitter;
import de.unijena.bioinf.jjobs.Partition;
import de.unijena.bioinf.ms.frontend.core.ApplicationCore;
import de.unijena.bioinf.ms.frontend.subtools.InstanceJob;
import de.unijena.bioinf.ms.frontend.utils.PicoUtils;
import de.unijena.bioinf.projectspace.FormulaResult;
import de.unijena.bioinf.projectspace.FormulaScoring;
import de.unijena.bioinf.projectspace.Instance;
import de.unijena.bioinf.projectspace.ProjectSpaceManagers;
import de.unijena.bioinf.rest.NetUtils;
import org.apache.commons.math3.util.Pair;
import org.jetbrains.annotations.NotNull;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * Subtooljob for CSI:FingerID (structure database search)
 * good example for how to create such a job
 */
public class FingerblastSubToolJob extends InstanceJob {

    public FingerblastSubToolJob(JobSubmitter submitter) {
        super(submitter);
        asSCHEDULER();
    }

    @Override
    public boolean isAlreadyComputed(@NotNull Instance inst) {
        return inst.loadCompoundContainer().hasResults() && inst.loadFormulaResults(FBCandidates.class).stream().map(SScored::getCandidate).anyMatch(c -> c.hasAnnotation(FBCandidates.class));
    }

    @Override
    protected void computeAndAnnotateResult(final @NotNull Instance inst) throws Exception {
        List<? extends SScored<FormulaResult, ? extends FormulaScore>> formulaResults =
                inst.loadFormulaResults(FormulaScoring.class, FTree.class, FingerprintResult.class, FBCandidates.class, CanopusResult.class);

        checkForInterruption();

        if (formulaResults == null || formulaResults.isEmpty()) {
            logInfo("Skipping instance \"" + inst.getExperiment().getName() + "\" because there are no trees computed.");
            return;
        }

        checkFingerprintCompatibilityOrThrow();

        checkForInterruption();

        // add CSIClientData to PS if it is not already there
        NetUtils.tryAndWait(() -> ProjectSpaceManagers.writeFingerIdDataIfMissing(inst.getProjectSpaceManager(), ApplicationCore.WEB_API), this::checkForInterruption);

        updateProgress(10);
        checkForInterruption();

        final @NotNull CSIPredictor csi = NetUtils.tryAndWait(() -> (CSIPredictor)
                        ApplicationCore.WEB_API.getStructurePredictor(
                                inst.getExperiment().getAnnotationOrThrow(PredictorTypeAnnotation.class)
                                        .toPredictors(inst.getExperiment().getPrecursorIonType().getCharge()).iterator().next()),
                this::checkForInterruption);

        updateProgress(15);
        checkForInterruption();

        final Map<FormulaResult, FingerIdResult> formulaResultsMap = formulaResults.stream().map(SScored::getCandidate)
                .filter(res -> res.hasAnnotation(FingerprintResult.class))
                .filter(res -> res.hasAnnotation(CanopusResult.class))
                .collect(Collectors.toMap(res -> res, res -> {
                    FingerIdResult idr = new FingerIdResult(res.getAnnotationOrThrow(FTree.class));
                    idr.setAnnotation(FingerprintResult.class, res.getAnnotationOrThrow(FingerprintResult.class));
                    return idr;
                }, (k, v) -> v,  LinkedHashMap::new));


        final Map<FormulaResult, CanopusResult> formulaCanopusResultsMap = formulaResultsMap.keySet().stream()
                .collect(Collectors.toMap(res -> res, res -> res.getAnnotationOrThrow(CanopusResult.class),
                        (k, v) -> v,  LinkedHashMap::new));

        updateProgress(20);
        {
            final FingerblastJJob job = new FingerblastJJob(csi, ApplicationCore.WEB_API, inst.getExperiment(), new ArrayList<>(formulaResultsMap.values()), new ArrayList<>(formulaCanopusResultsMap.values()));

            checkForInterruption();
            // do computation and await results -> objects are already in formulaResultsMap
            submitSubJob(job).awaitResult();
        }

        updateProgress(50);
        checkForInterruption();

        {
            //calculate and annotate tanimoto scores
            List<Pair<ProbabilityFingerprint, FingerprintCandidate>> tanimotoJobs = new ArrayList<>();
            updateProgress(55);
            formulaResultsMap.values().stream().filter(it -> it.hasAnnotation(FingerprintResult.class) && it.hasAnnotation(FingerblastResult.class)).forEach(it -> {
                final ProbabilityFingerprint fp = it.getPredictedFingerprint();
                it.getFingerprintCandidates().stream().map(SScored::getCandidate).forEach(candidate ->
                        tanimotoJobs.add(Pair.create(fp, candidate))
                );
            });

            updateProgress(60);
            checkForInterruption();

            List<BasicJJob<Boolean>> jobs = Partition.ofNumber(tanimotoJobs, 2 * SiriusJobs.getCPUThreads())
                    .stream().map(l -> new BasicJJob<Boolean>(JobType.CPU) {
                        @Override
                        protected Boolean compute() {
                            l.forEach(p -> p.getSecond().setTanimoto(
                                    Tanimoto.nonProbabilisticTanimoto(p.getSecond().getFingerprint(), p.getFirst())));
                            return Boolean.TRUE;
                        }
                    }).collect(Collectors.toList());

            updateProgress(65);
            jobs.forEach(this::submitJob);

            updateProgress(70);
            jobs.forEach(JJob::getResult);

            updateProgress(85);
            checkForInterruption();
        }


        //annotate FingerIdResults to FormulaResult
        inst.saveStructureSearchResults(formulaResultsMap);
        updateProgress(97);

    }

    @Override
    public String getToolName() {
        return PicoUtils.getCommand(FingerblastOptions.class).name();
    }
}
