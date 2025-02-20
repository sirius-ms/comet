/*
 *
 *  This file is part of the SIRIUS library for analyzing MS and MS/MS data
 *
 *  Copyright (C) 2013-2020 Kai Dührkop, Markus Fleischauer, Marcus Ludwig, Martin A. Hoffman and Sebastian Böcker,
 *  Chair of Bioinformatics, Friedrich-Schilller University.
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 3 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License along with SIRIUS. If not, see <https://www.gnu.org/licenses/lgpl-3.0.txt>
 */

package de.unijena.bioinf.sirius;

import de.unijena.bioinf.ChemistryBase.algorithm.scoring.FormulaScore;
import de.unijena.bioinf.ChemistryBase.math.Statistics;
import de.unijena.bioinf.ChemistryBase.ms.AnnotatedPeak;
import de.unijena.bioinf.ChemistryBase.ms.Deviation;
import de.unijena.bioinf.ChemistryBase.ms.ft.*;
import de.unijena.bioinf.FragmentationTreeConstruction.computation.FragmentationPatternAnalysis;
import de.unijena.bioinf.sirius.plugins.IsotopePatternInMs1Plugin;
import de.unijena.bioinf.sirius.scores.IsotopeScore;
import de.unijena.bioinf.sirius.scores.SiriusScore;
import de.unijena.bioinf.sirius.scores.TreeScore;
import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import org.jetbrains.annotations.NotNull;
import org.slf4j.LoggerFactory;

import java.util.HashSet;
import java.util.List;
import java.util.Optional;
import java.util.Set;

public class FTreeMetricsHelper {
    protected final FTree tree;
    protected final FragmentAnnotation<Score> fragmentScoring;
    protected final LossAnnotation<Score> lossScoring;
    protected final Fragment measuredIonRoot;

    public FTreeMetricsHelper(FTree tree) {
        this.tree = tree;
        this.fragmentScoring = tree.getOrCreateFragmentAnnotation(Score.class);
        this.lossScoring = tree.getOrCreateLossAnnotation(Score.class);
        this.measuredIonRoot = IonTreeUtils.getMeasuredIonRoot(tree);
    }

    public Optional<TreeStatistics> asTreeStats() {
        return tree.getAnnotation(TreeStatistics.class);
    }

    public double getSiriusScore() {
        return getSiriusScore(tree);
    }

    public double getTreeScore() {
        return getSiriusScore() - getIsotopeMs1Score();
    }

    public double getRootScore() {
        return getRootScore(tree);
    }

    public double getExplainedPeaksRatio() {
        return getExplainedPeaksRatio(tree);
    }

    public int getNumOfExplainedPeaks() {
        return getNumOfExplainedPeaks(tree);
    }

    public double getExplainedIntensityRatio() {
        return getExplainedIntensityRatio(tree);
    }

    public int getNumberOfExplainablePeaks() {
        return getNumberOfExplainablePeaks(tree);
    }

    public double getIsotopeMs1Score() {
            Score score = fragmentScoring.get(measuredIonRoot);
            return score.get(FragmentationPatternAnalysis.getScoringMethodName(IsotopePatternInMs1Plugin.Ms1IsotopePatternScorer.class));
    }

    public double getBeautificationPenalty() {
        return fragmentScoring.get(measuredIonRoot).get(Beautified.PENALTY_KEY);
    }

    public double getRecalibrationPenalty() {
        return fragmentScoring.get(measuredIonRoot).get(Recalibrated.PENALTY_KEY);
    }

    /**
     *
     * @param useAbsoluteValues absolute values (Betrag)
     * @return
     */
    private Deviation getMedianMassDeviation(boolean useAbsoluteValues) {
        DoubleArrayList ppms = new DoubleArrayList(), mzs = new DoubleArrayList();
        for (Fragment f : tree) {
            final Deviation dev = tree.getMassError(f);
            if (dev != Deviation.NULL_DEVIATION) {
                ppms.add(useAbsoluteValues ? Math.abs(dev.getPpm()) : dev.getPpm());
                mzs.add(useAbsoluteValues ? Math.abs(dev.getAbsolute()) : dev.getAbsolute());
            }
        }

        return new Deviation(
                ppms.isEmpty() ? Double.NaN : Statistics.median(ppms),
                mzs.isEmpty() ? Double.NaN : Statistics.median(mzs)
        );

    }

    public Deviation getMedianMassDeviation() {
        return getMedianMassDeviation(false);
    }

    public Deviation getMedianAbsoluteMassDeviation() {
        return getMedianMassDeviation(true);
    }


    // static helper methods

    public static double getSiriusScore(FTree tree) {
        return tree.getTreeWeight();
    }

    public static double getRootScore(@NotNull FTree tree) {
        return tree.getRootScore();
    }

    public static double getExplainedPeaksRatio(@NotNull FTree tree) {
        return tree.getAnnotationOrThrow(TreeStatistics.class).getRatioOfExplainedPeaks();
    }

    public static double getExplainedIntensityRatio(@NotNull FTree tree) {
        return tree.getAnnotationOrThrow(TreeStatistics.class).getExplainedIntensity();
    }

    public static int getNumberOfExplainablePeaks(@NotNull FTree tree) {
        final int numberOfExplainedPeaks = getNumOfExplainedPeaks(tree);
        if (numberOfExplainedPeaks==0) return 0;
        return (int)Math.ceil(getNumOfExplainedPeaks(tree) / getExplainedPeaksRatio(tree));
    }

    /**
     * number of peaks explained by FTree, ignoring the root
     * @param tree
     * @return
     */
    public static int getNumOfExplainedPeaks(@NotNull FTree tree) {
        // only count fragments that belong to a peak
        FragmentAnnotation<AnnotatedPeak> peakAno = tree.getFragmentAnnotationOrNull(AnnotatedPeak.class);
        LossAnnotation<LossType> lossAno = tree.getLossAnnotationOrNull(LossType.class);
        if (peakAno==null) {
            LoggerFactory.getLogger(FTreeMetricsHelper.class).error("Fragmentation tree has no peak annotation.");
            return 0; // should never happen
        }

        return (int)tree.getFragmentsWithoutRoot().stream().filter(x->peakAno.get(x).isMeasured() && (lossAno==null || lossAno.get(x.getIncomingEdge(), LossType::regular)!=LossType.insource())).count();
    }

    public static Set<FormulaScore> getScoresFromTree(@NotNull final FTree tree) {
        final FTreeMetricsHelper helper = new FTreeMetricsHelper(tree);
        final Set<FormulaScore> scores = new HashSet<>(3);
        scores.add(new SiriusScore(helper.getSiriusScore()));
        try {
            scores.add(new IsotopeScore(helper.getIsotopeMs1Score()));
            scores.add(new TreeScore(helper.getTreeScore()));
        } catch (Throwable e) {
            //todo remove debug stuff
            System.out.println("DEBUG: Something with this tree is wrong? Cannot calculate isotope score" + tree.getRoot().getFormula().toString());
            e.printStackTrace();
        }

        return scores;
    }

    public static void computeNormalizedSiriusScores(List<IdentificationResult> idResults, double maxTreeWeight, double remainingCandidatesTreeWeightExpSumEstimate) {
        if (idResults == null || idResults.isEmpty()) return;

        assert maxTreeWeight == idResults.stream().mapToDouble(r -> r.getTree().getTreeWeight()).max().orElse(0);

        double resultsTreeWeightExpSumEstimate = idResults.stream().mapToDouble(r ->Math.exp(r.getTree().getTreeWeight() - maxTreeWeight)).sum();

        final double normalizationFactor = resultsTreeWeightExpSumEstimate + remainingCandidatesTreeWeightExpSumEstimate;
        idResults.forEach(r -> r.setNormalizedScore(Math.exp(r.getTree().getTreeWeight() - maxTreeWeight) / normalizationFactor));
    }
}
