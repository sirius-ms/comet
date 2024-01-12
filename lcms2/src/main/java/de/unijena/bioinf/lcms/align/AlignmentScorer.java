package de.unijena.bioinf.lcms.align;

public class AlignmentScorer {

    private final double intensityAwareness;

    public AlignmentScorer(double intensityAwareness) {
        this.intensityAwareness = intensityAwareness;
    }

    public static AlignmentScorer expectSimilarIntensity() {
        return new AlignmentScorer(1.5);
    }
    public static AlignmentScorer intensityMayBeDifferent() {
        return new AlignmentScorer(4);
    }

    public double score(AlignmentStatistics stats, MoI left, MoI right) {

        final double rtDiff = left.getRetentionTime()-right.getMz();
        final double mzDiff = left.getMz() - right.getMz();
        final double intensityDelta = Math.log((left.getIntensity()+0.05)/(right.getIntensity()+0.05));


        return gaussianRbfLog(rtDiff, stats.expectedRetentionTimeDeviation) +
                gaussianRbfLog(mzDiff, stats.expectedMassDeviationBetweenSamples.absoluteFor(Math.max(left.getMz(),right.getMz()))) +
                gaussianRbfLog(intensityDelta, 0.1);

    }

    private static double gaussianRbfLog(double delta, double std) {
        return -(delta*delta)/(std*std);
    }

}
