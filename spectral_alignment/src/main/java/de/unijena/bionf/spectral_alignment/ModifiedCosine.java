package de.unijena.bionf.spectral_alignment;

import de.unijena.bioinf.ChemistryBase.ms.Deviation;
import de.unijena.bioinf.ChemistryBase.ms.Peak;
import de.unijena.bioinf.ChemistryBase.ms.utils.OrderedSpectrum;
import de.unijena.bioinf.ChemistryBase.ms.utils.SimpleMutableSpectrum;
import de.unijena.bioinf.ChemistryBase.ms.utils.SimpleSpectrum;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;

import java.util.Arrays;
import java.util.BitSet;

/**
 * This algorithm requires that there is at most one pair of peaks (u,v) where the m/z of u
 * and v are within the allowed mass tolerance.
 */
public class ModifiedCosine implements SpectralAlignmentScorer {

    public static void main(String[] args) {

        final SimpleMutableSpectrum A = new SimpleMutableSpectrum();
        A.addPeak(1, 1);
        A.addPeak(5, 1);
        A.addPeak(8, 1);
        A.addPeak(15, 0.5);
        final SimpleMutableSpectrum B = new SimpleMutableSpectrum();
        B.addPeak(1, 1);
        B.addPeak(8, 1);
        B.addPeak(15, 1);
        B.addPeak(18, 1);
        B.addPeak(20, 1);
        B.addPeak(25, 1);

        final ModifiedCosine modifiedCosine = new ModifiedCosine();
        modifiedCosine.score(new SimpleSpectrum(A), new SimpleSpectrum(B), 20, 30, new Deviation(10));

        for (int i=0; i < modifiedCosine.assignment.length; i+=2) {
            System.out.println("ASSIGN " + A.getMzAt(modifiedCosine.assignment[i]) + " WITH " + B.getMzAt(modifiedCosine.assignment[i+1]));
        }
        System.out.println("SCORE = " + modifiedCosine.score);
    }

    // assigns peak from left to right
    private int[] assignment;
    private double score;

    @Override
    public SpectralSimilarity score(OrderedSpectrum<Peak> left, OrderedSpectrum<Peak> right, double precursorLeft, double precursorRight, Deviation deviation) {
        return score(left,right,precursorLeft,precursorRight,deviation,1.0d);
    }

    public SpectralSimilarity score(OrderedSpectrum<Peak> left, OrderedSpectrum<Peak> right, double precursorLeft, double precursorRight, Deviation deviation, double powerIntensity) {
        if (precursorLeft <= precursorRight) {
            final DP dp = new DP(left, right, precursorLeft, precursorRight, deviation, powerIntensity);
            dp.compute();
            this.score = dp.score;
            this.assignment = dp.assignments.toArray();
        } else {
            final DP dp = new DP(right, left, precursorRight, precursorLeft, deviation, powerIntensity);
            dp.compute();
            this.score = dp.score;
            this.assignment = dp.assignments.toArray();
            for (int k=0; k < this.assignment.length; k+=2) {
                int swap = assignment[k];
                assignment[k] = assignment[k+1];
                assignment[k+1] = swap;
            }
        }
        return getSimilarity();
    }

    public double getScore() {
        return score;
    }

    public int[] getAssignment() {
        return assignment;
    }

    public SpectralSimilarity getSimilarity() {
        return new SpectralSimilarity(score, assignment.length>>1);
    }

    protected static class DP {
        OrderedSpectrum<Peak> left, right;
        double precursorLeft, precursorRight;
        Deviation dev;
        final BitSet visited;
        TIntArrayList assignments;
        double score;
        final double delta;
        final double powerIntensity;

        final short[] matches, reverseMatches, backref;
        final TDoubleArrayList dp;
        final TIntArrayList dpi;


        public DP(OrderedSpectrum<Peak> left, OrderedSpectrum<Peak> right, double precursorLeft, double precursorRight, Deviation dev, double powerIntensity) {
            this.dev = dev;
            this.left = left;
            this.right = right;
            this.precursorLeft = precursorLeft;
            this.precursorRight = precursorRight;
            this.visited = new BitSet(left.size());
            this.assignments = new TIntArrayList();
            delta = precursorRight - precursorLeft;

            this.matches = new short[left.size()];
            Arrays.fill(matches, (short) -1);
            this.reverseMatches = matches.clone();
            this.backref = new short[right.size()];
            Arrays.fill(backref, (short) -1);

            dp = new TDoubleArrayList();
            dpi = new TIntArrayList();
            this.powerIntensity = powerIntensity;
        }

        public double compute() {
            peakMatching();
            if (delta <= dev.absoluteFor(Math.max(precursorLeft,precursorRight))) {
                simpleAssignment();
            } else {
                reversePeakMatching();
                optimalAssignment();
            }
            return this.score;
        }

        // when two spectra have the same precursor, we do only match peaks directly
        private void simpleAssignment() {
            for (int i=0; i < matches.length; ++i) {
                if (matches[i]>=0) {
                    this.score += scoreFor(i, matches[i]);
                    this.assignments.add(i);
                    this.assignments.add(matches[i]);
                }
            }
        }

        // match peak from left to right
        public void peakMatching() {
            int i = 0, j = 0;
            while (i < matches.length && j < backref.length) {
                final double mzl = left.getMzAt(i), mzr = right.getMzAt(j);
                if (dev.inErrorWindow(mzl, mzr)) {
                    matches[i] = (short) j;
                    backref[j] = (short) i;
                    ++i;
                    ++j;
                } else if (mzl > mzr) {
                    ++j;
                } else {
                    ++i;
                }
            }
        }

        public void reversePeakMatching() {
            int i = reverseMatches.length - 1, j = backref.length - 1;
            while (i >= 0 && j >= 0) {
                final double mzl = precursorLeft - left.getMzAt(i), mzr = precursorRight - right.getMzAt(j);
                if (dev.inErrorWindow(mzl, mzr) && backref[j]!=i /* this is the degenerated case when mass delta is too small */) {
                    reverseMatches[i] = (short) j;
                    --i;
                    --j;
                } else if (mzl > mzr) {
                    --j;
                } else {
                    --i;
                }
            }
        }

        public void optimalAssignment() {
            score = 0d;
            for (int k = 0; k < matches.length; ++k) {
                if (visited.get(k)) continue;
                final int directMatch = matches[k];
                final int reverseMatch = reverseMatches[k];
                if (reverseMatch < 0) {
                    // easy case. We can just use direct assignment
                    if (directMatch >= 0) {
                        assignments.add(k);
                        assignments.add(directMatch);
                        score += scoreFor(k, directMatch);
                    }
                } else {
                    // we have to check for conflicting matches
                    final int conflictingMatch = backref[reverseMatch];
                    if (conflictingMatch < 0) {
                        // again, easy case. We just have to decide between two
                        // possibilities
                        double reverseScore = scoreFor(k, reverseMatch);
                        double directScore = Double.NEGATIVE_INFINITY;
                        if (directMatch >= 0) {
                            directScore = scoreFor(k, directMatch);
                        }
                        assignments.add(k);
                        if (directScore >= reverseScore) {
                            assignments.add(directMatch);
                            score += directScore;
                        } else {
                            assignments.add(reverseMatch);
                            score += reverseScore;
                        }
                    } else {
                        // difficult case. We need dynamic programing to resolve
                        // the optimal assignment
                        dp(k, directMatch, reverseMatch, conflictingMatch);
                    }
                }

            }
        }

        private void dp(int leftNode, int directMatch, int reverseMatch, int conflictingMatch) {
            dpi.clearQuick();
            dp.clearQuick();
            dp.add(scoreFor(leftNode, directMatch));
            dp.add(scoreFor(leftNode, reverseMatch));
            dp.add(0);
            dpi.add(leftNode);
            visited.set(conflictingMatch);
            int u = conflictingMatch;
            while (u >= 0) {
                final int n = dp.size();
                final double no = dp.getQuick(n - 1), rev = dp.getQuick(n - 2), match = dp.getQuick(n - 3);
                // match
                double matchScore = scoreFor(u, matches[u]);
                dp.add(Math.max(no, match) + matchScore);
                // reverse match
                double revmatchScore = scoreFor(u, reverseMatches[u]);
                final double maxall = Math.max(rev, Math.max(no, match));
                dp.add(maxall + revmatchScore);
                // no match
                dp.add(maxall);
                // next node
                dpi.add(u);
                visited.set(u);
                if (reverseMatches[u] >= 0) {
                    u = backref[reverseMatches[u]];
                } else u = -1;
            }
            backtrace();
        }

        private void backtrace() {
            int column = dp.size() - 3;
            int nodeIndex = dpi.size();
            int previousAssignment = 2; // nothing
            loop: while (true) {
                switch (previousAssignment) {
                    case 0: { // match
                        final int node = dpi.getQuick(nodeIndex);
                        assignments.add(node);
                        assignments.add(matches[node]);
                        if (column < 0) break loop;
                        // last entry can be either match or nothing
                        previousAssignment = dp.getQuick(column) > dp.getQuick(column + 2) ? 0 : 2;
                        break;
                    }
                    case 1: // reverse match
                    {
                        final int node = dpi.getQuick(nodeIndex);
                        assignments.add(node);
                        assignments.add(reverseMatches[node]);
                        if (column < 0) break loop;
                        // last entry can be either anything
                        previousAssignment = 0;
                        for (int r = 0; r < 3; ++r) {
                            if (dp.getQuick(column + r) > dp.getQuick(column + previousAssignment)) {
                                previousAssignment = r;
                            }
                        }
                        break;
                    }
                    case 2: { // nothing
                        if (column < 0) break loop;
                        // last entry can be either anything
                        previousAssignment = 0;
                        for (int r = 0; r < 3; ++r) {
                            if (dp.getQuick(column + r) > dp.getQuick(column + previousAssignment)) {
                                previousAssignment = r;
                            }
                        }
                    }

                }

                --nodeIndex;
                column -= 3;
            }
            // score is maximum of the last three dp entries
            double mxScore = 0d;
            for (int k=dp.size()-3; k < dp.size(); ++k) {
                mxScore = Math.max(dp.getQuick(k), mxScore);
            }
            score += mxScore;
        }

        private double scoreFor(int a, int b) {
            if (a < 0 || b < 0) return Double.NEGATIVE_INFINITY;
            return Math.pow(left.getIntensityAt(a) * right.getIntensityAt(b), powerIntensity);
        }
    }
}
