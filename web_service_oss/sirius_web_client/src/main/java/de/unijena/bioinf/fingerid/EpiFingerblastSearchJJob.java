package de.unijena.bioinf.fingerid;

import de.unijena.bioinf.ChemistryBase.algorithm.scoring.SScored;
import de.unijena.bioinf.ChemistryBase.algorithm.scoring.Scored;
import de.unijena.bioinf.ChemistryBase.chem.MolecularFormula;
import de.unijena.bioinf.ChemistryBase.fp.ProbabilityFingerprint;
import de.unijena.bioinf.ChemistryBase.ms.ft.FTree;
import de.unijena.bioinf.ChemistryBase.ms.ft.Fragment;
import de.unijena.bioinf.chemdb.FingerprintCandidate;
import de.unijena.bioinf.fingerid.blast.BayesnetScoring;
import de.unijena.bioinf.fingerid.blast.FingerblastResult;
import de.unijena.bioinf.fragmenter.*;
import de.unijena.bioinf.jjobs.BasicJJob;
import de.unijena.bioinf.jjobs.JJob;
import de.unijena.bioinf.jjobs.Partition;
import de.unijena.bioinf.ms.properties.PropertyManager;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;

public class EpiFingerblastSearchJJob extends FingerblastSearchJJob {

    private final boolean rankWithEpimetheus;

    public EpiFingerblastSearchJJob(@NotNull CSIPredictor predictor, @Nullable BayesnetScoring bayesnetScoring, FingerIdResult fingerIdResult, boolean rankWithEpimetheus){
        this(predictor, bayesnetScoring, fingerIdResult.getSourceTree(), fingerIdResult.getMolecularFormula(), fingerIdResult.getPredictedFingerprint(), rankWithEpimetheus);
    }

    public EpiFingerblastSearchJJob(@NotNull CSIPredictor predictor, FingerIdResult fingerIdResult, boolean rankWithEpimetheus){
        this(predictor, null, fingerIdResult, rankWithEpimetheus);
    }

    public EpiFingerblastSearchJJob(@NotNull CSIPredictor predictor, FTree tree, ProbabilityFingerprint fp, MolecularFormula formula, boolean rankWithEpimetheus){
        this(predictor, null, tree, formula, fp, rankWithEpimetheus);
    }

    public EpiFingerblastSearchJJob(@NotNull CSIPredictor predictor, @Nullable BayesnetScoring bayesnetScoring, FTree fTree, MolecularFormula formula, ProbabilityFingerprint fp, boolean rankWithEpimetheus){
        super(predictor, bayesnetScoring);
        this.setFtree(fTree);
        this.setFormula(formula);
        this.setFingerprint(fp);
        this.rankWithEpimetheus = rankWithEpimetheus;
    }

    @Override
    protected FingerblastResult compute() throws Exception {
        final FingerblastResult result = super.compute();
        this.checkForInterruption();
        if(!this.rankWithEpimetheus) return result;
        this.checkForInterruption();

        final Collection<FingerprintCandidate> candidates = result.getResults().stream().map(SScored::getCandidate).toList();
        this.checkForInterruption();
        final List<JJob<List<Scored<FingerprintCandidate>>>> scoringJobs = this.makeScoringJobs(candidates);
        this.checkForInterruption();
        scoringJobs.forEach(this::submitSubJob);
        this.checkForInterruption();
        final List<Scored<FingerprintCandidate>> epiScoredCandidates = scoringJobs.stream().flatMap(j -> j.takeResult().stream()).sorted(Comparator.reverseOrder()).toList();
        this.checkForInterruption();
        return new FingerblastResult(epiScoredCandidates);
    }

    private List<JJob<List<Scored<FingerprintCandidate>>>> makeScoringJobs(Collection<FingerprintCandidate> candidates){
        final List<List<FingerprintCandidate>> batches = Partition.ofNumber(candidates, PropertyManager.getNumberOfThreads());
        final SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        final FTree fTree = this.ftree;
        return batches.stream().map(batch -> new BasicJJob<List<Scored<FingerprintCandidate>>>(JobType.CPU){
            @Override
            protected List<Scored<FingerprintCandidate>> compute(){
                final ArrayList<Scored<FingerprintCandidate>> results = new ArrayList<>(batch.size());
                for(final FingerprintCandidate candidate : batch){
                    try {
                        final double score = scoreCandidate(candidate, fTree, smilesParser);
                        final Scored<FingerprintCandidate> scoredCandidate = new Scored<>(candidate, score);
                        results.add(scoredCandidate);
                    }catch(InvalidSmilesException e){
                        e.printStackTrace();
                        EpiFingerblastSearchJJob.this.logError("Invalid SMILES string of candidate. Corresponding candidate gets score Double.NEGATIVE_INFINITY.");
                        results.add(new Scored<>(candidate, Double.NEGATIVE_INFINITY));
                    }
                }
                results.sort(Comparator.reverseOrder());
                return results;
            }
        }).collect(Collectors.toList());
    }

    /*
    // #####################################
    // #### This is the original method.####
    // #####################################
    private static double scoreCandidate(FingerprintCandidate candidate, FTree fTree, SmilesParser smiParser) throws InvalidSmilesException {
        // INITIALISATION:
        final MolecularGraph molecule = new MolecularGraph(fTree.getRoot().getFormula(), smiParser.parseSmiles(candidate.getSmiles()));
        final EMFragmenterScoring2 scoring = new EMFragmenterScoring2(molecule, fTree);
        final CriticalPathSubtreeCalculator subtreeCalc = new CriticalPathSubtreeCalculator(fTree, molecule, scoring, true);
        subtreeCalc.setMaxNumberOfNodes(50000);

        // 'subtreeCalc' has to be initialized. This encompasses:
        // - the combinatorial fragmentation of 'molecule' to construct a CombinatorialGraph object with added terminal nodes (peaks)
        // - the initialisation procedures to execute the CriticalPath heuristic
        final HashSet<MolecularFormula> fset = new HashSet<>();
        for (final Fragment ft : fTree.getFragmentsWithoutRoot()) {
            fset.add(ft.getFormula());
            fset.add(ft.getFormula().add(MolecularFormula.getHydrogen()));
            fset.add(ft.getFormula().add(MolecularFormula.getHydrogen().multiply(2)));
            if (ft.getFormula().numberOfHydrogens() > 0)
                fset.add(ft.getFormula().subtract(MolecularFormula.getHydrogen()));
            if (ft.getFormula().numberOfHydrogens() > 1)
                fset.add(ft.getFormula().subtract(MolecularFormula.getHydrogen().multiply(2)));
        }
        subtreeCalc.initialize((node, nnodes, nedges) -> {
            if (fset.contains(node.getFragment().getFormula())) return true;
            return (node.getTotalScore() > -5f);
        });

        // Compute the CombinatorialSubtree using the CriticalPath heuristic and return its score:
        subtreeCalc.computeSubtree();
        return subtreeCalc.getScore();
    }
     */

    private static double scoreCandidate(FingerprintCandidate candidate, FTree fTree, SmilesParser smiParser) throws InvalidSmilesException {
        // INITIALISATION:
        final MolecularGraph molecule = new MolecularGraph(fTree.getRoot().getFormula(), smiParser.parseSmiles(candidate.getSmiles()));
        final MyScoring scoring = new MyScoring();
        final CriticalPathSubtreeCalculator subtreeCalc = new CriticalPathSubtreeCalculator(fTree, molecule, scoring, true);
        subtreeCalc.setMaxNumberOfNodes(50000);

        // 'subtreeCalc' has to be initialized. This encompasses:
        // - the combinatorial fragmentation of 'molecule' to construct a CombinatorialGraph object with added terminal nodes (peaks)
        // - the initialisation procedures to execute the CriticalPath heuristic
        final HashSet<MolecularFormula> fset = new HashSet<>();
        for (final Fragment ft : fTree.getFragmentsWithoutRoot()) {
            fset.add(ft.getFormula());
            fset.add(ft.getFormula().add(MolecularFormula.getHydrogen()));
            fset.add(ft.getFormula().add(MolecularFormula.getHydrogen().multiply(2)));
            if (ft.getFormula().numberOfHydrogens() > 0)
                fset.add(ft.getFormula().subtract(MolecularFormula.getHydrogen()));
            if (ft.getFormula().numberOfHydrogens() > 1)
                fset.add(ft.getFormula().subtract(MolecularFormula.getHydrogen().multiply(2)));
        }
        subtreeCalc.initialize((node, nnodes, nedges) -> {
            if (fset.contains(node.getFragment().getFormula())) return true;
            return (node.getDepth() < 3);
        });

        // Compute the CombinatorialSubtree using the CriticalPath heuristic and return its score:
        subtreeCalc.computeSubtree();
        return subtreeCalc.getScore();
    }

    @Override
    public String identifier(){
        String filePath = this.ftree == null ? "Null" : this.ftree.getAnnotation(de.unijena.bioinf.ChemistryBase.data.DataSource.class).orElse(new de.unijena.bioinf.ChemistryBase.data.DataSource(new File(""))).toString();
        return super.identifier() + " | " + filePath + " | WithEpimetheus";
    }

    private static class MyScoring implements CombinatorialFragmenterScoring{

        private final static double p = 0.2;
        private final static double q = 0.1;

        @Override
        public double scoreBond(IBond bond, boolean direction) {
            return 0;
        }

        @Override
        public double scoreFragment(CombinatorialNode fragment) {
            return 0;
        }

        @Override
        public double scoreEdge(CombinatorialEdge edge) {
            final CombinatorialNode sourceNode = edge.getSource();
            final CombinatorialNode targetNode = edge.getTarget();
            if(targetNode.getFragment().isInnerNode()) return 0d;

            // traverse to the root:
            double fragmentationProb = 1;
            CombinatorialNode currentNode = sourceNode;
            while(!currentNode.getIncomingEdges().isEmpty()){
                // find parent with shortest depth:
                for(CombinatorialEdge inEdge : currentNode.getIncomingEdges()){
                    final CombinatorialNode parent = inEdge.getSource();
                    if(parent.getDepth()+1 == currentNode.getDepth()){
                        fragmentationProb *= inEdge.getCut2() == null ? p : q;
                        currentNode = parent;
                        break;
                    }
                }
            }

            return fragmentationProb * targetNode.getFragment().getPeakIntensity();
        }

    }
}
