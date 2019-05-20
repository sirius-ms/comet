package de.unijena.bioinf.ChemistryBase.ms.ft;

import de.unijena.bioinf.ms.annotations.TreeAnnotation;

/**
 * This annotation is used when a tree is "beautiful" ;)
 * either it explains enough peaks or we already maxed out its tree size score
 *
 * The beautificationScoreAddition is the additional score a tree gains due to the beautification. It is
 * canceled in the scoring, such that the original score (without beautification) is preserved
 */
public final class Beautified implements TreeAnnotation  {

    public static enum State {
        UGLY, IN_PROCESS, BEAUTIFIED
    }

    public double getNodeBoost() {
        return nodeBoost;
    }

    public double getBeautificationPenalty() {
        return beautificationPenalty;
    }

    public final static Beautified IS_UGGLY = new Beautified(State.UGLY,0d,0d);

    protected final State state;
    protected final double nodeBoost;
    // difference between score of beautified instance and original instance
    protected final double beautificationPenalty;

    public static final String PENALTY_KEY = "BeautificationPenalty";

    private Beautified(State state, double nodeBoost, double beautificationPenalty) {
        this.state = state;
        this.nodeBoost = nodeBoost;
        this.beautificationPenalty = beautificationPenalty;
    }

    public static Beautified ugly() {
        return IS_UGGLY;
    }

    public static Beautified beautified(double nodeBoost, double beautificationPenalty) {
        return new Beautified(State.BEAUTIFIED, nodeBoost, beautificationPenalty);
    }
    public static Beautified inProcess(double nodeBoost) {
        return new Beautified(State.IN_PROCESS, nodeBoost, 0d);
    }

    public boolean isBeautiful() {
        return state==State.BEAUTIFIED;
    }

    public boolean isInProcess() {
        return state == State.IN_PROCESS;
    }

}
