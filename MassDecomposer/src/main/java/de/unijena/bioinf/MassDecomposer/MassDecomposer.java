package de.unijena.bioinf.MassDecomposer;

import java.util.*;

/**
 * Decomposes a given mass over an alphabet, returning all decompositions which mass equals the given mass
 * considering a given deviation.
 *
 * @see @see <a href="http://bio.informatik.uni-jena.de/bib2html/downloads/2005/BoeckerLiptak_MoneyChangingProblem_COCOON_2005.pdf">The Money Changing Problem revisited by Sebastian Böcker and Zsuzsanna Lipták</a>
 * @param <T> type of the alphabet's characters
 */
public class MassDecomposer<T> {

    protected long[][] ERT;
    protected double precision;
    protected final double errorPpm;
    protected final double absMassError;
    protected final List<Weight<T>> weights;
    protected double minError, maxError;
    protected final Alphabet<T> alphabet;
    protected final int[] orderedCharacterIds;
    protected DecompositionValidator<T> validator;

    /**
     * @param precision mass precision. A precision of 1e-3 means that three positions after decimal point are
     *                  considered for input masses
     * @param errorPPM relative mass error in parts per million. 10 ppm means that we have to consider a deviation
     *                 of 0.001% of the input mass
     * @param absMassError absolute mass error. The same as errorPPM but with absolute value. The considered
     *                     deviation is always the maximum of the relative and absolute error
     * @param alphabet the alphabet the mass is decomposed over
     */
    public MassDecomposer(double precision, double errorPPM, double absMassError, Alphabet<T> alphabet) {
        this.errorPpm = errorPPM;
        this.absMassError = absMassError;
        this.precision = precision;
        final int n = alphabet.size();
        this.weights = new ArrayList<Weight<T>>(n);
        for (int i=0; i < n; ++i) {
            weights.add(new Weight<T>(alphabet.get(i), alphabet.weightOf(i)));
        }
        Collections.sort(weights);
        this.alphabet = alphabet;
        this.orderedCharacterIds = new int[alphabet.size()];
        for (int i=0; i < alphabet.size(); ++i) {
            orderedCharacterIds[i] = alphabet.indexOf(weights.get(i).getOwner());
        }
    }

    public double getErrorPpm() {
        return errorPpm;
    }

    public double getAbsMassError() {
        return absMassError;
    }

    public Alphabet<T> getAlphabet() {
        return alphabet;
    }

    /**
     * @return the active validator for the decomposer
     */
    public DecompositionValidator<T> getValidator() {
        return validator;
    }

    /**
     * set a validator which removes invalid decompositions as early as possible before returning them
     * in {@link #decompose(double)}}
     * @param validator
     */
    public void setValidator(DecompositionValidator<T> validator) {
        this.validator = validator;
    }

    /**
     * Check if a mass is decomposable. This is done in constant time (especially: it is very very very fast!).
     * But it doesn't check if there is a valid decomposition. Therefore, even if the method returns true,
     * all decompositions may be invalid for the given validator or given bounds.
     * #decompose(mass) uses this function before starting the decomposition, therefore this method should only
     * be used if you don't want to start the decomposition algorithm.
     * @param mass
     * @return
     */
    public boolean maybeDecomposable(double mass) {
        init();
        final Interval range = integerBound(mass, getErrorForMass(mass));
        final long a = weights.get(0).getIntegerMass();
        for (long i = range.getMin(); i <= range.getMax(); ++i) {
            final int r = toInt(i % a);
            if (i >= ERT[r][weights.size()-1]) return true;
        }
        return false;
    }

    /**
     * returns the indizes of the characters in the same order as the compomeres
     * Don't modify this array, because the decomposer use it internally.
     */
    public int[] getCharacterIndizes() {
        return orderedCharacterIds;
    }

    /**
     * @see #decompose(double, java.util.Map)
     */
    public List<int[]> decompose(double mass) {
        return decompose(mass, null);
    }

    /**
     * computes all decompositions for the given mass. The runtime depends only on the number of characters and the
     * number of decompositions. Therefore this method is very fast as long as the number of decompositions is low.
     * Unfortunately, the number of decompositions increases nearly exponential in the number of characters and in the
     * input mass.
     *
     * This function can be called in multiple threads in parallel, because it does not modify the decomposer
     * @param mass input mass which should be decomposed
     * @param boundaries a map which maps some characters of the alphabet to an interval. The decomposer will discard
     *                   all decompositions for which the amound of the characters is not in the given interval. This is
     *                   done as early as possible, resulting a better performance than checking the boundary during the
     *                   validation step
     * @return list of decompositions. Use the {@link #getCharacterIndizes} and {@link Alphabet#get(int)} method to map the
     *         indizes of compomere to the characters.
     */
    public List<int[]> decompose(final double mass, Map<T, Interval> boundaries){
        init();
        if (mass == 0d) return Collections.emptyList();
        if (mass < 0d) throw new IllegalArgumentException("Expect positive mass for decomposition");
        double calcMass = mass;
        final double absError = getErrorForMass(mass);
        final int[] minValues = new int[weights.size()];
        final int[] boundsarray = new int[weights.size()];
        boolean minAllZero = true;
        Arrays.fill(boundsarray, Integer.MAX_VALUE);
        if (boundaries!=null && !boundaries.isEmpty()) {
            for (int i = 0; i < boundsarray.length; i++) {
                T el = weights.get(i).getOwner();
                Interval range = boundaries.get(el);
                if (range != null) {
                    boundsarray[i] = toInt(range.getMax()-range.getMin());
                    minValues[i] = toInt(range.getMin());
                    if (minValues[i] > 0) {
                        minAllZero = false;
                        calcMass -= weights.get(i).getMass() * range.getMin();
                    }
                }
            }
        }
        final ArrayList<int[]> results = new ArrayList<int[]>();
        ArrayList<int[]> rawDecompositions = null;
        final Interval interval = integerBound(calcMass, absError);
        if (!minAllZero && (Math.abs(calcMass) <= absError)) results.add(minValues);
        for (long m = interval.getMin(); m < interval.getMax(); ++m) {
            rawDecompositions = integerDecompose(m, boundsarray);
            for (int i=0; i < rawDecompositions.size(); ++i) {
                final int[] decomp = rawDecompositions.get(i);
                if (!minAllZero) {
                    for (int j=0; j < minValues.length; ++j) {
                        decomp[j] += minValues[j];
                    }
                }
                if (Math.abs(calcMass(decomp) - mass) > absError) continue;
                if ((validator == null) || validator.validate(decomp, orderedCharacterIds, alphabet)) {
                    results.add(decomp);
                }
            }
        }
        return results;
    }

    protected ArrayList<int[]> integerDecompose(long mass, int[] bounds){
        // Find compomers
        ArrayList<int[]> result = new ArrayList<int[]>();
        int k = weights.size();
        int[] c = new int[k], deepCopy;
        long[] j = new long[k], m = new long[k], lbound = new long[k], r = new long[k];
        boolean flagWhile = false; // flag wether we are in the while-loop or not
        final long a = weights.get(0).getIntegerMass();
        // Init
        for (int i=1; i<k; ++i){
            lbound[i] = Long.MAX_VALUE; // this is just to ensure, that lbound < m in the first iteration
        }

        int i = k-1;
        m[i] = mass; // m[i] corresponds to M, m[i-1] ^= m
        while (i != k){
            if (i == 0){
                deepCopy = new int[weights.size()];
                for (int index=0; index<c.length; ++index) deepCopy[index] = c[index];
                deepCopy[0] = (int) (m[i]/a);
                if (deepCopy[0] <= bounds[0]) result.add(deepCopy);
                ++i; // "return" from recursion
                flagWhile = true; // in this recursion-depth we are in the while-loop, cause the next recursion (the one we just exited) was called
                m[i-1] -= weights.get(i).getLcm(); // execute the rest of the while
                c[i] += weights.get(i).getL();
            } else {
                if (flagWhile){
                    //de.fsu.ms2tool.application.io.output.AlgorithmLogging.get.info("lbound: " +lbound[i]+" m: "+m[i-1]);
                    if (m[i-1] >= lbound[i] && c[i] <= bounds[i]){ //currently in while loop
                        //de.fsu.ms2tool.application.io.output.AlgorithmLogging.get.info("i: "+(i-1)+" m: "+m[i-1]);
                        --i; // "do" recursive call
                    } else {
                        flagWhile = false; //
                    }
                } else { //we are in the for-loop
                    if (j[i] < weights.get(i).getL() && m[i]-j[i]*weights.get(i).getIntegerMass()>=0){
                        c[i] = (int) j[i];
                        m[i-1] = m[i]-j[i]*weights.get(i).getIntegerMass();
                        r[i] = m[i-1]%a;
                        lbound[i] = ERT[(int)r[i]][i-1];
                        flagWhile = true; // call the while loop
                        ++j[i];
                    } else { //exit for loop
                        // reset "function variables"
                        lbound[i] = Long.MAX_VALUE;
                        j[i] = 0;
                        c[i] = 0;
                        ++i; // "return" from recursion
                        if (i != k) { // only if we are not done
                            flagWhile = true; // in this recursion-depth we are in the while-loop, cause the next recursion was called
                            m[i-1] -= weights.get(i).getLcm(); // execute the rest of the while
                            c[i] += weights.get(i).getL();
                        }
                    }
                }
            } // end if i == 0
        } // end while
        return result;
    } // end function

    /**
     * Initializes the decomposer. Computes the extended residue table. This have to be done only one time for
     * a given alphabet, independently from the masses you want to decompose. This method is called automatically
     * if you compute the decompositions, so call it only if you want to control the time of the initialisation.
     */
    public void init() {
        if (ERT != null) return;
        discretizeMasses();
        divideByGCD();
        computeLCMs();
        calcERT();
        computeErrors();
    }

    protected double calcMass(int[] input){
        double result = 0.0;
        for (int i = 0; i < input.length; ++i){
            result += input[i]*weights.get(i).getMass();
        }
        return result;
    }

    protected void calcERT(){
        long firstLongVal = weights.get(0).getIntegerMass();
        ERT = new long[(int)firstLongVal][weights.size()];
        long d, r, n, argmin;

        //Init
        ERT[0][0] = 0;
        for (int i = 1; i < ERT.length; ++i){
            ERT[i][0] = Long.MAX_VALUE; // should be infinity
        }

        //Filling the Table, j loops over columns
        for (int j = 1; j < ERT[0].length; ++j){
            ERT[0][j] = 0; // Init again
            d = gcd(firstLongVal, weights.get(j).getIntegerMass());
            for (long p = 0; p < d; p++){ // Need to start d Round Robin loops
                if (p == 0) {
                    n = 0; // 0 is the min in the complete RT or the first p-loop
                } else {
                    n = Long.MAX_VALUE; // should be infinity
                    argmin = p;
                    for (long i = p; i<ERT.length; i += d){ // Find Minimum in specific part of ERT
                        if (ERT[(int)i][j-1] < n){
                            n = ERT[(int)i][j-1];
                            argmin = i;
                        }
                    }
                    ERT[(int)argmin][j]= n;
                }
                if (n == Long.MAX_VALUE){ // Minimum of the specific part of ERT was infinity
                    for (long i = p; i<ERT.length; i += d){ // Fill specific part of ERT with infinity
                        ERT[(int)i][j] = Long.MAX_VALUE;
                    }
                } else { // Do normal loop
                    for (long i = 1; i < ERT.length/d; ++i){ // i is just a counter
                        n += weights.get(j).getIntegerMass();
                        r = n % firstLongVal;
                        if (ERT[(int)r][j-1] < n) n = ERT[(int)r][j-1]; // get the min
                        ERT[(int)r][j] = n;
                    }
                }
            } // end for p
        } // end for j
    }

    protected void discretizeMasses() {
        // compute integer masses
        for (int i=0; i  < weights.size(); ++i) {
            final Weight<T> weight = weights.get(i);
            weight.setIntegerMass((long)(weight.getMass() / precision));
        }
    }

    protected void divideByGCD() {
        if (weights.size() > 0) {
            long d = gcd(weights.get(0).getIntegerMass(), weights.get(1).getIntegerMass());
            for (int i=2; i < weights.size(); ++i) {
                d = gcd(d, weights.get(i).getIntegerMass());
                if (d == 1) return;
            }
            precision *= d;
            for (Weight<T> weight : weights) {
                weight.setIntegerMass(weight.getIntegerMass() / d);
            }
        }
    }

    protected void computeLCMs() {
        final Weight<T> first = weights.get(0);
        first.setL(1);
        first.setLcm(first.getIntegerMass());

        for(int i=1; i<weights.size();i++){
            final Weight<T> weight = weights.get(i);
            long temp = first.getIntegerMass() / gcd(first.getIntegerMass(), weight.getIntegerMass());
            weight.setL(temp);
            weight.setLcm(temp * weight.getIntegerMass());
        }
    }

    protected static long gcd(long u, long v) {
        long r = 0;

        while (v != 0) {
            r = u % v;
            u = v;
            v = r;
        }
        return u;
    }

    protected void computeErrors() {
        this.minError = 0d;
        this.maxError = 0d;
        for (Weight<T> weight : weights) {
            final double error = (precision * weight.getIntegerMass() - weight.getMass()) / weight.getMass();
            minError = Math.min(minError, error);
            maxError = Math.max(maxError, error);
        }
    }

    protected Interval integerBound(double mass, double error) {
        final double absError = error;
        return new Interval(
                Math.max(0, (long)Math.ceil((1 + minError) * (mass - absError) / precision)),
                Math.max(0, (long)Math.floor((1 + maxError) * (mass + absError) / precision))
        );
    }

    protected double getErrorForMass(double mass){
        return Math.max(errorPpm*1e-6*mass, absMassError);
    }

    protected static int toInt(long value) {
        if (value > Integer.MAX_VALUE) throw new ArithmeticException("Can't cast " + value + " to Integer");
        return (int)value;
    }

}
