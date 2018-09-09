/*
 * AntibioticModel.java
 *
 * Copyright © 2017-2018
 * University of Stirling
 * All Rights Reserved
 */

package uk.ac.stir.cs.antibiotic;

import java.util.*;
import java.util.concurrent.*;

/**
 * A stochastic model for antibiotic treatment. Candidate solutions may be
 * evaluated by calling the {@link #evaluate(int...)} method with an array of
 * integers. This implementation is consistent with the MATLAB reference
 * implementation. Use {@link uk.ac.stir.cs.antibiotic.AntibioticProblem}
 * wrapper for problem with specified objective or objectives.
 *
 * @author Lee A. Christie
 * @author Iona K. Paterson
 * @author Andy Hoyle
 *
 * @version 4.0.6
 */
public final class AntibioticModel {

    // Valued used to match MATLAB reference implementation
    static final double MATLAB_REF_W1 = 1.0;
    static final double MATLAB_REF_W2 = 0.1;
    static final int MATLAB_REF_V_MAX = 184;
    static final int MATLAB_REF_MAX_CONC = 60;
    static final double MATLAB_REF_PENALTY_V_MAX = Math.pow(10.0, 10.0);
    static final double MATLAB_REF_PENALTY_CONC = Math.pow(10.0, 10.0);

    private static final int DEFAULT_INITIAL_LOAD_1 = 900;
    private static final int DEFAULT_INITIAL_LOAD_2 = 100;

    // model settings
    private final int runs;
    private final int initialLoad1;
    private final int initialLoad2;
    private final Random random;
    private final int targetFailures;
    private final int maximumRuns;

    // Biological parameters
    private static final double r = 2.7726; // Reproduction rate of Susceptible
    private static final double c1 = 0.2; // Cost of carrying resistance
    private static final int K = 1000; // Carrying Capacity
    private static final double ms = 0.2; // Mortality Rate of Susceptibles
    private static final double mr = 0.2; // Mortality Rate of Resistants
    private static final double MaxS = r - ms; // Max net growth rate in absence
                                               // of AB for susceptibles
    private static final double MinS = -2.1; // Min net growth rate at high
                                             // AB concentrations for
                                             //susceptibles
    private static final int MICS = 16; // Pharmocodynamic MIC for Suceptible
    private static final int kS = 4; // Hill Coefficient
    private final double MaxR1 = r * (1 - c1) - mr; // Max net growth rate in
                                                    // absence of AB for
                                                    // resistants
    private static final double MinR1 = -2.1; // Min net growth rate at high AB
                                              // concentrations for resistants
    private static final int MICR1 = 32; // Pharmocodynamic MIC for Resistants
    private static final int kR1 = 4; // Hill Coefficient
    private static final double a = 0.48; // Degradation rate of AB

    /**
     * Creates an instance of the antibiotic model which runs for the specified
     * number of iterations (samples). The default thread-local random number
     * generator is used. The default initial loads of 900 for strain 1 and 100
     * for strain 2 will be used.
     *
     * @param runs the number of iterations (samples) to run, at least 1.
     * @return the antibiotic model.
     */
    public static AntibioticModel fixedSampleSize(
            final int runs) {
        if (runs < 1) {
            throw new IllegalArgumentException(
                    "runs = " + runs + ", expected: >= 1.");
        }
        return new AntibioticModel(runs,
                                   DEFAULT_INITIAL_LOAD_1,
                                   DEFAULT_INITIAL_LOAD_2,
                                   null,
                                   -1,
                                   -1);
    }

    /**
     * Creates an instance of the antibiotic model which runs for the specified
     * number of iterations (samples) using the given random number generator.
     * The default initial loads of 900 for strain 1 and 100 for strain 2 will
     * be used.
     *
     * @param runs the number of iterations (samples) to run, at least 1.
     * @param random the random number generator, not null.
     * @return the antibiotic model.
     */
    public static AntibioticModel fixedSampleSize(
            final int runs,
            final Random random) {
        if (runs < 1) {
            throw new IllegalArgumentException(
                    "runs = " + runs + ", expected: >= 1.");
        }
        if (random == null) {
            throw new NullPointerException("random");
        }
        return new AntibioticModel(runs,
                                   DEFAULT_INITIAL_LOAD_1,
                                   DEFAULT_INITIAL_LOAD_2,
                                   random,
                                   -1,
                                   -1);
    }

    /**
     * Creates an instance of the antibiotic model which runs for the specified
     * number of iterations (samples). The default thread-local random number
     * generator is used. Total bacterial load must be between 0 and 1000,
     * inclusive.
     *
     * @param runs the number of iterations (samples) to run, at least 1.
     * @param initialLoad1 the initial bacterial load for strain 1. Between 0
     * and 1000, inclusive.
     * @param initialLoad2 the initial bacterial load for strain 2. Between 0
     * and 1000, inclusive.
     * @return the antibiotic model.
     */
    public static AntibioticModel fixedSampleSize(
            final int runs,
            final int initialLoad1,
            final int initialLoad2) {
        if (runs < 1) {
            throw new IllegalArgumentException(
                    "runs = " + runs + ", expected: >= 1.");
        }
        if (initialLoad1 < 0 || initialLoad1 > 1000) {
            throw new IllegalArgumentException(
                    "initialLoad1 = " + initialLoad1 + ", expected: [0, 1000]");
        }
        if (initialLoad2 < 0 || initialLoad2 > 1000) {
            throw new IllegalArgumentException(
                    "initialLoad2 = " + initialLoad2 + ", expected: [0, 1000]");
        }
        final int totalLoad = initialLoad1 + initialLoad2;
        if (totalLoad < 0 || totalLoad > 1000) {
            throw new IllegalArgumentException(
                    "initialLoad1 + initialLoad2 = "
                    + (initialLoad1 + initialLoad2) + ", expected: [0, 1000]");
        }
        return new AntibioticModel(runs,
                                   initialLoad1,
                                   initialLoad2,
                                   null,
                                   -1,
                                   -1);
    }

    /**
     * Creates an instance of the antibiotic model which runs for the specified
     * number of iterations (samples) using the given random number generator.
     * Total bacterial load must be between 0 and 1000, inclusive.
     *
     * @param runs the number of iterations (samples) to run, at least 1.
     * @param initialLoad1 the initial bacterial load for strain 1. Between 0
     * and 1000, inclusive.
     * @param initialLoad2 the initial bacterial load for strain 2. Between 0
     * and 1000, inclusive.
     * @param random the random number generator, not null.
     * @return the antibiotic model.
     */
    public static AntibioticModel fixedSampleSize(
            final int runs,
            final int initialLoad1,
            final int initialLoad2,
            final Random random) {
        if (runs < 1) {
            throw new IllegalArgumentException(
                    "runs = " + runs + ", expected: >= 1.");
        }
        if (initialLoad1 < 0 || initialLoad1 > 1000) {
            throw new IllegalArgumentException(
                    "initialLoad1 = " + initialLoad1 + ", expected: [0, 1000]");
        }
        if (initialLoad2 < 0 || initialLoad2 > 1000) {
            throw new IllegalArgumentException(
                    "initialLoad2 = " + initialLoad2 + ", expected: [0, 1000]");
        }
        final int totalLoad = initialLoad1 + initialLoad2;
        if (totalLoad < 0 || totalLoad > 1000) {
            throw new IllegalArgumentException(
                    "initialLoad1 + initialLoad2 = "
                    + (initialLoad1 + initialLoad2) + ", expected: [0, 1000]");
        }
        if (random == null) {
            throw new NullPointerException("random");
        }
        return new AntibioticModel(runs,
                                   initialLoad1,
                                   initialLoad2,
                                   random,
                                   -1,
                                   -1);
    }

    /**
     * Creates an instance of the antibiotic model which runs either the maximum
     * number of iterations (samples) specified, or until the specified number
     * of failures has been reached, whichever happens first. The default
     * thread-local random number generator is used. The default initial loads
     * of 900 for strain 1 and 100 for strain 2 will be used.
     *
     * @param targetFailures the number of failures to reach before stopping the
     *                       simulation
     * @param maximumRuns the maximum number of times to run the model, at least
     *                    1, must be greater or equal to the targetFailures.
     * @return the antibiotic model.
     */
    public static AntibioticModel dynamicSampleSize(
            final int targetFailures,
            final int maximumRuns) {
        if (targetFailures < 1) {
            throw new IllegalArgumentException(
                    "targetFailures = " + targetFailures + ", expected: >= 1.");
        }
        if (maximumRuns < 1) {
            throw new IllegalArgumentException(
                    "maximumRuns = " + maximumRuns + ", expected: >= 1.");
        }
        if (targetFailures > maximumRuns) {
            throw new IllegalArgumentException(
                    "targetFailures = " + targetFailures
                  + ", maximumRuns = " + maximumRuns
                  + ", expected: targetFailures <= maximumRuns.");
        }
        return new AntibioticModel(-1,
                                   DEFAULT_INITIAL_LOAD_1,
                                   DEFAULT_INITIAL_LOAD_2,
                                   null,
                                   targetFailures,
                                   maximumRuns);
    }

    /**
     * Creates an instance of the antibiotic model which runs either the maximum
     * number of iterations (samples) specified, or until the specified number
     * of failures has been reached, whichever happens first using the given
     * random number generator. The default initial loads of 900 for strain 1
     * and 100 for strain 2 will be used.
     *
     * @param targetFailures the number of failures to reach before stopping the
     *                       simulation
     * @param maximumRuns the maximum number of times to run the model, at least
     *                    1, must be greater or equal to the targetFailures.
     * @param random the random number generator, not null.
     * @return the antibiotic model.
     */
    public static AntibioticModel dynamicSampleSize(
            final int targetFailures,
            final int maximumRuns,
            final Random random) {
        if (targetFailures < 1) {
            throw new IllegalArgumentException(
                    "targetFailures = " + targetFailures + ", expected: >= 1.");
        }
        if (maximumRuns < 1) {
            throw new IllegalArgumentException(
                    "maximumRuns = " + maximumRuns + ", expected: >= 1.");
        }
        if (targetFailures > maximumRuns) {
            throw new IllegalArgumentException(
                    "targetFailures = " + targetFailures
                  + ", maximumRuns = " + maximumRuns
                  + ", expected: targetFailures <= maximumRuns.");
        }
        if (random == null) {
            throw new NullPointerException("random");
        }
        return new AntibioticModel(-1,
                                   DEFAULT_INITIAL_LOAD_1,
                                   DEFAULT_INITIAL_LOAD_2,
                                   random,
                                   targetFailures,
                                   maximumRuns);
    }


    /**
     * Creates an instance of the antibiotic model which runs either the maximum
     * number of iterations (samples) specified, or until the specified number
     * of failures has been reached, whichever happens first. The default
     * thread-local random number generator is used. Total bacterial load must
     * be between 0 and 1000, inclusive.
     *
     * @param targetFailures the number of failures to reach before stopping the
     *                       simulation
     * @param maximumRuns the maximum number of times to run the model, at least
     *                    1, must be greater or equal to the targetFailures.
     * @param initialLoad1 the initial bacterial load for strain 1. Between 0
     * and 1000, inclusive.
     * @param initialLoad2 the initial bacterial load for strain 2. Between 0
     * and 1000, inclusive.
     * @return the antibiotic model.
     */
    public static AntibioticModel dynamicSampleSize(
            final int targetFailures,
            final int maximumRuns,
            final int initialLoad1,
            final int initialLoad2) {
        if (targetFailures < 1) {
            throw new IllegalArgumentException(
                    "targetFailures = " + targetFailures + ", expected: >= 1.");
        }
        if (maximumRuns < 1) {
            throw new IllegalArgumentException(
                    "maximumRuns = " + maximumRuns + ", expected: >= 1.");
        }
        if (targetFailures > maximumRuns) {
            throw new IllegalArgumentException(
                    "targetFailures = " + targetFailures
                  + ", maximumRuns = " + maximumRuns
                  + ", expected: targetFailures <= maximumRuns.");
        }
        if (initialLoad1 < 0 || initialLoad1 > 1000) {
            throw new IllegalArgumentException(
                    "initialLoad1 = " + initialLoad1 + ", expected: [0, 1000]");
        }
        if (initialLoad2 < 0 || initialLoad2 > 1000) {
            throw new IllegalArgumentException(
                    "initialLoad2 = " + initialLoad2 + ", expected: [0, 1000]");
        }
        final int totalLoad = initialLoad1 + initialLoad2;
        if (totalLoad < 0 || totalLoad > 1000) {
            throw new IllegalArgumentException(
                    "initialLoad1 + initialLoad2 = "
                    + (initialLoad1 + initialLoad2) + ", expected: [0, 1000]");
        }
        return new AntibioticModel(-1,
                                   initialLoad1,
                                   initialLoad2,
                                   null,
                                   targetFailures,
                                   maximumRuns);
    }

    /**
     * Creates an instance of the antibiotic model which runs either the maximum
     * number of iterations (samples) specified, or until the specified number
     * of failures has been reached, whichever happens first using the given
     * random number generator. Total bacterial load must be between 0 and 1000,
     * inclusive.
     *
     * @param targetFailures the number of failures to reach before stopping the
     *                       simulation
     * @param maximumRuns the maximum number of times to run the model, at least
     *                    1, must be greater or equal to the targetFailures.
     * @param initialLoad1 the initial bacterial load for strain 1. Between 0
     * and 1000, inclusive.
     * @param initialLoad2 the initial bacterial load for strain 2. Between 0
     * and 1000, inclusive.
     * @param random the random number generator, not null.
     * @return the antibiotic model.
     */
    public static AntibioticModel dynamicSampleSize(
            final int targetFailures,
            final int maximumRuns,
            final int initialLoad1,
            final int initialLoad2,
            final Random random) {
        if (targetFailures < 1) {
            throw new IllegalArgumentException(
                    "targetFailures = " + targetFailures + ", expected: >= 1.");
        }
        if (maximumRuns < 1) {
            throw new IllegalArgumentException(
                    "maximumRuns = " + maximumRuns + ", expected: >= 1.");
        }
        if (targetFailures > maximumRuns) {
            throw new IllegalArgumentException(
                    "targetFailures = " + targetFailures
                  + ", maximumRuns = " + maximumRuns
                  + ", expected: targetFailures <= maximumRuns.");
        }
        if (initialLoad1 < 0 || initialLoad1 > 1000) {
            throw new IllegalArgumentException(
                    "initialLoad1 = " + initialLoad1 + ", expected: [0, 1000]");
        }
        if (initialLoad2 < 0 || initialLoad2 > 1000) {
            throw new IllegalArgumentException(
                    "initialLoad2 = " + initialLoad2 + ", expected: [0, 1000]");
        }
        final int totalLoad = initialLoad1 + initialLoad2;
        if (totalLoad < 0 || totalLoad > 1000) {
            throw new IllegalArgumentException(
                    "initialLoad1 + initialLoad2 = "
                    + (initialLoad1 + initialLoad2) + ", expected: [0, 1000]");
        }
        if (random == null) {
            throw new NullPointerException("random");
        }
        return new AntibioticModel(-1,
                                   initialLoad1,
                                   initialLoad2,
                                   random,
                                   targetFailures,
                                   maximumRuns);
    }
    
    private AntibioticModel(
            final int runs,
            final int initialLoad1,
            final int initialLoad2,
            final Random random,
            final int targetFailures,
            final int maximumRuns) {
        this.runs = runs;
        this.initialLoad1 = initialLoad1;
        this.initialLoad2 = initialLoad2;
        this.random = random;
        this.targetFailures = targetFailures;
        this.maximumRuns = maximumRuns;
    }

    /**
     * Evaluates the candidate solution using a fitness evaluation based on a
     * weighting of the number of runs (samples) where bacteria is eradicated
     * and amount of antibiotic used. The fitness is returned, the input array
     * is not modified by the call. This implementation is consistent with the
     * MATLAB reference implementation.
     *
     * @param solution the dosage schedule, not null, length 10, all
     * non-negative.
     * @return the fitness from the stochastic model
     */
    public double evaluate(final int... solution) {

        // Validation and defensive copying of the input array
        final int[] x = validate(solution);

        // Calls evaluate implmentation for the copied solution
        final double[] result = __evaluate(true, x);

        // Unpack the returned values
        final double v_tot = result[0];
        final double highestConc = result[1];
        final double S1_Extinct = result[2];
        final double S2_Extinct = result[3];

        // Compute the fitness
        final double S1_Pen = 1 - S1_Extinct / runs;
        final double S2_Pen = 1 - S2_Extinct / runs;
        final double fitness = MATLAB_REF_W1 * 0.5 * (S1_Pen + S2_Pen)
                + MATLAB_REF_W2 * v_tot
                / MATLAB_REF_V_MAX;

        // Return penalty or fitness
        if (v_tot > MATLAB_REF_V_MAX) {
            return MATLAB_REF_PENALTY_V_MAX;
        } else if (highestConc > MATLAB_REF_MAX_CONC) {
            return MATLAB_REF_PENALTY_CONC;
        } else {
            return fitness;
        }

    }

    // Package-private implmentation of evaluate methods,
    // bypasses input validation
    // shortCircuit = skip evaluation if origonal constraints violated
    // x = the candidate
    double[] __evaluate(boolean shortCircuit, final int[] x) {

        // The total dosage
        final long v_tot = sum(x);

        // The duration of the treatment in days
        final int duration = duration(x);

        // Gets a random number generator
        final Random rng = (this.random != null)
                ? this.random : ThreadLocalRandom.current();

        final int tf = 15; // End Time

        // Time AB applied (vector)
        final int[] tint = new int[]{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, tf};

        // Concentration of AB applied (vector)
        final int[] cint = new int[]{0, x[0], x[1], x[2], x[3], x[4],
            x[5], x[6], x[7], x[8], x[9], 0};

        //System.err.println("cint = " + Arrays.toString(cint));
        
        //////// Deterministic for Concentration ////////
        //  Checks if AB concentration goes abve 60 at any point in treatment.
        final double[] conc2 = new double[11];
        for (int i = 1; i <= 10; i++) {
            conc2[i] = cint[i] + conc2[i - 1] * Math.exp(-a);
        }
        double highestConc = 0.0;
        for (int i = 2; i <= 11; i++) {
            if (conc2[i - 1] > highestConc) { // maxconc = 60 here
                highestConc = conc2[i - 1];
            }
        }
        //System.err.println("conc2 = " + Arrays.toString(conc2));
        //System.err.println("highestConc = " + highestConc);
        if (shortCircuit && ((highestConc > MATLAB_REF_MAX_CONC)
                || (v_tot > MATLAB_REF_V_MAX))) {
            return new double[]{v_tot,
                highestConc,
                Double.NaN,
                Double.NaN,
                Double.NaN,
                0.0,
                duration};
        }

        //////// Stochastic ////////

        int S1_Extinct = 0; // Initial Extinction for Strain 1
        int S2_Extinct = 0; // Initial Extinction for Strain 2
        int Both_Extinct = 0; // Initial Extinction for Both Strains
        
        final int[] bacterial_load = new int[2];
        
        int samplesTaken = 0;
        
        // Fixed sample size
        if (runs > -1) {
            
            for (int n = 0; n < runs; n++) {
                
                simulate(tint, cint, rng, bacterial_load);
                
                if (bacterial_load[0] < 1) {
                    S1_Extinct++;
                }
                if (bacterial_load[1] < 1) {
                    S2_Extinct++;
                }
                if (bacterial_load[0] < 1 && bacterial_load[1] < 1) {
                    Both_Extinct++;
                }
                
            }
            
            samplesTaken = runs;
        
        // Dynamic sample size
        } else {
            
            int failuresObserved = 0;
            while (samplesTaken < this.maximumRuns
                    && failuresObserved < this.targetFailures) {
                
                simulate(tint, cint, rng, bacterial_load);
                
                if (bacterial_load[0] < 1) {
                    S1_Extinct++;
                }
                if (bacterial_load[1] < 1) {
                    S2_Extinct++;
                }
                if (bacterial_load[0] < 1 && bacterial_load[1] < 1) {
                    Both_Extinct++;
                } else {
                    failuresObserved++;
                }
            
                samplesTaken++;
                
            }
            
        }

        // Returns the values needed to compute objectives
        return new double[]{v_tot,
            highestConc,
            S1_Extinct,
            S2_Extinct,
            Both_Extinct,
            samplesTaken,
            duration};

    }

    // Single simulation run, result returned in fourth parameter
    private void simulate(final int[] tint, final int[] cint, final Random rng,
                          final int[] result) {
        
        // Initial Conditions
        int S1 = this.initialLoad1; // Strain 1 least resistant
        int S2 = this.initialLoad2; // Strain 2 med resistant
        double C = 0; // Concentration of AB
        double time = 0.0; // Start Time

        // Does 1 run from day 0 to day 15
        for (int j = 1; j <= tint.length - 1; j++) {

            if (S1 < 1 && S2 < 1) {
                break;
            }

            C = C + cint[j - 1];
            final double C0 = C;
            final double timeAB = time;

            // pauses model each day to add next AB dose.
            while (time <= tint[j + 1 - 1]) {

                // Susceptible reproduction rate
                double rate0 = r * S1
                        * (1 - ((double) (S1 + S2) / K));

                // resistant reproduction rate
                double rate1 = r * S2
                        * (1 - ((double) (S1 + S2) / K)) * (1 - c1);

                // Susceptible death rate incl AB death
                double rate2 = ms
                        * S1 + ((MaxS - MinS)
                        * (Math.pow(C / MICS, kS)))
                        / ((Math.pow(C / MICS, kS))
                        - MinS / MaxS) * S1;

                // Resistant death rate incl AB death
                double rate3 = mr
                        * S2 + ((MaxR1 - MinR1)
                        * (Math.pow(C / MICR1, kR1)))
                        / ((Math.pow(C / MICR1, kR1))
                        - MinR1 / MaxR1) * S2;

                // Total rates
                double rateSum = rate0 + rate1 + rate2 + rate3;

                // Events
                final double x1 = rng.nextDouble();
                if (x1 <= rate0 / rateSum) {

                    // Strain 1 reproduces
                    S1 = S1 + 1;

                } else if (x1 <= (rate0 + rate1) / rateSum) {

                    // Strain 2 reproduces
                    S2 = S2 + 1;

                } else if (x1 <= (rate0 + rate1 + rate2) / rateSum) {

                    // Strain 1 dies naturally and AB
                    S1 = S1 - 1;

                } else {

                    // Strain 2 dies naturally and AB
                    S2 = S2 - 1;

                }

                // Generate time step for event
                time -= Math.log(rng.nextDouble()) / rateSum;

                // Concentration of AB
                C = C0 * Math.exp(-a * (time - timeAB));

                if (S1 < 1 && S2 < 1) {
                    break;
                }

            }

        }

        result[0] = S1;
        result[1] = S2;

    }

    // validation and defensive copying of the input array
    private int[] validate(final int[] solution) {
        if (solution == null) {
            throw new NullPointerException("solution");
        }
        if (solution.length != 10) {
            throw new IllegalArgumentException(
                    "candidate.length = " + solution.length
                    + ", expected: 10.");
        }
        final int[] x = solution.clone();
        for (int i = 0; i < x.length; i++) {
            if (x[i] < 0) {
                throw new IllegalArgumentException(
                        "solution[" + i + "] = " + x[i] + ", expected: >= 0.");
            }
        }
        return x;
    }

    // calculates the sum of an array of integers
    private static long sum(final int[] x) {
        long rv = 0; // long, to prevent overflow
        for (int i = 0; i < x.length; i++) {
            rv += x[i];
        }
        return rv;
    }

    // calculates the duration of the treatment in days
    private static int duration(final int[] x) {
        int lastNonZeroIndex = -1;
        for (int i = 0; i < x.length; i++) {
            if (x[i] != 0) {
                lastNonZeroIndex = i;
            }
        }
        return lastNonZeroIndex + 1;
    }

}
