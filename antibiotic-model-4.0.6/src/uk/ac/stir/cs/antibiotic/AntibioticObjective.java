/*
 * AntibioticObjective.java
 *
 * Copyright © 2017-2018
 * University of Stirling
 * All Rights Reserved
 */

package uk.ac.stir.cs.antibiotic;

import java.util.function.DoublePredicate;

/**
 * Defines objectives to build single or multi-objective problem from the
 * antibiotic model.
 *
 * @author Lee A. Christie
 * 
 * @version 4.0.6
 */
@FunctionalInterface
public interface AntibioticObjective {
    
    static final AntibioticObjective UNCURED_PROPORTION_SINGLETON
    = result -> {
        final double Both_Extinct = result[4];
        final double runs = result[5];
        return 1.0 - (Both_Extinct / runs);
    };
    
    /**
     * Computes the objective from the values given by the model.
     * 
     * @param result the values returned by the stochastic model
     * @return the objective value
     */
    double compute(double[] result);
    
    /**
     * The proportion of the patients which are uncured in the stochastic
     * simulation. Objective will return real values between 0.0 and 1.0,
     * inclusive.
     *
     * @return the objective, not null
     */
    public static AntibioticObjective uncuredProportion() {
        return UNCURED_PROPORTION_SINGLETON;
    }
    
    /**
     * The amount by which the actual maximum concentration of antibiotic at any
     * one time exceeds the specified safe limit, or 0.0 if the limit is never
     * exceeded. Objective will return non-negative real values.
     *
     * @param limit the maximum safe concentration of antibiotic at any one
     *              time.
     * @return the objective, not null
     */
    public static AntibioticObjective overdoseAmount(final double limit) {
        return result -> {
            final double highestConc = result[1];
            if (highestConc < limit) {
                return 0.0;
            } else {
                return highestConc - limit;
            }
        };
    }
    
    /**
     * The maximum concentration of antibiotic at any one time exceeds.
     * Objective will return non-negative real values. This is equivalent to
     * {@link #overdoseAmount(double)} for {@code limit = 0}. 
     *
     * @return the objective, not null
     */
    public static AntibioticObjective maximumConcentration() {
        return result -> result[1];
    }
    
    /**
     * The number of days for the complete treatment. For example, the treatment
     * {@code 0, 0, 0, 10, 10, 10, 10, 0, 0, 0} is 7 days in duration, because
     * the last non-zero dosage is at day 7. Objective will return non-negative
     * real values.
     *
     * @return the objective, not null
     */
    public static AntibioticObjective treatmentDuration() {
        return result -> result[6];
    }
    
    /**
     * The total amount of antibiotic administered. Objective will return
     * non-negative integers.
     *
     * @return the objective, not null
     */
    public static AntibioticObjective totalAntibiotic() {
        return result -> result[0];
    }
    
    /**
     * The fitness function used by the MATLAB reference implementation.
     *
     * @return the objective, not null
     */
    public static AntibioticObjective weighting() {
        return weighting(AntibioticModel.MATLAB_REF_W1,
                         AntibioticModel.MATLAB_REF_W2);
    }
    
    /**
     * The fitness function used by the MATLAB reference implementation, but
     * with different specified weighting value.
     * 
     * @param w1 the weight on the first part of the function
     * @param w2 the weight on the second part of the function
     * @return the objective, not null
     */
    public static AntibioticObjective weighting(final double w1,
                                                final double w2) {
        return result -> {
            final double v_tot = result[0];
            final double highestConc = result[1];
            final double S1_Extinct = result[2];
            final double S2_Extinct = result[3];
            final double runs = result[5];
            final double S1_Pen = 1 - S1_Extinct / runs;
            final double S2_Pen = 1 - S2_Extinct / runs;
            final double fitness
                    = w1 * 0.5 * (S1_Pen + S2_Pen)
                    + w2 * v_tot
                    / AntibioticModel.MATLAB_REF_V_MAX;
            if (v_tot > AntibioticModel.MATLAB_REF_V_MAX) {
                return AntibioticModel.MATLAB_REF_PENALTY_V_MAX;
            } else if (highestConc > AntibioticModel.MATLAB_REF_MAX_CONC) {
                return AntibioticModel.MATLAB_REF_PENALTY_CONC;
            } else {
                return fitness;
            }
        };
    }

    /**
     * The fitness is specified by summing other given sub-objectives.
     *
     * @param objectives the sub-objectives, not null, no null elements
     * @return the objective, not null
     */
    public static AntibioticObjective sum(
            final AntibioticObjective... objectives) {
        if (objectives == null) {
            throw new NullPointerException("obejctives");
        }
        final AntibioticObjective[] clone = objectives.clone();
        for (int i = 0; i < clone.length; i++) {
            if (clone[i] == null) {
                throw new NullPointerException("obejctives[" + i + "]");
            }
        }
        return result -> {
            double rv = 0.0;
            for (AntibioticObjective objective : clone) {
                rv += objective.compute(result);
            }
            return rv;
        };
    }

    /**
     * The fitness is specified by multiplying the specified sub-objectives by a
     * scalar.
     *
     * @param objective the sub-objective, not null
     * @param factor the scalar multiple
     * @return the objective, not null
     */
    public static AntibioticObjective multiply(
            final AntibioticObjective objective,
            final double factor) {
        if (objective == null) {
            throw new NullPointerException("obejctive");
        }
        return result -> objective.compute(result) * factor;
    }

    /**
     * Adds a second objective if the first matches a predicate.
     *
     * @param initialObjective the sub-objective, not null
     * @param predicate the predicate to test on the first sub-objective, not
     *                  null
     * @param objectiveToAdd the sub-objective to add, not null
     * @return the objective, not null
     */
    public static AntibioticObjective conditionalAdd(
            final AntibioticObjective initialObjective,
            final DoublePredicate predicate,
            final AntibioticObjective objectiveToAdd) {
        if (initialObjective == null) {
            throw new NullPointerException("initialObjective");
        }
        if (predicate == null) {
            throw new NullPointerException("predicate");
        }
        if (objectiveToAdd == null) {
            throw new NullPointerException("objectiveToAdd");
        }
        return result -> (
                predicate.test(initialObjective.compute(result))
              ? initialObjective.compute(result)
                      + objectiveToAdd.compute(result)
              : initialObjective.compute(result));
    }
    
}
