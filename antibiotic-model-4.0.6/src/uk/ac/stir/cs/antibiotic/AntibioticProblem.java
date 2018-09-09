/*
 * AntibioticProblem.java
 *
 * Copyright © 2017-2018
 * University of Stirling
 * All Rights Reserved
 */

package uk.ac.stir.cs.antibiotic;

import java.util.*;
import org.uma.jmetal.problem.ConstrainedProblem;
import org.uma.jmetal.problem.impl.*;
import org.uma.jmetal.solution.*;
import org.uma.jmetal.solution.impl.*;
import org.uma.jmetal.util.solutionattribute.impl.NumberOfViolatedConstraints;
import org.uma.jmetal.util.solutionattribute.impl.OverallConstraintViolation;

/**
 * jMetal adapter class for AntibioticModel.
 *
 * @author Lee A. Christie
 * 
 * @version 4.0.6
 */
public class AntibioticProblem
        extends AbstractIntegerProblem {

    private static final int MAX_LENGTH = 10;
    private static final long serialVersionUID = 04_00_00L;
    
    private final AntibioticModel model;
    private final int maxIndividualDosage;
    private final int minStart;
    private final double maxFailure;
    private final int length;
    private final AntibioticObjective[] objectives;
    
    public static AntibioticProblemBuilder builder(
            final AntibioticModel model,
            final AntibioticObjective firstObjective) {
        return new AntibioticProblemBuilder(model, firstObjective);
    }
    
    /**
     * Adapts an antibiotic model to a single- or multi-objective optimisation
     * problem using the jMetal IntegerProblem interface.
     *
     * @param model the antibiotic model, not null.
     * @param maxIndividualDosage the maximum dosage, which specified the limit
     *                            on variable values, &gt;= 1;
     * @param length the maximum length of treatment 1 to 10
     * @param minStart the minimum initial dosage, 0 to maxIndividualDosage
     * @param maxFailure the maximum failure rate 0.0 to 1.0, maxFailure rate of
     *                   1.0 disables the constraint
     * @param objective0 objective[0], not null.
     * @param additionalObjectives objective[1] ... objective[n-1], not null,
     *                             these additional objectives are not specified
     *                             if creating a single objective problem.
     */
    private AntibioticProblem(
            final AntibioticModel model,
            final int maxIndividualDosage,
            final int length,
            final int minStart,
            final double maxFailure,
            final AntibioticObjective objective0,
            final AntibioticObjective... additionalObjectives) {
        if (model == null) {
            throw new NullPointerException("model");
        }
        if (maxIndividualDosage < 1) {
            throw new IllegalArgumentException("maxIndividualDosage = "
                    + maxIndividualDosage + ", expected >= 1.");
        }
        if (minStart < 0 || minStart > maxIndividualDosage) {
            throw new IllegalArgumentException("minStart = "
                    + minStart
                    + ", expected 0 to maxIndividualDosage.");
        }
        if (length < 1 || length > MAX_LENGTH) {
            throw new IllegalArgumentException("length = "
                    + length + ", expected 1 to " + MAX_LENGTH + ".");
        }
        if (!Double.isFinite(maxFailure) || maxFailure < 0 || maxFailure > 1) {
            throw new IllegalArgumentException("maxFailure = "
                    + maxFailure
                    + ", expected 0.0 to 1.0.");
        }
        if (objective0 == null) {
            throw new NullPointerException("objective 0");
        }
        if (additionalObjectives == null) {
            throw new NullPointerException("objective 1");
        }
        this.model = model;
        this.maxIndividualDosage = maxIndividualDosage;
        this.minStart = minStart;
        this.maxFailure = maxFailure;
        this.length = length;
        this.objectives = clone(objective0, additionalObjectives);
    }
    
    /**
     * Returns the upper bound on a solution variable (specified by
     * {@code maxindividualDosage} in the constructor).
     *
     * @param index the index, between 0 and (length-1) inclusive
     * @return the upper bound
     */
    @Override
    public Integer getUpperBound(int index) {
        if (index < 0 || index >= length) {
            throw new NoSuchElementException(
                    "index = " + index + ", expected 0 to "
                            + (length - 1) + ".");
        }
        return maxIndividualDosage;
    }

    /**
     * Returns the upper bound on a solution variable ({@code 0}).
     *
     * @param index the index, between 0 and (length-1) inclusive
     * @return 0
     */
    @Override
    public Integer getLowerBound(int index) {
        if (index < 0 || index >= length) {
            throw new NoSuchElementException(
                    "index = " + index + ", expected 0 to "
                            + (length - 1) + ".");
        }
        if (index == 0) {
            return minStart;
        }
        return 0;
    }

    /**
     * Not supported.
     *
     * @param lowerLimit not used
     */
    @Override
    protected void setLowerLimit(List<Integer> lowerLimit) {
        throw new UnsupportedOperationException("Not supported.");
    }

    /**
     * Not supported.
     *
     * @param upperLimit not used
     */
    @Override
    protected void setUpperLimit(List<Integer> upperLimit) {
        throw new UnsupportedOperationException("Not supported.");
    }

    /**
     * Creates a solution object for this problem.
     *
     * @return a new solution object
     */
    @Override
    public IntegerSolution createSolution() {
        return new DefaultIntegerSolution(this);
    }
    
    /**
     * Returns the name of the problem.
     *
     * @return the name of the problem
     */
    @Override
    public String getName() {
        return this.getClass().getSimpleName();
    }
    
    /**
     * Returns the number of constraints (0).
     *
     * @return 0
     */
    @Override
    public int getNumberOfConstraints() {
        return (this.maxFailure < 1.0) ? 1 : 0;
    }
    
    /**
     * Returns the number of objectives, specified by the constructor.
     *
     * @return the number of objectives
     */
    @Override
    public int getNumberOfObjectives() {
        return objectives.length;
    }
    
    /**
     * Returns the number of variables (length).
     *
     * @return the length
     */
    @Override
    public int getNumberOfVariables() {
        return length;
    }

    /**
     * Evaluates the candidate solution using this problem's specified
     * objectives. The fitness will be written to objective zero in the given
     * candidate solution.
     *
     * @param solution the dosage schedule, not null, length must be equal to
     *                 the problem length, all non-negative.
     * @return the objective fitness values
     */
    public double[] evaluateMultiObjective(int... solution) {
    
        // Validation, defensive copying, and padding of the input array
        final int[] x = validateCopyAndPad(solution);
        
        // Calls evaluate implmentation for the copied solution
        final double[] result = model.__evaluate(false, x);
        
        // Assign the fitness for each obejctive
        final double[] rv = new double[objectives.length];
        for (int i = 0; i < objectives.length; i++) {
            rv[i] = objectives[i].compute(result);
        }
        return rv;
        
    }

    /**
     * Evaluates the candidate solution using this problem's specified
     * objectives. The fitness will be written to objective zero in the given
     * candidate solution.
     *
     * @param solution the dosage schedule, not null, length must be equal to
     *                 the problem length, all non-negative.
     * @return the first objective fitness value
     */
    public double evaluateFirstObjective(int... solution) {
    
        // Validation, defensive copying, and padding of the input array
        final int[] x = validateCopyAndPad(solution);
        
        // Calls evaluate implmentation for the copied solution
        final double[] result = model.__evaluate(false, x);
        
        // Assign the fitness for each obejctive
        return objectives[0].compute(result);
        
    }

    /**
     * Evaluates the candidate solution using this problem's specified
     * objectives. The fitness will be written to objective zero in the given
     * candidate solution.
     *
     * @param solution the dosage schedule, not null, length must be equal to
     *                 the problem length, all non-negative.
     */
    @Override
    public void evaluate(final IntegerSolution solution) {
    
        // Validation and defensive copying of the input array
        final int[] x = validateCopyAndPad(solution);
        
        // Calls evaluate implmentation for the copied solution
        final double[] result = model.__evaluate(false, x);
        
        // Assign the fitness for each obejctive
        for (int i = 0; i < objectives.length; i++) {
            solution.setObjective(i, objectives[i].compute(result));
        }
        
        // Record the actual number of samples taken
        solution.setAttribute("samples", (int) result[5]);
        
        // Set solution as evaluated
        solution.setAttribute("evaluated", true);
        
    }
    
    // validation, defensive copying, and padding of the input array
    private int[] validateCopyAndPad(int[] solution) {
        if (solution == null) {
            throw new NullPointerException("solution");
        }
        if (solution.length != length) {
            throw new IllegalArgumentException(
                    "solution.length() = "
                    + solution.length + ", expected: " + length + ".");
        }
        final int[] x = new int[MAX_LENGTH];
        for (int i = 0; i < solution.length; i++) {
            x[i] = solution[i];
            if (x[i] < 0) {
                throw new IllegalArgumentException(
                        "solution[" + i + "] = " + x[i] + ", expected: >= 0.");
            }
        }
        return x;
    }
    
    // validation and defensive copying of the input jMetal solution
    private int[] validateCopyAndPad(IntegerSolution solution) {
        if (solution == null) {
            throw new NullPointerException("solution");
        }
        final int numVars = solution.getNumberOfVariables();
        if (numVars != length) {
            throw new IllegalArgumentException(
                    "solution.getNumberOfVariables() = "
                    + numVars + ", expected: " + length + ".");
        }
        final int[] x = new int[MAX_LENGTH];
        for (int i = 0; i < solution.getNumberOfVariables(); i++) {
            x[i] = solution.getVariableValue(i);
            if (x[i] < 0) {
                throw new IllegalArgumentException(
                        "solution[" + i + "] = " + x[i] + ", expected: >= 0.");
            }
        }
        return x;
    }

    // clone the objective array from var-args input
    private AntibioticObjective[] clone(AntibioticObjective obj0,
                                        AntibioticObjective[] obj) {
        AntibioticObjective[] rv
                = new AntibioticObjective[obj.length + 1];
        rv[0] = obj0;
        for (int i = 0; i < obj.length; i++) {
            rv[i+1] = obj[i];
            if (rv[i+1] == null) {
                throw new NullPointerException("objectives[" + (i + 1) + "]");
            }
        }
        return rv;
    }

    /**
     * Builder object for AntibioticProblem. Use
     * {@link #builder(uk.ac.stir.cs.antibiotic.AntibioticModel, uk.ac.stir.cs.antibiotic.AntibioticObjective) }
     * to create a builder instance.
     */
    public static final class AntibioticProblemBuilder {

        private final AntibioticModel model;
        private final List<AntibioticObjective> objectives;
        
        private int treatmentDuration = 10;
        private int minimumInitialDosage = 0;
        private int maximumIndividualDosage = 60;
        private double maxFailure = 1.0;
        
        private AntibioticProblemBuilder(
            final AntibioticModel model,
            final AntibioticObjective firstObjective) {
            this.model = model;
            objectives = new ArrayList<>(3);
            objectives.add(firstObjective);
        }
        
        /**
         * Sets the maximum individual dosage. If not called, the default is
         * 60.
         * 
         * @param maximumIndividualDosage the maximum individual dosage
         * @return this builder
         */
        public AntibioticProblemBuilder
                setMaximumIndividualDosage(int maximumIndividualDosage) {
            this.maximumIndividualDosage = maximumIndividualDosage;
            return this;
        }
                
        /**
         * Sets the treatment duration. If not called, the default is 10.
         *
         * @param treatmentDuration the treatment duration.
         * @return this builder
         */
        public AntibioticProblemBuilder
                setMaximumTreatmentDuration(int treatmentDuration) {
            this.treatmentDuration = treatmentDuration;
            return this;
        }
                
        /**
         * Sets the minimum initial dosage. If not called, the default is 0.
         *
         * @param minimumInitialDosage the minimum initial dosage
         * @return this builder
         */
        public AntibioticProblemBuilder
                setMinimumInitialDosage(int minimumInitialDosage) {
            this.minimumInitialDosage = minimumInitialDosage;
            return this;
        }
                
        /**
         * Sets the maximum failure rate. If not called, the default is 1.0.
         *
         * @param maxFailure the maximum failure rate
         * @return this builder
         */
        public AntibioticProblemBuilder
                setMaximumFailureRate(double maxFailure) {
            this.maxFailure = maxFailure;
            return this;
        }
                
        /**
         * Adds another objective.
         *
         * @param objective the objective
         * @return this builder
         */
        public AntibioticProblemBuilder
                appendObjective(AntibioticObjective objective) {
            this.objectives.add(objective);
            return this;
        }
        
        /**
         * Builds the problem.
         *
         * @return the problem
         */
        public AntibioticProblem build() {
            AntibioticObjective[] restOfObjectives
                    = new AntibioticObjective[objectives.size()-1];
            for (int i = 1; i < objectives.size(); i++) {
                restOfObjectives[i-1] = objectives.get(i);
            }
            if (maxFailure < 1.0) {
                return new ConstrainedAntibioticProblem(
                        model,
                        maximumIndividualDosage,
                        treatmentDuration,
                        minimumInitialDosage,
                        maxFailure,
                        objectives.get(0),
                        restOfObjectives);
            } else {
                return new AntibioticProblem(
                        model,
                        maximumIndividualDosage,
                        treatmentDuration,
                        minimumInitialDosage,
                        maxFailure,
                        objectives.get(0),
                        restOfObjectives);
            }
        }
                
    }
    
    private static class ConstrainedAntibioticProblem
            extends AntibioticProblem
            implements ConstrainedProblem<IntegerSolution> {
        
        private OverallConstraintViolation<IntegerSolution>
                overallConstraintViolationDegree
                = new OverallConstraintViolation<IntegerSolution>();
        private NumberOfViolatedConstraints<IntegerSolution>
                numberOfViolatedConstraints
                = new NumberOfViolatedConstraints<IntegerSolution>();
        
        private int failureRateObjectiveIndex = -1;
        
        private static final long serialVersionUID = 04_00_00L;
        
        private ConstrainedAntibioticProblem(
                final AntibioticModel model,
                final int maxIndividualDosage,
                final int length,
                final int minStart,
                final double maxFailure,
                final AntibioticObjective objective0,
                final AntibioticObjective... additionalObjectives) {
            super(model,
                    maxIndividualDosage,
                  length,
                  minStart,
                  maxFailure,
                  objective0,
                  additionalObjectives);
            if (objective0 ==
                    AntibioticObjective.UNCURED_PROPORTION_SINGLETON) {
                failureRateObjectiveIndex = 0;
            }
            if (failureRateObjectiveIndex == -1) {
                for (int i = 0; i < additionalObjectives.length; i++) {
                    if (additionalObjectives[i] ==
                            AntibioticObjective.UNCURED_PROPORTION_SINGLETON) {
                        failureRateObjectiveIndex = i + 1;
                    }
                }
            }
        }        

        @Override
        public void evaluateConstraints(IntegerSolution solution) {
            
            double failureRate = Double.NaN;
            
            // Check for evaluated value
            if (failureRateObjectiveIndex != -1) {
                Object evaluated = null;
                try {
                    evaluated = solution.getAttribute("evaluated");
                } catch (RuntimeException ignored) {}
                if (evaluated != null & evaluated.equals(true)) {
                    try {
                        failureRate = solution.getObjective(
                                failureRateObjectiveIndex);
                    } catch (RuntimeException ignored) {}
                } else {
                    throw new IllegalStateException(
                            "Solution objectives must be evaluated "
                                    + "to evaluate constraints.");
                }
            }
            
            // Evaluated only in the case that failure rate is not an objective
            if (!Double.isFinite(failureRate)
                    || failureRate < 0.0 || failureRate > 1.0) {
    
                // Validation and defensive copying of the input array
                final int[] x = super.validateCopyAndPad(solution);

                // Calls evaluate implmentation for the copied solution
                final double[] result = super.model.__evaluate(false, x);
                
                // Gets the failure rate
                failureRate = AntibioticObjective
                        .UNCURED_PROPORTION_SINGLETON.compute(result);
        
            }
            
            // Sets the constraints
            overallConstraintViolationDegree.setAttribute(
                    solution, Math.min(super.maxFailure - failureRate, 0.0));
            numberOfViolatedConstraints.setAttribute(
                    solution, (failureRate > super.maxFailure) ? 1 : 0);
            
        }
        
    }
    
}
