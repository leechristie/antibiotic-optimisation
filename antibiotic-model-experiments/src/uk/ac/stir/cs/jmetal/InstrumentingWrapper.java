/*
 * InstrumentingWrapper.java
 *
 * Copyright © 2017-2018
 * University of Stirling
 * All Rights Reserved
 */

package uk.ac.stir.cs.jmetal;

import org.uma.jmetal.problem.ConstrainedProblem;
import org.uma.jmetal.problem.IntegerProblem;
import org.uma.jmetal.solution.IntegerSolution;
import org.uma.jmetal.util.solutionattribute.impl.NumberOfViolatedConstraints;
import org.uma.jmetal.util.solutionattribute.impl.OverallConstraintViolation;

/*
 * A wrapper class to count evlauations and samples. To be used with
 * {@link uk.ac.stir.cs.antibiotic.AntibioticProblem}. May be used with
 * constrained or unconstrained problems.
 *
 * @author Lee A. Christie
 *
 * @version 4.0.6
 */
public final class InstrumentingWrapper {

    private InstrumentingWrapper() {
        throw new AssertionError("Static utility class constructor.");
    }
    
    /**
     * Adds a wrapper on the integer problem which counts evaluations and
     * samples for the antibiotic problem.
     *
     * @param problem the problem
     * @param expectedLimit the expected evaluation limit
     * @return the wrapped problem
     */
    public static IntegerProblem instrument(
            final IntegerProblem problem,
            final int expectedLimit) {
        return instrument(problem, expectedLimit, false);
    }
    
    /**
     * Adds a wrapper on the integer problem which counts evaluations and
     * samples for the antibiotic problem.
     *
     * @param problem the problem
     * @param expectedLimit the expected evaluation limit
     * @param verbose return extra information
     * @return the wrapped problem
     */
    public static IntegerProblem instrument(
            final IntegerProblem problem,
            final int expectedLimit,
            final boolean verbose) {
        if (problem instanceof ConstrainedProblem<?>) {
            return new InstrumentedConstrainedProblem(problem, expectedLimit, verbose);
        } else {
            return new InstrumentedProblem(problem, expectedLimit, verbose);
        }
    }
    
    private static class InstrumentedConstrainedProblem
            implements IntegerProblem, ConstrainedProblem<IntegerSolution> {
        private static final long serialVersionUID = 04_00_00L;
        private final IntegerProblem problem;
        private final int expectedLimit;
        private int evals = 0;
        private long totalSamples = 0;
        private boolean verbose;
        public InstrumentedConstrainedProblem(final IntegerProblem problem,
                                              final int expectedLimit,
                                              final boolean verbose) {
            if (!(problem instanceof ConstrainedProblem<?>)) {
                throw new AssertionError();
            }
            this.problem = problem;
            this.expectedLimit = expectedLimit;
            this.verbose = verbose;
        }   
        @Override
        public int getNumberOfVariables() {
                return problem.getNumberOfVariables();
        }
        @Override
        public int getNumberOfObjectives() {
                return problem.getNumberOfObjectives();
        }
        @Override
        public int getNumberOfConstraints() {
                return problem.getNumberOfConstraints();
        }
        @Override
        public String getName() {
                return problem.getName();
        }
        @Override
        public void evaluate(IntegerSolution solution) {
            problem.evaluate(solution);
            if (verbose) {
                System.out.println("Variables :");
                for (int i = 0; i < solution.getNumberOfVariables(); i++) {
                    System.out.print("\t" + solution.getVariableValue(i));
                }
                System.out.println();
            }
            int samples = (Integer) solution.getAttribute("samples");
            totalSamples += samples;
            System.out.println("evaluation #"
                    + (++evals) + " of " + expectedLimit
                    + "\t samples so far : " + totalSamples);
            if (verbose) {
                System.out.print("Objectives :");
                for (int i = 0; i < solution.getNumberOfObjectives(); i++) {
                    System.out.print("\t" + solution.getObjective(i));
                }
                System.out.println();
            }
        }
        @Override
        public IntegerSolution createSolution() {
                return problem.createSolution();
        }
        @Override
        public Integer getLowerBound(int index) {
            return problem.getLowerBound(index);
        }
        @Override
        public Integer getUpperBound(int index) {
            return problem.getUpperBound(index);
        }
        @Override
        public void evaluateConstraints(IntegerSolution solution) {
            ((ConstrainedProblem<IntegerSolution>) problem)
                    .evaluateConstraints(solution);
            if (verbose) {
                System.out.print("OverallConstraintViolation :");
                System.out.println("\t" + solution.getAttribute(OverallConstraintViolation.class));
                System.out.print("NumberOfViolatedConstraints :");
                System.out.println("\t" + solution.getAttribute(NumberOfViolatedConstraints.class));
                System.out.println();
            }
        }
    }
    
    private static class InstrumentedProblem
            implements IntegerProblem {
        private static final long serialVersionUID = 04_00_00L;
        private final IntegerProblem problem;
        private final int expectedLimit;
        private int evals = 0;
        private long totalSamples = 0;
        private boolean verbose;
        public InstrumentedProblem(final IntegerProblem problem,
                                   final int expectedLimit,
                                   final boolean verbose) {
            if (problem instanceof ConstrainedProblem<?>) {
                throw new AssertionError();
            }
            this.problem = problem;
            this.expectedLimit = expectedLimit;
            this.verbose = verbose;
        }   
        @Override
        public int getNumberOfVariables() {
                return problem.getNumberOfVariables();
        }
        @Override
        public int getNumberOfObjectives() {
                return problem.getNumberOfObjectives();
        }
        @Override
        public int getNumberOfConstraints() {
                return problem.getNumberOfConstraints();
        }
        @Override
        public String getName() {
                return problem.getName();
        }
        @Override
        public void evaluate(IntegerSolution solution) {
            problem.evaluate(solution);
            if (verbose) {
                System.out.println("Variables :");
                for (int i = 0; i < solution.getNumberOfVariables(); i++) {
                    System.out.print("\t" + solution.getVariableValue(i));
                }
                System.out.println();
            }
            int samples = (Integer) solution.getAttribute("samples");
            totalSamples += samples;
            System.out.println("evaluation #"
                    + (++evals) + " of " + expectedLimit
                    + "\t samples so far : " + totalSamples);
            if (verbose) {
                System.out.print("Objectives :");
                for (int i = 0; i < solution.getNumberOfObjectives(); i++) {
                    System.out.print("\t" + solution.getObjective(i));
                }
                System.out.println();
                System.out.println();
            }
        }
        @Override
        public IntegerSolution createSolution() {
                return problem.createSolution();
        }
        @Override
        public Integer getLowerBound(int index) {
            return problem.getLowerBound(index);
        }
        @Override
        public Integer getUpperBound(int index) {
            return problem.getUpperBound(index);
        }
    }
    
}
