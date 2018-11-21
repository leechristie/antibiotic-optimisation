package uk.ac.stir.cs.antibiotic.experiment;

import java.util.*;
import org.uma.jmetal.algorithm.*;
import org.uma.jmetal.algorithm.multiobjective.nsgaii.*;
import org.uma.jmetal.operator.*;
import org.uma.jmetal.operator.impl.crossover.*;
import org.uma.jmetal.operator.impl.mutation.*;
import org.uma.jmetal.operator.impl.selection.*;
import org.uma.jmetal.problem.IntegerProblem;
import org.uma.jmetal.problem.Problem;
import org.uma.jmetal.solution.*;
import org.uma.jmetal.util.*;
import org.uma.jmetal.util.comparator.*;
import uk.ac.stir.cs.antibiotic.*;
import static org.uma.jmetal.runner.AbstractAlgorithmRunner.*;
import static uk.ac.stir.cs.antibiotic.AntibioticObjective.*;
//import static uk.ac.stir.cs.jmetal.InstrumentingWrapper.instrument;
import org.uma.jmetal.util.fileoutput.SolutionListOutput;
import org.uma.jmetal.util.fileoutput.impl.DefaultFileOutputContext; 

public class Main {

    public static void main(String[] args) throws Exception {
        
        //System.out.println(Arrays.toString(args));    
    
        int length = Integer.parseInt(arg(args, "length", "10"));
        System.err.println("length = " + length);
    
        int runs = Integer.parseInt(arg(args, "runs", "1000"));
        System.err.println("runs = " + runs);
    
        String seedString = arg(args, "seed", "none");
        long seed = -1;
        if ("none".equalsIgnoreCase(seedString)) {
            System.err.println("seed = " + seedString);
        } else {
            seed = Long.parseLong(seedString);
            System.err.println("seed = " + seed);
        }
        
        int maxEvals = Integer.parseInt(arg(args, "evals", "25000"));
        System.err.println("evals = " + maxEvals);
    
        int popSize = Integer.parseInt(arg(args, "popsize", "100"));
        System.err.println("popsize = " + popSize);
    
        int strain1 = Integer.parseInt(arg(args, "strain1", "700"));
        System.err.println("strain1 = " + strain1);
    
        int strain2 = Integer.parseInt(arg(args, "strain2", "100"));
        System.err.println("strain2 = " + strain2);
    
        String sampletype = arg(args, "sampletype", "fixed");
        System.err.println("sampletype = " + sampletype);
    
        int minstart = Integer.parseInt(arg(args, "minstart", "0"));
        System.err.println("minstart = " + minstart);
    
        int maxIndividualDosage = Integer.parseInt(arg(args, "dosage", "60"));
        System.err.println("dosage = " + maxIndividualDosage);
        
        // Define the model
        final AntibioticModel model;
        if ("fixed".equalsIgnoreCase(sampletype)) {
            if ("none".equalsIgnoreCase(seedString)) {
                model = AntibioticModel.fixedSampleSize(runs, strain1, strain2);
            } else {
                model = AntibioticModel.fixedSampleSize(runs, strain1, strain2, new Random(seed));
            }
        } else if ("dynamic".equalsIgnoreCase(sampletype)) {
            int failures = Integer.parseInt(arg(args, "failures", "100"));
            System.err.println("failures = " + failures);
            if ("none".equalsIgnoreCase(seedString)) {
                model = AntibioticModel.dynamicSampleSize(failures, runs, strain1, strain2);
            } else {
                model = AntibioticModel.dynamicSampleSize(failures, runs, strain1, strain2, new Random(seed));
            }
        } else {
            throw new IllegalArgumentException("sampletype expected fixed or dynamic");
        }
        
        String[] secondobjective = arg(args, "secondobjective", "maximumconcentration").split(",");
        List<String> secondObjectiveList = Arrays.asList(secondobjective);
        System.err.println("secondobjective = " + secondObjectiveList);
    
        double maxFail = Double.parseDouble(arg(args, "maxfail", "1.0"));
        System.err.println("maxfail = " + maxFail);

        // Define the jMetal problem to optimise the model
        AntibioticProblem.AntibioticProblemBuilder builder
                = AntibioticProblem.builder(model, uncuredProportion())
                                   .setMaximumIndividualDosage(maxIndividualDosage)
                                   .setMaximumTreatmentDuration(length)
                                   .setMaximumFailureRate(maxFail)
                                   .setMinimumInitialDosage(minstart);
        for (String o : secondObjectiveList) {
            if (o.equalsIgnoreCase("maximumconcentration")) {
                builder.appendObjective(maximumConcentration());
            } else if (o.equalsIgnoreCase("totalantibiotic")) {
                builder.appendObjective(totalAntibiotic());
            } else if (o.equalsIgnoreCase("treatmentduration")) {
                builder.appendObjective(treatmentDuration());
            } else {
                throw new IllegalArgumentException("Unknown objective : " + o);
            }
        }
        
        //final Problem<IntegerSolution> problem
        //        = instrument(builder.build(), maxEvals);
        final Problem<IntegerSolution> problem = builder.build();

        // Define the jMetal Algorithm
        CrossoverOperator<IntegerSolution> crossover;
        double mutationProbability = 1.0 / problem.getNumberOfVariables();
        double mutationDistributionIndex = 20.0;
        MutationOperator<IntegerSolution> mutation
                = new IntegerPolynomialMutation(mutationProbability,
                        mutationDistributionIndex);
        SelectionOperator<List<IntegerSolution>, IntegerSolution> selection
                = new BinaryTournamentSelection<>(
                new RankingAndCrowdingDistanceComparator<>());
        double crossoverProbability = 0.9;
        double crossoverDistributionIndex = 20.0;
        crossover = new IntegerSBXCrossover(
                crossoverProbability, crossoverDistributionIndex);
        final NSGAIIBuilder<IntegerSolution> nsagabuilder
                = new NSGAIIBuilder<>(problem, crossover, mutation)
                        .setSelectionOperator(selection)
                        .setMaxEvaluations(maxEvals)
                        .setPopulationSize(popSize);
//        builder.setSolutionListEvaluator(
//                builder.getSolutionListEvaluator());
        final Algorithm<List<IntegerSolution>> algorithm = nsagabuilder.build();
        
        // Run the algorithm
        AlgorithmRunner algorithmRunner
                = new AlgorithmRunner.Executor(algorithm)
                                     .execute() ;
        
        // Print the results
        List<IntegerSolution> population = algorithm.getResult() ;
        long computingTime = algorithmRunner.getComputingTime() ;
        JMetalLogger.logger.info(
                "Total execution time: " + computingTime + "ms");
        printFinalSolutionSet(population);
        
        new SolutionListOutput(population)
        .setSeparator("\t")
        .setVarFileOutputContext(new DefaultFileOutputContext("VAR.tsv"))
        .setFunFileOutputContext(new DefaultFileOutputContext("FUN.tsv"))
        .print();
        JMetalLogger.logger.info("Objectives values have been written to file FUN.tsv");
        JMetalLogger.logger.info("Variables values have been written to file VAR.tsv"); 
        
    }

    private static String arg(String[] args, String name) {
        return arg(args, name, null);
    }
    
    private static String arg(String[] args, String name, String defaultString) {
        int index = -1;
        for (int i = 0; i < args.length; i++) {
            if (args[i].equalsIgnoreCase("--" + name)
                    || args[i].equalsIgnoreCase("-" + name)
                    || args[i].equalsIgnoreCase(name)) {
                index = i + 1;
                break;
            }
        }
        if (index == -1) {
            if (defaultString == null) {
                System.err.println("Missing argument: --" + name);
                System.exit(1);
            } else {
                return defaultString;
            }
        }
        if (index >= args.length) {
            System.err.println("Missing argument value for --" + name);
            System.exit(1);
        }
        return args[index];
    }

}
