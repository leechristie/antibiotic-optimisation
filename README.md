# Antibiotic Optimisation

## Antibiotic Model 4.0.6

A stochastic model for antibiotic resistance.

<details><summary>Click to see documentation...</summary>

### Example Usage (Model)

`AntibioticModel` can be used directly, giving results to compare to the MATLAB
reference implementation. The model constructor accepts additional parameters to
specify initial bacterial load and random number generator for the simulation.

    int samples = 1000;
    AntibioticModel model = new AntibioticModel(samples);
    
    int[] solution = new int[] {10, 10, 10, 10, 10, 10, 10, 10, 10, 10};

    double fitness = model.evaluate(solution);

### Example Usage (Problem)

`AntibioticProblem` is used to create an `IntegerProblem` instance for jMetal
algorithms. Many objectives may be specified.

    int maxIndividualDosage = 60;
    int maxConcentraition = 60;
    AntibioticProblem problem = new AntibioticProblem(
            model,
            maxIndividualDosage,
            AntibioticObjective.totalAntibiotic(),
            AntibioticObjective.overdoseAmount(maxConcentraition),
            AntibioticObjective.uncuredProportion());

    Algorithm<List<DoubleSolution>> algorithm = ... // jMetal Algorithm usage

</details>

## Antibiotic Model Experiments

Documentation to follow

## Results Post-Processing

Documentation to follow
