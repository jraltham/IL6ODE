
  # IL6ODE

**Ordinary Differential Equation Modeling of the IL-6–Tumor Microenvironment–Immune Checkpoint Inhibitor Axis**

## Overview

This repository contains the MATLAB code developed by James Altham to investigate the IL-6–tumor microenvironment (TME)–immune checkpoint inhibitor (ICI) axis. The framework uses a system of ordinary differential equations (ODEs) to investigate treatment dynamics, optimize anti-IL-6 dosing strategies, characterize tumor and toxicity outcomes, and evaluate parameter sensitivity.

The codebase includes the ODE model and simulation framework, parameter definitions, global sensitivity analysis (GSA), virtual population simulations, toxicity parameter calibration, and manuscript figure generation.

## Software Requirements

The code was developed and tested using **MATLAB R2026b (26.2.0.3386108)**.

The following MATLAB toolboxes were used in the project:

* Global Optimization Toolbox
* Image Processing Toolbox
* Optimization Toolbox
* Parallel Computing Toolbox
* Partial Differential Equation Toolbox
* Signal Processing Toolbox
* Statistics and Machine Learning Toolbox

Specific toolbox requirements may vary by script. The Global Optimization Toolbox is required for the particle swarm optimization used in toxicity parameter calibration.

## File Descriptions

### Core ODE Model and Simulation

* **`IL6DosingODE.m`**: Defines the ODE system, including the core model equations and supplementary equations governing the modeled biological and treatment dynamics.
* **`IL6DosingParameters.m`**: Defines model parameter values and, where applicable, parameter ranges used for simulation, calibration, and sensitivity analysis.
* **`IL6ParamUnpack.m`**: Extracts and organizes parameters defined in `IL6DosingParameters.m` into variables accessible to the model and associated analysis functions.
* **`IL6RunSim.m`**: Executes the ODE simulation using specified initial conditions, simulation time spans, and model parameters. The system is solved using MATLAB's `ode15s` solver to accommodate potentially stiff dynamics arising from processes operating on different timescales.

### Global Sensitivity Analysis and Figure Generation

* **`IL6DosingSobolGSA.m`**: Implements global sensitivity analysis using Saltelli sampling to evaluate the influence of model parameters on selected outcomes. This script also generates the complete Figure 3 presented in the manuscript.
* **`IL6DosingFigures.m`**: Generates the core manuscript figures, including Figures 2, 4, 5, and 6. Figure 3 is generated separately by `IL6DosingSobolGSA.m`.

### Virtual Population and Parameter Calibration

* **`populationLoss.m`**: Calculates the calibration loss for the toxicity parameters \(L_{\mathrm{tox}}\), \(E_{\mathrm{tox}}\), \(r_x\), and \(r_{x,\mathrm{res}}\) by comparing simulated population-level toxicity outcomes with literature-derived calibration targets.
* **`populationSimulation.m`**: Simulates a virtual population of \(N_P\) patients for a given candidate parameter set. Patient-specific parameters are sampled uniformly from their specified ranges, and individual toxicity outcomes are aggregated into population-level statistics.
* **`runPopulationCalibration.m`**: Runs particle swarm optimization to estimate the toxicity parameters by minimizing the population-level calibration loss. It then evaluates the optimized parameters across a virtual patient population and visualizes and saves the results.
* **`toxicityFeatures.m`**: Extracts individual-patient toxicity features from the ODE solution, including normalized toxicity onset time, resolution time, and maximum toxicity.
* **`tumorOutcome.m`**: Extracts a tumor outcome metric from the ODE solution. If a local minimum in tumor burden is identified, the function returns the minimum tumor burden; otherwise, it returns the final tumor burden at the end of the simulation.

## Model and Analysis Overview

The repository supports the following computational workflows:

1. **ODE simulation:** Simulates the coupled biological and treatment dynamics under specified initial conditions, parameter values, and dosing schedules.
2. **Dosing analysis:** Evaluates tumor and toxicity outcomes under different treatment timing and dosing strategies.
3. **Global sensitivity analysis:** Quantifies the influence of model parameters on selected simulation outcomes using Saltelli-based sampling.
4. **Virtual population simulation:** Evaluates model behavior across heterogeneous virtual patients.
5. **Parameter calibration:** Estimates toxicity-related parameters by comparing simulated population outcomes with literature-derived targets.
6. **Figure generation:** Produces the manuscript figures associated with the model simulations and analyses.

## Reproducibility

For reproducible simulations and analyses, use the parameter definitions, initial conditions, dosing schedules, parameter ranges, and random seeds specified in the relevant scripts. Where applicable, record the MATLAB release, toolbox versions, and configuration used to generate each result.
