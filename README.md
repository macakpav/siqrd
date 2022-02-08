## SIQRD
A pandemic prediction model of 5 ODEs and 5 coefficients/parameters - infection rate (beta), rate of immunity loss (mu), recovery rate (gamma), rate of dicovery and isolating infected people (delta) and death rate (alpha).
Solves initial value problem of SIQRD equations as they were shown during course of Scientific Software on KU Leuven in 2020/21. Implements 3 ddt schemes for intial value problem of ODEs - forward Euler, backward Euler and Heun's method.
Implemented using UBLAS library from BOOST.

### Core code concepts

#### Debugging 
Assert statements and debug levels are used to ease debugging of code. Compilation differs for debugging purposes and final, optimized code (-DNDEBUG flag) is not slowed down by these additional bits. 

#### Modularity
Is kept in mind to reduce run-time overhead but also to increase reusability of code.

#### Generics over Inheritance
The code prefers Generic programming (templating, concepts) over inheritance to maximize efficiency in terms of both memory and speed.

### Executables
Compilable executable files '.cpp' are located in 'cpp/src' folder.

#### Solvertest
Tests the ODE solvers on special system of ODEs for which analytical solution is known
dot{x}_n(t) = − 10 (x_n − (n-1)/10.0)^3 for n in 1:50
Initial condition [0.01 0.02 0.03 ... 0.5]

#### Simulation
Simulates SIQRD equations with all three ddt methods. Demonstrates the effect of delta coefficient on the results.
Variation of delta parameter (3 values coresponding to different strength of counter-measures) and output file names are hardcoded. Initial condition of infected (I0) and susceptible (S0) people is read together with other model parameters from 'inputs/parameters.in'. 

#### Estimation
Optimizes parameter values of SIQRD equations against an arbitrary target function - e.g. least square error (LSE) of the simulation and experimental data. Uses gradient optimisation methods: Conjugate gradient method (CGM) and 
Broyden–Fletcher–Goldfarb–Shanno (BFGS) using line search with Wolfe conditions. For LSE, two example input setups are present in input directory 'observations1.in' and 'observations2.in', each with their corresponding initial estimate of parameters in a different file 'parameters_observations*.in'. The data of second observation have a certain perturbation added to values of S, I and R counts, in contrast to first observation set with smooth data.

##### Estimation1
Uses Heun's scheme, to optimize parameters against both input observations. Uses both CGM and BFGS with tolerance 1e-12, which is checked against method-dependent residual.

##### Estimation2
Uses all three schemes, to optimize parameters against both input observations, but only with BFGS using tolerance 1e-7.

#### Benchmarking
Files 'bench_mem.cpp' and 'bench_time.cpp' are used to check memory issues and time efficiency (speed) respectively. Runnable using provided Makefile from 'cpp/' folder using make mem and make time.

## Usage
Allrun and Allclean scripts.


