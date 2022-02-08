/*
    Name:     simulation
    Purpose:  Simulates SIQRD equations with all methods. Variation of delta parameter and output file names are hardcoded.
    Author:   pavel.macak@fs.cvut.cz
    Compilation: Makefile is provided - make run1 to run after compilation, make simulation to only compile.
    Command line arguments: (2) Number of time steps and final simulation time.
    Input files: 'parameters.in'
    Output files: Yes
*/

#include "debug_levels.hpp"

#include <cassert>
#include <iostream>
#include <iomanip>
#include <string>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/assignment.hpp>

namespace ublas = boost::numeric::ublas;

#include "siqrd/odeSys_siqrd.hpp"
#include "siqrd/outputFileName.hpp"
#include "ode/odeSolver.hpp"
#include "ode/eulerForward.hpp"
#include "ode/heun.hpp"
#include "ode/eulerBackward.hpp"
#include "saving/saveResults.hpp"

int main(int argc, char const *argv[])
{
    typedef float working_precision;

#ifndef NINFO
    std::cout << "Program started." << std::endl;
#endif
    assert(argc == 3);

#ifdef DLVL0
    std::cout << "Command line arguments: " << std::endl;
    for (int i = 1; i < argc; i++)
    {
        std::cout << argv[i] << std::endl;
    }
    std::cout << std::endl;
#endif
    int N = atoi(argv[1]);
    working_precision T = atof(argv[2]);

    assert(N > 0);
    assert(T > 0);

    siqrd::OdeSys_SIQRD<working_precision> eqns("inputs/parameters.in");
    typedef typename ode::EulerForward<decltype(eqns)> fwe;
    typedef typename ode::EulerBackward<decltype(eqns)> bwe;
    typedef typename ode::Heun<decltype(eqns)> heun;

    ublas::matrix<working_precision, ublas::column_major> scratch_space(decltype(eqns)::dim, N + 1);
    auto parameters = eqns.parameters();

    ode::OdeSolver<fwe> fwe_solver(N, T);
    ode::OdeSolver<bwe> bwe_solver(N, T);
    ode::OdeSolver<heun> heun_solver(N, T);

    parameters[3] = 0.0;
    eqns.set_parameters(parameters);
    fwe_solver.solve(eqns, scratch_space);
#ifdef DLVL1
    std::cout << "Last values: " << std::endl
              << "Suspicable:  " << scratch_space(0, N) << std::endl
              << "Infected:    " << scratch_space(1, N) << std::endl
              << "Quarantined: " << scratch_space(2, N) << std::endl
              << "Recovered:   " << scratch_space(3, N) << std::endl
              << "Dead:        " << scratch_space(4, N) << std::endl
              << std::endl;
#endif
    saving::saveResults(T / N, scratch_space, "outputs/fwe_no_measures.out");

    parameters[3] = 0.2;
    eqns.set_parameters(parameters);

    // this was a problematic case, probably due to permuation matrix bugfix it now works
    // parameters <<= -0.944722, 1.69566, -1.49697, -1.50073, 0.0986043;
    // eqns.set_parameters(parameters);
    // parameters <<= 15000, 5, 0, 0, 0;
    // eqns.set_initial_condition(parameters);

    bwe_solver.solve(eqns, scratch_space);
#ifdef DLVL1
    std::cout << "Last values: " << std::endl
              << "Suspicable:  " << scratch_space(0, N) << std::endl
              << "Infected:    " << scratch_space(1, N) << std::endl
              << "Quarantined: " << scratch_space(2, N) << std::endl
              << "Recovered:   " << scratch_space(3, N) << std::endl
              << "Dead:        " << scratch_space(4, N) << std::endl
              << std::endl;
#endif
    saving::saveResults(T / N, scratch_space, "outputs/bwe_quarantine.out");

    parameters[3] = 0.9;
    eqns.set_parameters(parameters);
    heun_solver.solve(eqns, scratch_space);
#ifdef DLVL1
    std::cout << "Last values: " << std::endl
              << "Suspicable:  " << scratch_space(0, N) << std::endl
              << "Infected:    " << scratch_space(1, N) << std::endl
              << "Quarantined: " << scratch_space(2, N) << std::endl
              << "Recovered:   " << scratch_space(3, N) << std::endl
              << "Dead:        " << scratch_space(4, N) << std::endl
              << std::endl;
#endif
    saving::saveResults(T / N, scratch_space, "outputs/heun_lockdown.out");

#ifndef NINFO
    std::cout << "Program finished." << std::endl;
#endif
    return 0;
}
