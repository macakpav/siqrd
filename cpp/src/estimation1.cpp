/*
    Name:     estimation1
    Purpose:  Runs CGM and BFGS with Heun's method on two example cases of observations.
    Author:   pavel.macak@fs.cvut.cz
    Compilation: Makefile is provided - make run3 to run after compilation, make estimation1 to only compile.
    Command line arguments: None
    Input files: 'parameters_observations?.in', 'observations?.in'
    Output files: Yes
*/

#include "debug_levels.hpp"

#include <iostream>
#include <string>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
namespace ublas = boost::numeric::ublas;

#include "ode/eulerForward.hpp"
#include "ode/heun.hpp"
#include "ode/eulerBackward.hpp"

#include "siqrd/runParamSearch.hpp"
#include "siqrd/odeSys_siqrd.hpp"

int main()
{
    typedef long double working_precision;
    const working_precision tol = 1e-12;

#ifndef NINFO
    std::cout << "Program started." << std::endl
              << std::endl;
#endif

    const std::string observations1 = "observations1",
                      starting_guess1 = "parameters_" + observations1,
                      observations2 = "observations2",
                      starting_guess2 = "parameters_" + observations2;

    // way to run using other solvers
    // typedef typename ode::EulerForward<siqrd::OdeSys_SIQRD<working_precision>> fwe;
    // typedef typename ode::EulerBackward<siqrd::OdeSys_SIQRD<working_precision>> bwe;
    // siqrd::runCGM<fwe>(observations, parameters, tol);
    // siqrd::runBFGS<bwe>(observations, parameters, tol);

    typedef typename ode::Heun<siqrd::OdeSys_SIQRD<working_precision>> heun;
    siqrd::runCGM<heun>(observations1, starting_guess1, tol);
    siqrd::runBFGS<heun>(observations1, starting_guess1, tol);
    siqrd::runCGM<heun, optimization::FR_formula>(observations2, starting_guess2, tol);
    siqrd::runBFGS<heun>(observations2, starting_guess2, tol);

#ifndef NINFO
    std::cout << "Program finished." << std::endl;
#endif
    return 0;
}
