/*
    Name:     estimation2
    Purpose:  Runs BFGS with all methods on first example case of observations.
    Author:   pavel.macak@fs.cvut.cz
    Compilation: Makefile is provided - make run4 to run after compilation, make estimation2 to only compile
    Command line arguments: None
    Input files: 'parameters_observations1.in', 'observations?.in'
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
    typedef double working_precision;

#ifndef NINFO
    std::cout << "Program started." << std::endl
              << std::endl;
#endif

    const working_precision tol = 1e-7;
    const std::string observations1 = "observations1",
                      starting_guess1 = "parameters_" + observations1;

    typedef typename ode::EulerForward<siqrd::OdeSys_SIQRD<working_precision>> fwe;
    typedef typename ode::EulerBackward<siqrd::OdeSys_SIQRD<working_precision>> bwe;
    typedef typename ode::Heun<siqrd::OdeSys_SIQRD<working_precision>> heun;

    siqrd::runBFGS<heun>(observations1, starting_guess1, tol);
    siqrd::runBFGS<fwe>(observations1, starting_guess1, tol);
    siqrd::runBFGS<bwe>(observations1, starting_guess1, tol);


#ifndef NINFO
    std::cout << "Program finished." << std::endl;
#endif
    return 0;
}
