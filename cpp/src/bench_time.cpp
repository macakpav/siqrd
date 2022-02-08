/*
    Name:     bench_time
    Purpose:  Timing program, switching between CGM and BFGS hardcoded.
    Author:   pavel.macak@fs.cvut.cz
    Compilation: Makefile is provided - make time to run after compilation, make bench_time to only compile
    Command line arguments: None
    Input files: 'parameters_observations1.in'
    Output files: Yes

    My configuration    (average, standard deviation)
    heun CGM Time(s):   0.23418   0.0111761
    heun BFGS Time(s):  0.0376051 0.00203765
*/

#include "debug_levels.hpp"

#include <iostream>
#include <string>
#include <chrono>

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
    // typedef typename ode::EulerForward<siqrd::OdeSys_SIQRD<working_precision>> fwe;
    // typedef typename ode::EulerBackward<siqrd::OdeSys_SIQRD<working_precision>> bwe;
    typedef typename ode::Heun<siqrd::OdeSys_SIQRD<working_precision>> heun;

#ifndef NINFO
    std::cout << "Program started." << std::endl
              << std::endl;
#endif

    const working_precision tol = 1e-7;
    const std::string observations = "observations1",
                      starting_guess = "parameters_" + observations;

    int number_exp = 100;
    int discard = 5;
    double elapsed_time = 0.;
    double average_time = 0.;
    double squared_time = 0.;
    double time_diff = 0.;
    for (int exp = 0; exp < number_exp + discard; exp++)
    {
        auto t_start = std::chrono::high_resolution_clock::now();



/****************************TIMED BLOCK****************************/


        // siqrd::runCGM<heun>(observations, starting_guess, tol);
        siqrd::runBFGS<heun>(observations, starting_guess, tol);


/*******************************************************************/



        auto t_end = std::chrono::high_resolution_clock::now();
        if (exp >= discard)
        {
            elapsed_time = std::chrono::duration<double>(t_end - t_start).count();
            time_diff = elapsed_time - average_time;
            average_time += time_diff / (exp - discard + 1);
            squared_time += time_diff * time_diff;
        }
    }
    std::cout << "Time(s): " << average_time << " " << std::sqrt(squared_time / (number_exp - 1)) << std::endl;

#ifndef NINFO
    std::cout << "Program finished." << std::endl;
#endif
    return 0;
}
