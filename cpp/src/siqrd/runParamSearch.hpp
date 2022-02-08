#ifndef RUNPARAMSEARCH_HPP
#define RUNPARAMSEARCH_HPP
/*
    Wrapper functions to running CGM and BFGS methods.
*/

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
namespace ublas = boost::numeric::ublas;

#include "../saving/saveResults.hpp"
#include "lse_siqrd.hpp"
#include "../optimization/cgm.hpp"
#include "../optimization/bfgs.hpp"

namespace siqrd
{
    template <typename scheme, typename nu_k_formula = optimization::FR_formula>
    void runCGM(std::string observations, std::string parameters, typename scheme::value_type tol)
    {
        const std::string in_folder = "inputs/",
                          out_folder = "outputs/",
                          out_file = out_folder + scheme::method_name + "_cgm_" + observations + ".out",
                          observ_file = in_folder + observations + ".in",
                          param_file = in_folder + parameters + ".in";

        typedef typename scheme::value_type working_precision;

        siqrd::LSE_siqrd<scheme>
            target_evaluator(observ_file, param_file);
        // get information about SIQRD eqns from LSE object
        auto eqns = target_evaluator.get_eqns(); //copy constructor? (should be correct, only float-type members)
        const int N = target_evaluator.get_N();
        const working_precision T = target_evaluator.get_T();
        const auto starting_parameters = eqns.parameters();
        ublas::matrix<working_precision, ublas::column_major> scratchSpace(decltype(eqns)::dim, N + 1);

        // run the search, simulate again, write results
        eqns.set_parameters(starting_parameters);
        auto final_params = optimization::CGM<nu_k_formula>(target_evaluator, starting_parameters, tol);
        eqns.set_parameters(final_params);
        ode::OdeSolver<scheme> solver(N, T);
        solver.solve(eqns, scratchSpace);
        // std::cout << scheme::method_name << " found with CGM parameters: " << final_params << std::endl;
        saving::saveResults(T / N, scratchSpace, out_file);
    }

    template <typename scheme>
    void runBFGS(std::string observations, std::string parameters, typename scheme::value_type tol)
    {
        const std::string in_folder = "inputs/",
                          out_folder = "outputs/",
                          out_file = out_folder + scheme::method_name + "_bfgs_" + observations + ".out",
                          observ_file = in_folder + observations + ".in",
                          param_file = in_folder + parameters + ".in";

        typedef typename scheme::value_type working_precision;
        ublas::matrix<working_precision, ublas::column_major> identity_matrix = ublas::identity_matrix<working_precision>(5);

        siqrd::LSE_siqrd<scheme>
            target_evaluator(observ_file, param_file);
        // get information about SIQRD eqns from LSE object
        auto eqns = target_evaluator.get_eqns(); //copy constructor? (should be correct, only float-type members)
        const int N = target_evaluator.get_N();
        const working_precision T = target_evaluator.get_T();
        const auto starting_parameters = eqns.parameters();
        ublas::matrix<working_precision, ublas::column_major> scratchSpace(decltype(eqns)::dim, N + 1);

        // run the search, simulate again, write results
        eqns.set_parameters(starting_parameters);
        auto final_params = optimization::BFGS(target_evaluator, starting_parameters, tol);
        eqns.set_parameters(final_params);
        ode::OdeSolver<scheme> solver(N, T);
        solver.solve(eqns, scratchSpace);
        // std::cout << scheme::method_name << " found with BFGS parameters: " << final_params << std::endl;
        saving::saveResults(T / N, scratchSpace, out_file);
    }

} // namespace siqrd

#endif