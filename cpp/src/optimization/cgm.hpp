#ifndef CGM_HPP
#define CGM_HPP
/*
    CGM algorithm optimizes variables against some target function minimum.
    Also contains formulas for choosing nu_k parameter.
*/

#include <iostream>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
namespace ublas = boost::numeric::ublas;

#include "lineSearch.hpp"

namespace optimization
{
    /*
////Uses concepts:
SchemeType 
    constructor:
        SchemeType(const value_type steps, const value_type final_time)
    member types:
        size_type, value_type
    member functions:
        void time_step(const OdeSystem &system, const vector_type &old_time, vector_type &new_time)
    static variables:
        size_type dim
        char[] method_name
    */

    struct FR_formula
    {
        template <typename v1, typename v2>
        static inline auto nu_k(const v1 &grad_target_k, const v2 &grad_target_old)
        {
            return pow(ublas::norm_2(grad_target_k) / ublas::norm_2(grad_target_old), 2);
        }
    };
    struct PR_formula
    {
        template <typename v1, typename v2>
        static inline auto nu_k(const v1 &grad_target_k, const v2 &grad_target_old)
        {
            return ublas::inner_prod(grad_target_k, (grad_target_k - grad_target_old)) / pow(ublas::norm_2(grad_target_old), 2);
        }
    };
    /*
////Uses concepts:
target_functor 
    member types:
        size_type, value_type
    member functions:
        value_type target_fun(const& variables_vector)
        void gradient(const& variables_vector, const& value_type, & gradient_vect_out)
    static variables:
        size_type dim
    */

    template <typename formula = FR_formula, typename target_functor, typename vector_type, typename scalar_type>
    typename std::enable_if<std::is_floating_point<typename target_functor::value_type>::value &&
                                std::is_integral<typename target_functor::size_type>::value &&
                                std::is_arithmetic<scalar_type>::value,
                            vector_type>::type
    CGM(target_functor &target_fun, const vector_type &starting_variables, const scalar_type tolerance)
    {
#ifdef DLVL1
        std::cout << "Starting CGM" << std::endl;
#endif
        typedef typename target_functor::value_type value_type;
        typedef typename target_functor::size_type size_type;

        const size_type dim = starting_variables.size();
        const size_type max_iters = 1000;
        const value_type max_step_size = 0.01;

        bool converged = false;
        value_type target_k, step_size, res, nu_k, p_norm2;
        vector_type variables(dim), direction(dim), grad_target_k(dim), grad_target_old(dim);
        LineSearch<vector_type> line_search(dim, tolerance);

#ifndef NDEBUG
        for (size_type i = 0; i < dim; i++)
        {
            variables[i] = direction[i] = grad_target_k[i] = grad_target_old[i] = std::numeric_limits<value_type>::quiet_NaN();
        }
        target_k = step_size = res = nu_k = p_norm2 = std::numeric_limits<value_type>::quiet_NaN();
#endif

        direction.clear();
        nu_k = 0;
        variables.assign(starting_variables);

        size_type k;
        for (k = 0; k < max_iters; k++)
        {
            // norm of parameters in step_size k
            p_norm2 = ublas::norm_2(variables);

            // target values in this step_size
            target_k = target_fun(variables);
            grad_target_old.assign(grad_target_k);
            target_fun.gradient(variables, target_k, grad_target_k);

            // new direction
            if ((k % dim != 0)) // this is a constant pattern of if, should get optimized
            {
                // nu_k selection
                nu_k = formula::nu_k(grad_target_k, grad_target_old);
                direction.assign((-1) * grad_target_k + nu_k * direction);
            }
            else
            {
                direction.assign((-1) * grad_target_k);
            }

            // line search for step_size size
            step_size = line_search(variables, direction, target_k, grad_target_k, target_fun, max_step_size);

            // check for convergence
            res = step_size * ublas::norm_2(direction) / p_norm2;
            if ((res < tolerance) && (k % dim != 0))
            {
                // if seemingly converged, try again without nu_k orhogonalization (exactly in gradient direction)
                direction.assign((-1) * grad_target_k);
                step_size = line_search(variables, direction, target_k, grad_target_k, target_fun, max_step_size);
                res = step_size * ublas::norm_2(direction) / p_norm2;
            }

            // real convergence check
            if (res < tolerance)
            {
                converged = true;
                break;
            }
#ifdef DLVL1
            std::cout << "CGM residual in step_size " << k << ": " << res << std::endl;
#endif
            // new values for parameters
            variables.assign(variables + step_size * direction);
#ifdef DLVL2
            std::cout << "CGM new variables " << variables << ": " << res << std::endl
                      << std::endl;
#endif
        }

        if (converged)
        {
#ifndef NINFO
            std::cout << "CGM converged in " << k << " iterations. " << std::endl
                      << "Final variables:" << variables << std::endl;
#endif
        }
        else
        {
            std::cerr << std::endl
                      << "CGM did NOT converge in iteration limit(" << max_iters << ")!" << std::endl
                      << "Variables are:" << variables << std::endl
                      << std::endl;
        }

        return variables;
    }

} // namespace optimization

#endif