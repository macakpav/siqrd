#ifndef BFGS_HPP
#define BFGS_HPP
/*
    BFGS algorithm optimizes variables against some target function minimum.
*/

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
namespace ublas = boost::numeric::ublas;

#include "lineSearch.hpp"

namespace optimization
{

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

    template <typename target_functor, typename vector_type, typename scalar_type, typename matrix_type = ublas::matrix<typename target_functor::value_type, ublas::column_major>>
    typename std::enable_if<std::is_floating_point<typename target_functor::value_type>::value &&
                                std::is_integral<typename target_functor::size_type>::value &&
                                std::is_arithmetic<scalar_type>::value,
                            vector_type>::type
    BFGS(target_functor &target_fun,
         const vector_type &starting_variables,
         const scalar_type tolerance,
         matrix_type hessian = ublas::identity_matrix<typename target_functor::value_type>(target_functor::dim)) //copy of hessian matrix is on purpose
    {
#ifdef DLVL1
        std::cout << "Starting BFGS" << std::endl;
#endif
        typedef typename target_functor::value_type value_type;
        typedef typename target_functor::size_type size_type;

        const size_type dim = starting_variables.size();
        const size_type max_iters = 1000;
        const value_type max_step_size = 1.0;

        assert(hessian.size1() == hessian.size2());
        assert(hessian.size1() == dim);

        bool converged = false;
        value_type target_k, step_size, res;
        vector_type direction(dim), variables(dim), grad_target_k(dim), grad_target_old(dim),
            y(dim), hs(dim), sh(dim);
        ublas::permutation_matrix<int> pm(dim);
        const ublas::permutation_matrix<int> pm_default(dim);
        matrix_type temp = hessian;
        LineSearch<vector_type> line_search(dim, tolerance * 100);

#ifndef NDEBUG
        for (size_type i = 0; i < dim; i++)
        {
            direction[i] = variables[i] = grad_target_k[i] = grad_target_old[i] = y[i] = hs[i] = sh[i] = std::numeric_limits<value_type>::quiet_NaN();
        }
        target_k = step_size = res = std::numeric_limits<value_type>::quiet_NaN();
#endif

        direction.clear();
        variables.assign(starting_variables);
        target_k = target_fun(variables);
        target_fun.gradient(variables, target_k, grad_target_k);

        size_type k;
        for (k = 0; k < max_iters; k++)
        {
#ifdef DLVL2
            std::cout << "BFGS, iteration " << k << ", variables: " << variables << std::endl
                      << "gradient: " << grad_target_k << std::endl
                      << "Target value: " << target_k << std::endl;
#endif

            //solve inversion to get new direction
            direction.assign(-grad_target_k);
            temp.assign(hessian);
            pm.assign(pm_default); // pm.assign(pm_default);
            ublas::lu_factorize(temp, pm);
            ublas::lu_substitute(temp, pm, direction);
#ifdef DLVL2
            std::cout << "Direction: " << direction << std::endl;
#endif
            //line search for optimal step size
            step_size = line_search(variables, direction, target_k, grad_target_k, target_fun, max_step_size);
#ifdef DLVL2
            std::cout << "Step size: " << step_size << std::endl;
#endif
            //convergence check
            res = step_size * ublas::norm_2(direction) / ublas::norm_2(variables);
            if (res < tolerance)
            {
                converged = true;
                break;
            }

            // update Hessian matrix
            variables.assign(variables + step_size * direction);
            grad_target_old.assign(grad_target_k);
            target_k = target_fun(variables);
            target_fun.gradient(variables, target_k, grad_target_k);

            y.assign(grad_target_k - grad_target_old); // grad_k+1 - grad_k
            // could not make it work by inserting ublas::prod(direction, hessian) straight into the update expresion
            hs.assign(ublas::prod(hessian, direction));
            sh.assign(ublas::prod(direction, hessian));

            noalias(hessian) += -ublas::outer_prod(hs, sh) / (ublas::inner_prod(sh, direction)) +
                                ublas::outer_prod(y, y) / (ublas::inner_prod(direction, y) * step_size);
            //B_k = B_k - (B_k * s x s * B_k) / (s * Bk * s) + (y x y) / (s * y);
#ifdef DLVL1
            std::cout << "BFGS residual in step_size " << k << ": " << res << std::endl;
#endif
        }
        if (converged)
        {
#ifndef NINFO
            std::cout << "BFGS converged in " << k << " iterations. " << std::endl
                      << "Final variables:" << variables << std::endl;
#endif
        }
        else
        {
            std::cerr << std::endl
                      << "BFGS did NOT converge in iteration limit(" << max_iters << ")!" << std::endl
                      << "Variables are:" << variables << std::endl
                      << std::endl;
        }

        return variables;
    };

} // namespace optimization

#endif