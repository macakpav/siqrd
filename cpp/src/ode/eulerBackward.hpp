#ifndef EULERBACKWARD_HPP
#define EULERBACKWARD_HPP
/*
    Euler backward method for solving ODE system.
*/

#include <cassert>
#include <iostream>
#include <numeric>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
namespace ublas = boost::numeric::ublas;

namespace ode
{
    /*
    Satisfies concepts:
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


    Uses concepts:
OdeSystem 
    member types:
        size_type, value_type
    member functions:
        vector_type initial_condition() const
        vector_type operator()( vector_type )
        void operator()(const v1 &variables_vector, v2 &return_vector)
        void jacobian()( variables, &output_matrix )
    static variables:
        size_type dim
*/

    template <typename OdeSystem>
    class EulerBackward
    {
    public:
        typedef typename std::enable_if<std::is_floating_point<typename OdeSystem::value_type>::value,
                                        typename OdeSystem::value_type>::type value_type;
        typedef typename std::enable_if<std::is_integral<typename OdeSystem::size_type>::value,
                                        typename OdeSystem::size_type>::type size_type;

    private:
        value_type dT_;
        ublas::vector<value_type> temp_;
        ublas::matrix<value_type> jac_;
        ublas::permutation_matrix<int> pm_, pm_default;

        // method private settings
        static const value_type constexpr tolerance = std::numeric_limits<value_type>::epsilon() * 100.0;
        static const size_type constexpr max_iter = 1000;

    public:
        static const char constexpr method_name[] = "bwe";
        static const value_type constexpr dim = OdeSystem::dim;

    public:
        EulerBackward() : pm_(0), pm_default(0){};
        EulerBackward(const value_type steps, const value_type final_time)
            : dT_(final_time / steps), temp_(dim), jac_(dim, dim), pm_(dim), pm_default(dim){};
        ~EulerBackward(){};

    public:
        template <typename v1, typename v2>
        void time_step(const OdeSystem &system, const v1 &old_time, v2 &new_time)
        {
            assert((decltype(dim))old_time.size() == dim);
            assert((decltype(new_time.size()))old_time.size() == new_time.size());
            value_type res = std::numeric_limits<value_type>::infinity();
            value_type norm = ublas::norm_1(old_time);
#ifdef DMETHODS
            std::cout << "Old time: " << old_time << std::endl;
            std::cout << "Normalizer: " << norm << std::endl;
            bool converged = false;
#endif
            new_time.assign(old_time);
            size_type i;
            for (i = 0; i < max_iter; i++)
            {
                // check convergence
                temp_.assign((old_time - new_time) + dT_ * system(new_time));
                res = ublas::norm_inf(temp_) / norm;
                if (res < tolerance)
                {
                    #ifdef DMETHODS
                        bool converged = true;
                    #endif
                    break;
                }

                // jacobian matrix
                jacG(new_time, jac_, system);

                // compute inversion
                pm_.assign(pm_default);
                ublas::lu_factorize(jac_, pm_);
                ublas::lu_substitute(jac_, pm_, temp_);

                // new iteration
                new_time.assign(new_time - temp_);
            }

#ifdef DMETHODS
            if (converged)
            {

                std::cout << "Euler backward converged in " << i + 1 << " steps." << std::endl
                          << std::endl;
            }
            else
            {
                std::cout << "Euler backward did not converge! " << std::endl;
                std::exit(1);
            }
#endif

#ifdef DMETHODS
            std::cout << "New time: " << new_time << std::endl;
            std::cout << "Residual: " << res << std::endl;
#endif
        }

    private:
        template <typename vector, typename matrix_type>
        inline void jacG(const vector &vars, matrix_type &matrix, const OdeSystem &system) const
        {
            system.jacobian(vars, matrix);
            matrix *= dT_;
            for (size_t i = 0; i < matrix.size1(); i++)
            {
                matrix(i, i) -= 1.0;
            }
        }
    };
} // namespace ode
#endif
