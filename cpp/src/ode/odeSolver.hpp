#ifndef ODESOLVERS_HPP
#define ODESOLVERS_HPP
/*
    OdeSolver class that uses method to solve a system of ODE until target time T with N steps.
*/

#include <cassert>
#include <iostream>
#include <algorithm>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

namespace ublas = boost::numeric::ublas;
namespace ode
{
    /*
////Satisfies concepts:
OdeSolver
    constructor:
        OdeSolver(const int noSteps, const value_type maxTime, OdeSystem &ode_sys)
    member types:
        size_type, value_type
    member functions:
        void solve(matrix_type &results_matrix)
    static variables:
        size_type dim


////Uses concepts:
SchemeType 
    constructor:
        SchemeType(const value_type steps, const value_type final_time)
    member types:
        size_type, value_type
    member functions:
        void time_step(const OdeSystem &ode_sys, const vector_type &old_time, vector_type &new_time)
    static variables:
        size_type dim
        char[] method_name

OdeSystem 
    member types:
        size_type, value_type
    member functions:
        vector_type initial_condition() const
        vector_type operator()( vector_type )
        return_vector operator()(const v1 &variables_vector)
        void operator()(const v1 &variables_vector, v2 &return_vector)
        void jacobian()( variables, &output_matrix )
    static variables:
        size_type dim
*/

    template <typename SchemeType>
    class OdeSolver
    {
    public:
        typedef typename std::enable_if<std::is_floating_point<typename SchemeType::value_type>::value,
                                        typename SchemeType::value_type>::type value_type;
        typedef typename std::enable_if<std::is_integral<typename SchemeType::size_type>::value,
                                        typename SchemeType::size_type>::type size_type;

    private:
        size_type N_;
        value_type T_;
        SchemeType method_;

    public:
        static const size_type constexpr dim = SchemeType::dim;

    public:
        //default constructor
        OdeSolver(){};
        //basic constructor
        OdeSolver(const int noSteps, const value_type maxTime)
            : N_(noSteps), T_(maxTime), method_(N_, T_){};

        // solve ode system using method_, put results to results_matrix, first column is assigned initial condition
        template <typename OdeSystem, typename matrix_type> //should be column major
        void solve(OdeSystem &ode_sys, matrix_type &results_matrix)
        {
#ifdef DODESOLVER
            std::cout << "Solving ODE ode_sys using " << SchemeType::method_name << std::endl;
#endif
            // std::cout << "Eqns parameters: " << ode_sys.parameters() << std::endl;
            assert(OdeSystem::dim == results_matrix.size1());
            assert(results_matrix.size2() == N_ + 1);
#ifndef NDEBUG
            for (decltype(results_matrix.size1()) i = 0; i < results_matrix.size1(); i++)
                for (decltype(results_matrix.size2()) j = 0; j < results_matrix.size2(); j++)
                    results_matrix(i, j) = std::numeric_limits<value_type>::quiet_NaN();
#endif

            // assign inital condition to first column
            auto init = ublas::column(results_matrix, 0);
            init = ode_sys.initial_condition();
#ifdef DODESOLVER
            std::cout << "Initial condition: " << std::endl
                      << init << std::endl;
#endif
            for (size_type step = 0; step < N_; step++)
            {
                // two columns of matricies
                auto old_time = ublas::column(results_matrix, step);
                auto new_time = ublas::column(results_matrix, step + 1);
                // compute new_time using method
                method_.time_step(ode_sys, old_time, new_time);

#ifdef DODESOLVER
                if (step % (N_ / 10) == N_ / 10 - 1)
                    std::cout << "Done " << step + 1 << " steps." << std::endl
                              << std::endl;
#endif
            }
#ifdef DODESOLVER
            std::cout << "Last values: " << std::endl
                      << "First variable:  " << results_matrix(0, N_) << std::endl
                      << "Last variable:   " << results_matrix(OdeSystem::dim - 1, N_)
                      << std::endl;
#endif
        };
    };
} // namespace ode
#endif
