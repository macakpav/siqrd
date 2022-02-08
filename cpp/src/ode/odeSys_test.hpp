#ifndef ODE_SYSTEM_SIQRD_HPP
#define ODE_SYSTEM_SIQRD_HPP
/*
    System of ODE dx/dt_n(t) = − 10 (x_n − (n-1)/10.0)^3 for n in 1:50
*/

#include <cassert>
#include <cmath>
#include <algorithm>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

namespace ublas = boost::numeric::ublas;

namespace ode
{
    /*
////Satisfies concepts:
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
    template <typename Type = double, typename SizeType = typename ublas::vector<Type>::size_type>
    class OdeSys_test
    {
    public:
        typedef typename std::enable_if<std::is_integral<SizeType>::value, SizeType>::type size_type;
        typedef typename std::enable_if<std::is_floating_point<Type>::value, Type>::type value_type;

        static const size_type dim = 50;

        OdeSys_test(){};
        ~OdeSys_test(){};

    public:
        //     Initial condition [0.01 0.02 0.03 ... 0.5]
        auto initial_condition() const
        {
            ublas::vector<value_type> ret(dim);
            int i = 0;
            auto fun = [&i]() { return ++i * 0.01; };
            std::generate(ret.begin(), ret.end(), fun);

            return ret;
        }

    public:
        auto analytic_solution(value_type t) const
        {
            ublas::vector<value_type> ret(dim);
            auto init_cond = initial_condition();
            value_type k, c, m;

            for (size_type i = 0; i < dim; i++)
            {
                k = 0.1 * i;
                m = init_cond[i];
                if (std::fabs(m - k) < std::numeric_limits<value_type>::epsilon())
                {
                    ret[i] = 0.0;
                }
                else
                {
                    c = 1 / std::pow(m - k, 2);
                    if (m > k)
                    {
                        ret[i] = k + std::sqrt(1 / (20 * t + c));
                    }
                    else
                    {
                        ret[i] = k - std::sqrt(1 / (20 * t + c));
                    }
                }
            }
            return ret;
        }

    public:
        template <typename vector_type>
        auto operator()(const vector_type &variables_vector) const
        {
            assert(variables_vector.size() == dim);
            ublas::vector<value_type> ret_vector(dim);
            for (size_type i = 0; i < dim; i++)
            {
                ret_vector[i] = -10 * std::pow((variables_vector[i] - 0.1 * i), 3);
            }

            return ret_vector;
        };

        template <typename v1, typename v2>
        void operator()(const v1 &variables_vector, v2 &return_vector) const
        {
            assert(variables_vector.size() == dim);
            assert(return_vector.size() == dim);
            for (size_type i = 0; i < dim; i++)
            {
                return_vector[i] = -10 * std::pow((variables_vector[i] - 0.1 * i), 3);
            }
        };

        template <typename vector_type, typename matrix_type>
        void jacobian(const vector_type &variables_vector, matrix_type &jac_matrix) const
        {
            assert(variables_vector.size() == dim);
            assert(jac_matrix.size1() == variables_vector.size());
            assert(jac_matrix.size1() == jac_matrix.size2());

            jac_matrix.clear();

            for (size_type i = 0; i < dim; i++)
            {
                jac_matrix(i, i) = -30 * std::pow((variables_vector[i] - 0.1 * i), 2);
            }
        };
    };
} // namespace ode
#endif
