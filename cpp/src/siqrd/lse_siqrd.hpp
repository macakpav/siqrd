#ifndef LSE_SIQRD_HPP
#define LSE_SIQRD_HPP
/*
    Least square error calculation of SIQRD equations. 
    Also LSE gradient approximation using finite difference. 
*/

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
namespace ublas = boost::numeric::ublas;

#include "odeSys_siqrd.hpp"
#include "../ode/odeSolver.hpp"

namespace siqrd
{

    /*
    pases SchemeType concept to OdeSolver template

////Satisfies concepts:
target_functor 
    member types:
        size_type, value_type
    member functions:
        value_type target_fun(const& variables_vector)
        void gradient(const& variables_vector, const& value_type, & gradient_vect_out)
    static variables:
        size_type dim
    */
    template <typename SchemeType>
    class LSE_siqrd
    {
    public:
        typedef typename std::enable_if<std::is_floating_point<typename SchemeType::value_type>::value,
                                        typename SchemeType::value_type>::type value_type;
        typedef typename std::enable_if<std::is_integral<typename SchemeType::size_type>::value,
                                        typename SchemeType::size_type>::type size_type;
        typedef SchemeType method;

    private:
        size_type no_days_;
        ublas::matrix<value_type, ublas::column_major> prediction_;
        ublas::matrix<value_type, ublas::column_major> scratch_space_;
        ublas::vector<value_type> params_temp_, init_cond_;
        value_type pop_size_squared_;

        static const size_type constexpr RATIO = 8;
        static const value_type constexpr EPS = 1e-5; // step value for finite difference

        OdeSys_SIQRD<value_type, size_type> eqns_;
        ode::OdeSolver<SchemeType> solver_;

    private:
        const static size_type constexpr eqns_dim = OdeSys_SIQRD<>::dim;

    public:
        const static size_type constexpr dim = OdeSys_SIQRD<>::no_params;

    public:
        LSE_siqrd(const std::string &observation_file, const std::string &parameter_file) : params_temp_(dim)
        {
            std::ifstream file(observation_file);
            size_type file_dim;
            file >> no_days_ >> file_dim;
            assert(file_dim == eqns_dim);
            prediction_ = ublas::matrix<value_type, ublas::column_major>(eqns_dim, no_days_);
            scratch_space_ = ublas::matrix<value_type, ublas::column_major>(eqns_dim, (no_days_ - 1) * RATIO + 1);

            eqns_ = decltype(eqns_)(parameter_file, false);
            solver_ = decltype(solver_)(scratch_space_.size2() - 1, (value_type)(no_days_ - 1));

            value_type unused;
            for (size_type i = 0; i < no_days_; i++)
            {
                file >> unused;
                for (size_type j = 0; j < eqns_dim; j++)
                {
                    file >> prediction_(j, i);
                }
            }
            // std::cout << prediction_ << std::endl;
            file.close();
            init_cond_ = ublas::column(prediction_, 0);
            eqns_.set_initial_condition(init_cond_);
            pop_size_squared_ = std::accumulate(init_cond_.begin(), init_cond_.end(), 0.0);
            pop_size_squared_ *= pop_size_squared_;
        };
        ~LSE_siqrd(){};

    public:
        auto get_eqns() { return eqns_; }
        auto get_N() { return no_days_ * RATIO; }
        auto get_T() { return no_days_; }

    public:
        template <typename vect>
        inline value_type operator()(vect const &p)
        {
            assert(p.size() == dim);
            return lse(p);
        }

    private:
        template <typename vect>
        value_type lse(vect const &params)
        {
            assert(params.size() == dim);
            assert((scratch_space_.size2() - 1) / (no_days_ - 1) == RATIO);

            eqns_.set_initial_condition(init_cond_);
            eqns_.set_parameters(params);
            solver_.solve(eqns_, scratch_space_);

            value_type lse = 0.0;
            size_type i = 0;
            for (size_type day = 0; day < no_days_; day++)
            {
                auto x_ip = ublas::column(scratch_space_, i);
                auto x_i = ublas::column(prediction_, day);
                lse += pow(ublas::norm_2(x_i - x_ip), 2);

                i += RATIO;
            }

            lse /= ((value_type)(no_days_)*pop_size_squared_);
#ifdef DLVL3
            std::cout << "LSE: " << lse << std::endl
                      << std::endl;
#endif
            return lse;
        };

    public:
        template <typename v1, typename v2>
        void gradient(v1 const &p, const value_type lse_0, v2 &grad)
        {
            assert(p.size() == dim);
            assert(grad.size() == dim);

            params_temp_.assign(p);
            params_temp_[0] += EPS;
            grad[0] = (lse(params_temp_) - lse_0) / EPS;
            for (decltype(p.size()) i = 1; i < p.size(); i++)
            {
                params_temp_[i-1]=p[i-1];
                params_temp_[i] += EPS;
                grad[i] = (lse(params_temp_) - lse_0) / EPS;
            }
#ifdef DLVL3
            std::cout << "gradient of LSE: " << std::endl
                      << grad << std::endl;
#endif
        }
    };
} // namespace siqrd

#endif