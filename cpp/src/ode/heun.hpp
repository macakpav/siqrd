#ifndef HEUN_HPP
#define HEUN_HPP
/*
    Heun method for solving ODE system.
*/

#include <cassert>
#include <iostream>

#include <boost/numeric/ublas/vector.hpp>
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
    class Heun
    {
    public:
        typedef typename std::enable_if<std::is_floating_point<typename OdeSystem::value_type>::value,
                                        typename OdeSystem::value_type>::type value_type;
        typedef typename std::enable_if<std::is_integral<typename OdeSystem::size_type>::value,
                                        typename OdeSystem::size_type>::type size_type;

    private:
        value_type dT_;
        ublas::vector<value_type> temp_;

    public:
        static const char constexpr method_name[] = "heun";
        static const size_type constexpr dim = OdeSystem::dim;

    public:
        Heun(){};
        Heun(const size_type steps, const value_type final_time)
            : dT_(final_time / (value_type)steps), temp_(dim){};
        ~Heun(){};

    public:
        template <typename v1, typename v2>
        inline void time_step(const OdeSystem &system, const v1 &old_time, v2 &new_time)
        {
            assert(old_time.size() == dim);
            assert(old_time.size() == new_time.size());
#ifdef DMETHODS
            std::cout << "Old time: " << old_time << std::endl;
#endif

            system(old_time, temp_);
            new_time.assign(old_time + dT_ * (0.5 * temp_ + 0.5 * system(old_time + dT_ * temp_)));

#ifdef DMETHODS
            std::cout << "New time: " << new_time << std::endl;
#endif
        }
    };
} // namespace ode
#endif
