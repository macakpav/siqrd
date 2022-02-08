#ifndef EULERFORWARD_HPP
#define EULERFORWARD_HPP
/*
    Euler forward method for solving ODE system.
*/

#include <cassert>
#include <iostream>

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
    class EulerForward
    {
    public:
        typedef typename std::enable_if<std::is_floating_point<typename OdeSystem::value_type>::value,
                                        typename OdeSystem::value_type>::type value_type;
        typedef typename std::enable_if<std::is_integral<typename OdeSystem::size_type>::value,
                                        typename OdeSystem::size_type>::type size_type;

    private:
        value_type dT_;

    public:
        static const char constexpr method_name[] = "fwe";
        static const size_type constexpr dim = OdeSystem::dim;

    public:
        EulerForward(){};
        EulerForward(const size_type steps, const value_type final_time)
            : dT_(final_time / steps){};
        ~EulerForward(){};

    public:
        template <typename v1, typename v2>
        inline void time_step(const OdeSystem &system, const v1 &old_time, v2 &new_time) const
        {
            assert(old_time.size() == dim);
            assert(old_time.size() == new_time.size());
            #ifdef DMETHODS
            std::cout << "Old time: " << old_time << std::endl;
            #endif

            new_time.assign(old_time + dT_ * system(old_time));
            
            #ifdef DMETHODS
            std::cout << "New time: " << new_time << std::endl;
            #endif
        }
    };
} // namespace ode
#endif
