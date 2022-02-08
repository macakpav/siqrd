#ifndef LINESEARCH_HPP
#define LINESEARCH_HPP
/*
    Approximate line search to determine optimal step size using Wolfe's conditions.
*/

#include <cassert>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
namespace ublas = boost::numeric::ublas;

namespace optimization
{
    template <typename vector_type>
    class LineSearch
    {
    public:
        typedef typename std::enable_if<std::is_floating_point<typename vector_type::value_type>::value,
                                        typename vector_type::value_type>::type value_type;
        typedef typename std::enable_if<std::is_integral<typename vector_type::size_type>::value,
                                        typename vector_type::size_type>::type size_type;

    private:
        const value_type min_step_;
        vector_type gradient_, p_step_;

        // constants for line search
        static const value_type constexpr C1 = 1e-4,
                                          C2 = 0.9;

    public:
        LineSearch(size_type dim, value_type tolerance) : min_step_(tolerance), gradient_(dim), p_step_(dim){};
        ~LineSearch(){};

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
        template <typename v1, typename v2, typename v3, typename s1, typename s2, typename target_functor>
        typename std::enable_if<std::is_floating_point<s1>::value && std::is_floating_point<s2>::value,
                                value_type>::type
        operator()(const v1 &pk, const v2 &dk, const s1 target_k, const v3 &grad_target_k, target_functor &target, s2 step_size)
        {
            assert(pk.size() == dk.size());
            assert(pk.size() == grad_target_k.size());

            const auto dk_grad_prod = ublas::inner_prod(dk, grad_target_k);

            // right hand sides of Wolfe's conditions
            auto rhs1 = [target_k, dk_grad_prod](auto const step_size) { return target_k + C1 * step_size * dk_grad_prod; };
            const auto rhs2 = -C2 * dk_grad_prod;

            p_step_ = pk + step_size * dk;
            auto target_after_step = target(p_step_);

            size_t i;
            for (i = 0; i < 100 && step_size > min_step_; i++)
            {
                target.gradient(p_step_, target_after_step, gradient_);
                // check Wolfe's conditions
                if ((target_after_step <= rhs1(step_size)) &&
                    ((-1.0) * ublas::inner_prod(dk, gradient_) <= rhs2))
                {
                    break;
                }
                step_size = step_size / 2;
                p_step_.assign(pk + step_size * dk);
                target_after_step = target(p_step_);
            }
#ifdef DLVL1
            std::cout << "\tChosen step size in " << i << " iterations: " << step_size << std::endl;
#endif
            return step_size;
        }
    };

} // namespace optimization

#endif