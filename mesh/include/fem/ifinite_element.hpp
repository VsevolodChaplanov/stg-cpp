#ifndef STG_IFINITE_ELEMENT_HPP
#define STG_IFINITE_ELEMENT_HPP

#include <array>
#include <complex>
#include <concepts>
#include <cstddef>
#include <functional>
#include <geometry/geometry.hpp>
#include <numeric>
#include <range/v3/range_concepts.hpp>
#include <ranges>
#include <vector>

namespace stg::mesh {

    template<std::floating_point T>
    class IFiniteElement {
    public:
        using value_type = T;

        virtual value_type lumped(std::size_t i) const = 0;
        virtual value_type jacobian(std::size_t i, std::size_t j) const = 0;
        virtual value_type inverse_j(std::size_t i, std::size_t j) const = 0;
        virtual constexpr std::size_t basis_functions_n() const noexcept = 0;
        virtual const std::vector<std::size_t>& global_indices() const noexcept = 0;
        virtual Point<value_type> real_space_point(Point<T> param_point) const noexcept = 0;
        /* Вытащить отсюда в логику аппроксиматора */
        virtual value_type integrate(const std::vector<T>& values) const = 0;
        virtual std::complex<value_type> integrate(const std::vector<std::complex<T>>& values) const = 0;
        virtual value_type interpolate(const Point<T>& param_point, const std::vector<value_type>& values) const = 0;
        virtual ~IFiniteElement() = default;
    };

    template<std::floating_point T, std::size_t np>
    class LagrangianFiniteElement : public IFiniteElement<T> {
    public:
        using value_type = IFiniteElement<T>::value_type;
        using IFiniteElement<T>::inverse_j;

        LagrangianFiniteElement(std::vector<size_t> global_indices,
                                std::array<value_type, np> lumped,
                                Point<T> base_vertex) noexcept
            : global_indices_{std::move(global_indices)},
              lumped_{std::move(lumped)},
              base_point_{std::move(base_vertex)} {}

        /*
        * Get lumped mass matrix element
        */
        value_type lumped(size_t i) const override { return lumped_[i]; }

        value_type integrate(const std::vector<value_type>& values) const override {
            return integrate(values.begin(), values.end());
        }

        std::complex<value_type> integrate(const std::vector<std::complex<T>>& values) const override {
            std::complex<value_type> result{0, 0};
            for (const std::size_t index: std::views::iota(0ul, np)) {
                result += lumped_[index] * values[index];
            }
            return result;
        }

        /*
        * Integrate values in finite element
        * The values in the passed array must be
        * in the same sequence as the vertices inside the element
        */
        template<std::forward_iterator Iter>
        value_type integrate(Iter begin, Iter end) const {
            const size_t values_amount = std::distance(begin, end);
            if (values_amount != np) {
                throw std::logic_error("Dimensions are not same");
            }
            const auto result = std::inner_product(lumped_.cbegin(), lumped_.cend(), begin, 0.);
            return result;
        }

        /*
        * Returns function gradients of values under iterator in any point
        * of parametric space in finite element
        */
        template<std::forward_iterator Iter>
        Vector<T> gradient_at_point(Iter begin, Iter end, Point<T> param_point) const {
            Vector<T> result{0., 0., 0.};
            const auto gradients = basis_gradients_at_point(param_point);

            for (const auto& gradient: gradients) {
                const auto val = *begin;
                const T x_elem = gradient.template get<0>() * val * inverse_j(0, 0);
                const T y_elem = gradient.template get<1>() * val * inverse_j(1, 1);
                const T z_elem = gradient.template get<2>() * val * inverse_j(2, 2);
                result += Vector<T>{x_elem, y_elem, z_elem};
                ++begin;
            }

            return result;
        }

        /*
        * Return gradients of function values under iterator
        * Assumes that function values are calculated in vertices
        */
        template<std::forward_iterator Iter>
        std::array<Vector<T>, 8> gradients(Iter begin, Iter end) const {
            return {
                    gradient_at_point(begin, end, {0, 0, 0}),
                    gradient_at_point(begin, end, {1, 0, 0}),
                    gradient_at_point(begin, end, {1, 1, 0}),
                    gradient_at_point(begin, end, {0, 1, 0}),
                    gradient_at_point(begin, end, {0, 0, 1}),
                    gradient_at_point(begin, end, {1, 0, 1}),
                    gradient_at_point(begin, end, {1, 1, 1}),
                    gradient_at_point(begin, end, {0, 1, 1}),
            };
        }

        /*
        * Get number of basis functions in finite element
        */
        constexpr size_t basis_functions_n() const noexcept override { return np; }

        /*
        * Return global indices of finite element
        */
        const std::vector<std::size_t>& global_indices() const noexcept override { return global_indices_; }

    protected:
        const std::vector<size_t> global_indices_;
        const std::array<value_type, np> lumped_;
        const Point<value_type> base_point_;

        virtual std::array<Vector<T>, np> basis_gradients_at_point(Point<T> point) const = 0;
    };


}// namespace stg::mesh

#endif//STG_IFINITE_ELEMENT_HPP
