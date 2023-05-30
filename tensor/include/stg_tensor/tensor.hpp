#ifndef STG_TENSOR_HPP
#define STG_TENSOR_HPP

#include <array>
#include <bits/ranges_algobase.h>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometry.hpp>
#include <compare>
#include <concepts>
#include <range/v3/range/conversion.hpp>
#include <range/v3/to_container.hpp>
#include <range/v3/view/tokenize.hpp>
#include <range/v3/view/transform.hpp>
#include <range/v3/view/zip.hpp>

namespace stg::tensor {
    namespace bg = boost::geometry;

    template<std::floating_point T>
    class Tensor {
    public:
        using value_type = T;
        using iterator = std::array<T, 9>::iterator;
        using const_iterator = std::array<T, 9>::const_iterator;

        Tensor() = default;

        template<std::input_iterator Iter>
        Tensor(Iter begin, Iter end) : values_{begin, end} {}

        explicit constexpr Tensor(std::array<T, 9> tensor) : values_(std::move(tensor)) {}

        Tensor(T m11, T m12, T m13, T m21, T m22, T m23, T m31, T m32, T m33)
            : values_{m11, m12, m13, m21, m22, m23, m31, m32, m33} {}

        value_type get(size_t i, size_t j) const { return values_[i * 3 + j]; }

        void set(size_t i, size_t j, value_type value) { values_[i * 3 + j] = value; }

        bg::model::point<value_type, 3, bg::cs::cartesian>
        operator*(const bg::model::point<value_type, 3, bg::cs::cartesian>& vector) const {
            const value_type x = vector.template get<0>() * values_[0] + vector.template get<1>() * values_[1] + vector.template get<2>() * values_[2];
            const value_type y = vector.template get<0>() * values_[3] + vector.template get<1>() * values_[4] + vector.template get<2>() * values_[5];
            const value_type z = vector.template get<0>() * values_[6] + vector.template get<1>() * values_[7] + vector.template get<2>() * values_[8];
            return {x, y, z};
        }

        [[nodiscard("Returns new tensor")]] Tensor operator*(value_type coeff) const {
            Tensor result = *this;
            std::for_each(result.begin(), result.end(), [coeff](auto& elem) { elem *= coeff; });
            return result;
        }

        [[nodiscard("Returns new tensor")]] Tensor operator/(value_type coeff) const {
            Tensor result = *this;
            std::for_each(result.begin(), result.end(), [coeff](auto& elem) { elem /= coeff; });
            return result;
        }

        Tensor operator+(const Tensor& other) const {
            return {other.values_[0] + values_[0], other.values_[1] + values_[1], other.values_[2] + values_[2],
                    other.values_[3] + values_[3], other.values_[4] + values_[4], other.values_[5] + values_[5],
                    other.values_[6] + values_[6], other.values_[7] + values_[7], other.values_[8] + values_[8]};
        }

        constexpr auto operator<=>(const Tensor&) const = default;

        [[nodiscard("Returns new lower triangular tensor")]] constexpr Tensor cholesky() const {
            const value_type a_11 = std::sqrt(values_[0]);
            const value_type a_21 = values_[3] / a_11;
            const value_type a_22 = std::sqrt(values_[4] - a_21 * a_21);
            const value_type a_31 = values_[6] / a_11;
            const value_type a_32 = (values_[7] - a_31 * a_21) / a_22;
            const value_type a_33 = std::sqrt(values_[8] - a_31 * a_31 - a_32 * a_32);
            return Tensor{std::array{a_11, 0., 0.,
                                     a_21, a_22, 0.,
                                     a_31, a_32, a_33}};
        }

        value_type determinant() const;

        [[nodiscard]] constexpr size_t size() const { return size_; }

        iterator begin() { return values_.begin(); }

        iterator end() { return values_.end(); }

        const_iterator cbegin() const { return values_.cbegin(); }

        const_iterator cend() const { return values_.cend(); }

    private:
        std::array<value_type, 9> values_;
        static constexpr inline size_t size_ = 9;
    };


    template<std::floating_point T>
    Tensor<T>::value_type Tensor<T>::determinant() const {
        return values_[0] * values_[4] * values_[8] + values_[1] * values_[5] * values_[6] + values_[2] * values_[3] * values_[7] - values_[2] * values_[4] * values_[6] - values_[0] * values_[5] * values_[7] - values_[1] * values_[3] * values_[8];
    }
}// namespace stg::tensor

#endif//STG_TENSOR_HPP
