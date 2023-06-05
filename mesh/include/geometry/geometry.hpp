#ifndef STG_GEOMETRY_HPP
#define STG_GEOMETRY_HPP
#include <boost/geometry/algorithms/area.hpp>
#include <boost/geometry/algorithms/assign.hpp>
#include <boost/geometry/core/cs.hpp>
#include <boost/geometry/geometries/concepts/point_concept.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <concepts>
#include <random>
#include <stg_tensor/tensor.hpp>

namespace bg = boost::geometry;
/* Contains geometric primitives based on boost::geometry */
namespace stg {
    template<std::floating_point T>
    using Point = bg::model::point<T, 3, bg::cs::cartesian>;

    template<std::floating_point T>
    using Vector = bg::model::point<T, 3, bg::cs::cartesian>;

    template<std::floating_point T>
    using Triangle = bg::model::polygon<Point<T>>;

    template<std::floating_point T>
    using Quad = bg::model::polygon<Point<T>>;

    template<std::floating_point T>
    void operator+=(Point<T>& one, const Point<T>& other) {
        const auto x = one.template get<0>();
        const auto y = one.template get<1>();
        const auto z = one.template get<2>();

        const auto dx = other.template get<0>();
        const auto dy = other.template get<1>();
        const auto dz = other.template get<2>();

        one.template set<0>(x + dx);
        one.template set<1>(y + dy);
        one.template set<2>(z + dz);
    }

    template<std::floating_point T>
    Point<T> operator+(const Point<T>& one, const Point<T>& other) {
        return Point<T>{
                one.template get<0>() + other.template get<0>(),
                one.template get<1>() + other.template get<1>(),
                one.template get<2>() + other.template get<2>(),
        };
    }

    template<std::floating_point T>
    Point<T> operator-(const Point<T>& one, const Point<T>& other) {
        return Point<T>{
                one.template get<0>() - other.template get<0>(),
                one.template get<1>() - other.template get<1>(),
                one.template get<2>() - other.template get<2>(),
        };
    }

    template<std::floating_point T>
    Point<T> operator/(const Point<T>& one, const T& value) {
        return Point<T>{
                one.template get<0>() / value,
                one.template get<1>() / value,
                one.template get<2>() / value,
        };
    }

    template<std::floating_point T>
    Point<T> operator*(const Point<T>& one, const T& value) {
        return Point<T>{
                one.template get<0>() * value,
                one.template get<1>() * value,
                one.template get<2>() * value,
        };
    }

    template<std::floating_point T>
    void operator*=(Point<T>& one, const T& value) {
        const double x = one.template get<0>() * value;
        const double y = one.template get<1>() * value;
        const double z = one.template get<2>() * value;

        one.template set<0>(x);
        one.template set<1>(y);
        one.template set<2>(z);
    }

    template<std::floating_point T>
    T dot_product(const Vector<T>& first, const Vector<T>& second) {
        const T xx = first.template get<0>() * second.template get<0>();
        const T yy = first.template get<1>() * second.template get<1>();
        const T zz = first.template get<2>() * second.template get<2>();
        return xx + yy + zz;
    }

    template<std::floating_point T>
    Vector<T> cross_product(const Vector<T>& first, const Vector<T>& second) {
        const T xx = first.template get<1>() * second.template get<2>() - first.template get<2>() * second.template get<1>();
        const T yy = first.template get<2>() * second.template get<0>() - first.template get<0>() * second.template get<2>();
        const T zz = first.template get<0>() * second.template get<1>() - first.template get<1>() * second.template get<0>();
        return {xx, yy, zz};
    }

    template<std::floating_point T>
    Vector<T> scale_to_length(const Vector<T>& vector, T length) {
        const auto ex_len = std::sqrt(dot_product(vector, vector));
        return vector * length / ex_len;
    }

    template<std::floating_point T>
    Vector<T> operator*(const tensor::Tensor<T>& matrix, const Vector<T>& vector) {
        const T first = matrix.template get(0, 0) + vector.template get<0>() + matrix.template get(0, 1) + vector.template get<1>() + matrix.template get(0, 2) + vector.template get<2>();
        const T second = matrix.template get(1, 0) + vector.template get<0>() + matrix.template get(1, 1) + vector.template get<1>() + matrix.template get(1, 2) + vector.template get<2>();
        const T third = matrix.template get(2, 0) + vector.template get<0>() + matrix.template get(2, 1) + vector.template get<1>() + matrix.template get(2, 2) + vector.template get<2>();
        return {first, second, third};
    }
}// namespace stg

#endif//STG_GEOMETRY_HPP
