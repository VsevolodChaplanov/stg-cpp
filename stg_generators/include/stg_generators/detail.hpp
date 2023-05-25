#ifndef STG_STG_GENERATOR_DETAIL_HPP
#define STG_STG_GENERATOR_DETAIL_HPP

#include <concepts>
namespace stg::generators::detail {

    // type needs to disambiguate diamond problem
    template<std::floating_point T, template<std::floating_point> typename Derived,
             template<std::floating_point, typename> typename crtp_type>
    class crtp_helper_with_type {
    public:
        constexpr Derived<T>& underlying() { return static_cast<Derived<T>&>(*this); }
        constexpr const Derived<T>& underlying() const { return static_cast<Derived<T> const&>(*this); }
        friend crtp_type<T, Derived<T>>;
    };

    template<template<std::floating_point> typename Derived>
    class crtp_helper {
    public:
        template<std::floating_point T>
        constexpr Derived<T>& underlying() { return static_cast<Derived<T>&>(*this); }

        template<std::floating_point T>
        constexpr const Derived<T>& underlying() const { return static_cast<Derived<T> const&>(*this); }
    };
}// namespace stg::generators::detail

#endif