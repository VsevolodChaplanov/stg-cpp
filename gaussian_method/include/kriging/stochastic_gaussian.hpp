#ifndef _STOCHASTIC_GAUSSIAN_HPP_
#define _STOCHASTIC_GAUSSIAN_HPP_

#include "space.hpp"
#include "varfun.hpp"
#include <armadillo>
#include <array>
#include <functional>
#include <string>
#include <utility>

namespace stg::kriging {
    template<size_t Dim>
    class StochasticGaussian {
    public:
        struct Params {
            size_t eigen_cut = 1000;
            double variance_cut = 0.05;
        };

        StochasticGaussian(
                PhysicalSpace space,
                const IVarFun<Dim>& varfun,
                Params params) : _space(std::move(space)), _params(params) { initialize(varfun); }

        [[nodiscard]] const PhysicalSpace& space() const { return _space; }
        const Params& params() const { return _params; }

        std::array<std::vector<double>, Dim> generate(size_t seed) const;

        [[nodiscard]] std::string dstr() const;

    private:
        PhysicalSpace _space;
        Params _params;
        arma::Col<double> _eigval;
        arma::Mat<double> _eigvec;

        void initialize(const IVarFun<Dim>& varfun);
    };
}// namespace stg::kriging

#endif
