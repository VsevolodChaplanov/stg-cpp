#ifndef STG_FOURIER_HPP
#define STG_FOURIER_HPP

#include <concepts>
#include <fftw3.h>
#include <geometry/geometry.hpp>
#include <memory>
#include <range/v3/view/iota.hpp>
#include <rtable/cube_relation_table.hpp>
#include <vector>

namespace stg::fourier {
    namespace rv = ranges::views;

    template<std::floating_point T>
    std::vector<T> fourier3(const std::shared_ptr<mesh::CubeRelationTable<T>>& space, const std::vector<T>& f) {
        std::size_t n_l = space->n();
        std::size_t n_vert = space->n_vertices();
        const auto fourier_space = space->make_fourier_space();

        auto* input = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n_vert);
        for (const std::size_t index: rv::iota(0ul, n_vert))
            input[index][0] = input[index][1] = 0;

        // mb from 0 to n_vert?
        for (const std::size_t k: rv::iota(0ul, n_l))
            for (const std::size_t j: rv::iota(0ul, n_l))
                for (const std::size_t i: rv::iota(0ul, n_l)) {
                    std::size_t ind = space->lin_index(i, j, k);
                    input[ind][0] = f[ind];
                }

        auto* output = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n_vert);

        fftw_plan plan = fftw_plan_dft_3d(static_cast<int>(n_l), static_cast<int>(n_l), static_cast<int>(n_l), input, output,
                                          FFTW_FORWARD, FFTW_ESTIMATE);

        fftw_execute(plan);

        std::vector<T> ret(n_vert, 0);
        const T nrm = std::pow(space->h() / (2 * std::numbers::pi_v<T>), 3);
        const T init_value = -space->l() / 2;
        for (const std::size_t k: rv::iota(0ul, n_l))
            for (const std::size_t j: rv::iota(0ul, n_l))
                for (const std::size_t i: rv::iota(0ul, n_l)) {
                    size_t ind = space->lin_index(i, j, k);

                    size_t i2 = ((i + (n_l - 1) / 2) % n_l);
                    size_t j2 = ((j + (n_l - 1) / 2) % n_l);
                    size_t k2 = ((k + (n_l - 1) / 2) % n_l);
                    size_t ind2 = i2 + j2 * n_l + k2 * n_l * n_l;

                    const auto fourier_vert = fourier_space->vertex(i2, j2, k2);
                    const T x_f = fourier_vert.template get<0>();
                    const T y_f = fourier_vert.template get<1>();
                    const T z_f = fourier_vert.template get<2>();

                    double num_arg = (x_f + y_f + z_f) * init_value;
                    double num_real = cos(num_arg);
                    double num_imag = -sin(num_arg);

                    double v_real = output[ind][0] * num_real - output[ind][1] * num_imag;
                    double v_imag = output[ind][0] * num_imag + output[ind][1] * num_real;

                    if (std::abs(v_imag) > 1e-8) {
                        throw std::runtime_error("nonzero imag part: " + std::to_string(v_imag));
                    }
                    ret[ind2] = v_real * nrm;
                }

        fftw_destroy_plan(plan);
        fftw_free(input);
        fftw_free(output);

        return ret;
    }
}// namespace stg::fourier

#endif//STG_FOURIER_HPP
