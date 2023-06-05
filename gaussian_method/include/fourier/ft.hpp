#ifndef _FT_HPP
#define _FT_HPP

#include "rtable/cube_relation_table.hpp"
#include "space.hpp"
#include <array>
#include <memory>
#include <mesh_builders/cube_fe_mesh.hpp>
#include <vector>

namespace stg::kriging {
    // N - partition in single direction. Total number of cells is N*N*N. Sould be odd.
    // f is given in physical space on the cube [-L/2, L/2] as cell_centered values
    // !! Only for real valued transformations
    //
    // => (1/2pi)^3 * int3{-L/2;L/2} ( f(x)*exp(-i*dot(x, k)*dx )
    std::vector<double> fourier3(PhysicalSpace ps, const std::vector<double>& f);

    // N - partition in single direction. Total number of cells is N*N*N. Sould be odd.
    // fk is given in fourier space on the cube [-Lk/2, Lk/2] as cell_centered values
    // !! Only for real valued transformations
    //
    // => int3{-Lk/2;Lk/2} ( fk(k)*exp(i*dot(x, k)*dk )
    std::vector<double> inverse_fourier3(FourierSpace fs, const std::vector<double>& fk);


        const auto fourier_space = ps.template make_fourier_space<mesh::CubeMeshBuilder<T>>();

        // fill input
        auto* input = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n3);
        for (size_t i = 0; i < n3; ++i) input[i][0] = input[i][1] = 0;

        for (size_t k = 0; k < n; ++k)
            for (size_t j = 0; j < n; ++j)
                for (size_t i = 0; i < n; ++i) {
                    size_t ind = ps.lin_index(i, j, k);
                    input[ind][0] = f[ind];
                }

        // allocate output
        auto* output = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n3);

        // build plan
        fftw_plan plan = fftw_plan_dft_3d(n, n, n, input, output, FFTW_FORWARD, FFTW_ESTIMATE);

        // exec
        fftw_execute(plan);

        // set values
        std::vector<T> ret(n3, 0);
        T nrm = std::pow(ps.h() / (2 * M_PI), 3);
        for (size_t k = 0; k < n; ++k)
            for (size_t j = 0; j < n; ++j)
                for (size_t i = 0; i < n; ++i) {
                    const auto& linear_coordinate_ft = fourier_space->linear_coordinates();
                    const auto& linear_coordinate = ps.linear_coordinates();
                    size_t ind = ps.lin_index(i, j, k);

                    size_t i2 = ((i + (n - 1) / 2) % n);
                    size_t j2 = ((j + (n - 1) / 2) % n);
                    size_t k2 = ((k + (n - 1) / 2) % n);
                    size_t ind2 = i2 + j2 * n + k2 * n * n;


                    T num_arg = (linear_coordinate_ft[i2] + linear_coordinate_ft[j2] + linear_coordinate_ft[k2]) * linear_coordinate[0];
                    T num_real = std::cos(num_arg);
                    T num_imag = -std::sin(num_arg);

                    T v_real = output[ind][0] * num_real - output[ind][1] * num_imag;
                    T v_imag = output[ind][0] * num_imag + output[ind][1] * num_real;

                    if (std::abs(v_imag) > 1e-4) {
                        throw std::runtime_error("nonzero imag part: " + std::to_string(v_imag));
                    }
                    ret[ind2] = v_real * nrm;
                }

        // free memory
        fftw_destroy_plan(plan);
        fftw_free(input);
        fftw_free(output);

        return ret;
    }

    template<std::floating_point T>
    std::vector<T> inverse_fourier3(const mesh::CubeRelationTable<T>& fs, const std::vector<T>& fk) {
        const std::size_t n = fs.n();
        const std::size_t n3 = fs.n_vertices();

        const auto real_space = fs.template make_real_space<mesh::CubeMeshBuilder<T>>();

        // fill input
        auto* input = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n3);
        for (size_t i = 0; i < n3; ++i) input[i][0] = input[i][1] = 0;

        for (size_t k = 0; k < n; ++k)
            for (size_t j = 0; j < n; ++j)
                for (size_t i = 0; i < n; ++i) {
                    size_t ind = fs.lin_index(i, j, k);
                    input[ind][0] = fk[ind];
                }

        // allocate output
        auto* output = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n3);

        // build plan
        fftw_plan plan = fftw_plan_dft_3d(n, n, n, input, output, FFTW_BACKWARD, FFTW_ESTIMATE);

        // exec
        fftw_execute(plan);

        // set values
        std::vector<T> ret(n3, 0);
        T nrm = std::pow(fs.h(), 3);
        for (size_t k = 0; k < n; ++k)
            for (size_t j = 0; j < n; ++j)
                for (size_t i = 0; i < n; ++i) {
                    const auto& linear_coordinate_ft = fs.linear_coordinates();
                    const auto& linear_coordinate = real_space->linear_coordinates();
                    size_t ind = i + j * n + k * n * n;

                    size_t i2 = ((i + (n - 1) / 2) % n);
                    size_t j2 = ((j + (n - 1) / 2) % n);
                    size_t k2 = ((k + (n - 1) / 2) % n);
                    size_t ind2 = i2 + j2 * n + k2 * n * n;


                    T num_arg = (linear_coordinate[i2] + linear_coordinate[j2] + linear_coordinate[k2]) * linear_coordinate_ft[0];
                    T num_real = std::cos(num_arg);
                    T num_imag = std::sin(num_arg);

                    T v_real = output[ind][0] * num_real - output[ind][1] * num_imag;
                    T v_imag = output[ind][0] * num_imag + output[ind][1] * num_real;

                    if (std::abs(v_imag) > 1e-4) {
                        throw std::runtime_error("nonzero imag part: " + std::to_string(v_imag));
                    }
                    ret[ind2] = v_real * nrm;
                }

        // free memory
        fftw_destroy_plan(plan);
        fftw_free(input);
        fftw_free(output);

        return ret;
    }
}// namespace stg::kriging

#endif
