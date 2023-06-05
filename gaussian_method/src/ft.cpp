#include "ft.hpp"
#include "mesh_builders/mesh_builders.hpp"
#include <cassert>
#include <cmath>
#include <fftw3.h>
#include <iostream>
#include <stdexcept>

namespace stg::kriging {
    std::vector<double> fourier3(PhysicalSpace ps, const std::vector<double>& f) {
        size_t N = ps.N;
        size_t n3 = N * N * N;
        if (f.size() != n3) {
            throw std::runtime_error("Invalid length of f array");
        }
        FourierSpace fs = ps.fourier_space();

        // fill input
        auto* input = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n3);
        for (size_t i = 0; i < n3; ++i) input[i][0] = input[i][1] = 0;

        for (size_t k = 0; k < N; ++k)
            for (size_t j = 0; j < N; ++j)
                for (size_t i = 0; i < N; ++i) {
                    size_t ind = ps.lin_index(i, j, k);
                    input[ind][0] = f[ind];
                }

        // allocate output
        auto* output = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n3);

        // build plan
        fftw_plan plan = fftw_plan_dft_3d(N, N, N, input, output, FFTW_FORWARD, FFTW_ESTIMATE);

        // exec
        fftw_execute(plan);

        // set values
        std::vector<double> ret(n3, 0);
        double nrm = std::pow(ps.h / (2 * M_PI), 3);
        for (size_t k = 0; k < N; ++k)
            for (size_t j = 0; j < N; ++j)
                for (size_t i = 0; i < N; ++i) {
                    size_t ind = ps.lin_index(i, j, k);

                    size_t i2 = ((i + (N - 1) / 2) % N);
                    size_t j2 = ((j + (N - 1) / 2) % N);
                    size_t k2 = ((k + (N - 1) / 2) % N);
                    size_t ind2 = i2 + j2 * N + k2 * N * N;

                    double num_arg = (fs.coo[i2] + fs.coo[j2] + fs.coo[k2]) * ps.coo[0];
                    double num_real = cos(num_arg);
                    double num_imag = -sin(num_arg);

                    double v_real = output[ind][0] * num_real - output[ind][1] * num_imag;
                    double v_imag = output[ind][0] * num_imag + output[ind][1] * num_real;

                    if (std::abs(v_imag) > 1e-8) {
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

    // perform inverse fourier transform
    std::vector<double> inverse_fourier3(FourierSpace fs, const std::vector<double>& f) {
        size_t N = fs.N;
        size_t n3 = N * N * N;
        if (f.size() != n3) {
            throw std::runtime_error("Invalid length of f array");
        }
        PhysicalSpace ps = fs.physical_space();

        // fill input
        auto* input = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n3);
        for (size_t i = 0; i < n3; ++i) input[i][0] = input[i][1] = 0;

        for (size_t k = 0; k < N; ++k)
            for (size_t j = 0; j < N; ++j)
                for (size_t i = 0; i < N; ++i) {
                    size_t ind = i + N * j + N * N * k;
                    input[ind][0] = f[ind];
                }

        // allocate output
        auto* output = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n3);

        // build plan
        fftw_plan plan = fftw_plan_dft_3d(N, N, N, input, output, FFTW_BACKWARD, FFTW_ESTIMATE);

        // exec
        fftw_execute(plan);

        // set values
        std::vector<double> ret(N * N * N, 0);
        double nrm = std::pow(fs.h, 3);
        for (size_t k = 0; k < N; ++k)
            for (size_t j = 0; j < N; ++j)
                for (size_t i = 0; i < N; ++i) {
                    size_t ind = i + j * N + k * N * N;

                    size_t i2 = ((i + (N - 1) / 2) % N);
                    size_t j2 = ((j + (N - 1) / 2) % N);
                    size_t k2 = ((k + (N - 1) / 2) % N);
                    size_t ind2 = i2 + j2 * N + k2 * N * N;

                    double num_arg = (ps.coo[i2] + ps.coo[j2] + ps.coo[k2]) * fs.coo[0];
                    double num_real = cos(num_arg);
                    double num_imag = sin(num_arg);

                    double v_real = output[ind][0] * num_real - output[ind][1] * num_imag;
                    double v_imag = output[ind][0] * num_imag + output[ind][1] * num_real;

                    if (std::abs(v_imag) > 1e-8) {
                        throw std::runtime_error("nonzero imag part: " + std::to_string(v_imag));
                    }
                    ret[ind2] = v_real * nrm;
                }

        // free memory
        fftw_destroy_plan(plan);
        fftw_free(input);
        fftw_free(output);

        return ret;
    };

    // std::vector<double> fourier3(const mesh::CubeRelationTable<double>& ps, const std::vector<double>& f) {
    //     const std::size_t n = ps.n();
    //     const std::size_t n3 = ps.n_vertices();

    //     const auto fourier_space = ps.make_fourier_space<mesh::CubeMeshBuilder<double>>();

    //     // fill input
    //     auto* input = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n3);
    //     for (size_t i = 0; i < n3; ++i) input[i][0] = input[i][1] = 0;

    //     for (size_t k = 0; k < n; ++k)
    //         for (size_t j = 0; j < n; ++j)
    //             for (size_t i = 0; i < n; ++i) {
    //                 size_t ind = ps.lin_index(i, j, k);
    //                 input[ind][0] = f[ind];
    //             }

    //     // allocate output
    //     auto* output = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n3);

    //     // build plan
    //     fftw_plan plan = fftw_plan_dft_3d(n, n, n, input, output, FFTW_FORWARD, FFTW_ESTIMATE);

    //     // exec
    //     fftw_execute(plan);

    //     // set values
    //     std::vector<double> ret(n3, 0);
    //     double nrm = std::pow(ps.h() / (2 * M_PI), 3);
    //     for (size_t k = 0; k < n; ++k)
    //         for (size_t j = 0; j < n; ++j)
    //             for (size_t i = 0; i < n; ++i) {
    //                 const auto& linear_coordinate_ft = fourier_space->linear_coordinates();
    //                 const auto& linear_coordinate = ps.linear_coordinates();
    //                 size_t ind = ps.lin_index(i, j, k);

    //                 size_t i2 = ((i + (n - 1) / 2) % n);
    //                 size_t j2 = ((j + (n - 1) / 2) % n);
    //                 size_t k2 = ((k + (n - 1) / 2) % n);
    //                 size_t ind2 = i2 + j2 * n + k2 * n * n;


    //                 double num_arg = (linear_coordinate_ft[i2] + linear_coordinate_ft[j2] + linear_coordinate_ft[k2]) * linear_coordinate[0];
    //                 double num_real = std::cos(num_arg);
    //                 double num_imag = -std::sin(num_arg);

    //                 double v_real = output[ind][0] * num_real - output[ind][1] * num_imag;
    //                 double v_imag = output[ind][0] * num_imag + output[ind][1] * num_real;

    //                 if (std::abs(v_imag) > 1e-8) {
    //                     throw std::runtime_error("nonzero imag part: " + std::to_string(v_imag));
    //                 }
    //                 ret[ind2] = v_real * nrm;
    //             }

    //     // free memory
    //     fftw_destroy_plan(plan);
    //     fftw_free(input);
    //     fftw_free(output);

    //     return ret;
    // }

    // std::vector<double> inverse_value_typefourier3(const mesh::CubeRelationTable<double>& fs, const std::vector<double>& fk) {
    //     const std::size_t n = fs.n();
    //     const std::size_t n3 = fs.n_vertices();

    //     const auto real_space = fs.make_real_space<mesh::CubeMeshBuilder<double>>();

    //     // fill input
    //     auto* input = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n3);
    //     for (size_t i = 0; i < n3; ++i) input[i][0] = input[i][1] = 0;

    //     for (size_t k = 0; k < n; ++k)
    //         for (size_t j = 0; j < n; ++j)
    //             for (size_t i = 0; i < n; ++i) {
    //                 size_t ind = fs.lin_index(i, j, k);
    //                 input[ind][0] = fk[ind];
    //             }

    //     // allocate output
    //     auto* output = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n3);

    //     // build plan
    //     fftw_plan plan = fftw_plan_dft_3d(n, n, n, input, output, FFTW_BACKWARD, FFTW_ESTIMATE);

    //     // exec
    //     fftw_execute(plan);

    //     // set values
    //     std::vector<double> ret(n3, 0);
    //     double nrm = std::pow(fs.h(), 3);
    //     for (size_t k = 0; k < n; ++k)
    //         for (size_t j = 0; j < n; ++j)
    //             for (size_t i = 0; i < n; ++i) {
    //                 const auto& linear_coordinate_ft = fs.linear_coordinates();
    //                 const auto& linear_coordinate = real_space->linear_coordinates();
    //                 size_t ind = i + j * n + k * n * n;

    //                 size_t i2 = ((i + (n - 1) / 2) % n);
    //                 size_t j2 = ((j + (n - 1) / 2) % n);
    //                 size_t k2 = ((k + (n - 1) / 2) % n);
    //                 size_t ind2 = i2 + j2 * n + k2 * n * n;


    //                 double num_arg = (linear_coordinate[i2] + linear_coordinate[j2] + linear_coordinate[k2]) * linear_coordinate_ft[0];
    //                 double num_real = std::cos(num_arg);
    //                 double num_imag = std::sin(num_arg);

    //                 double v_real = output[ind][0] * num_real - output[ind][1] * num_imag;
    //                 double v_imag = output[ind][0] * num_imag + output[ind][1] * num_real;

    //                 if (std::abs(v_imag) > 1e-8) {
    //                     throw std::runtime_error("nonzero imag part: " + std::to_string(v_imag));
    //                 }
    //                 ret[ind2] = v_real * nrm;
    //             }

    //     // free memory
    //     fftw_destroy_plan(plan);
    //     fftw_free(input);
    //     fftw_free(output);

    //     return ret;
    // }
}// namespace stg::kriging