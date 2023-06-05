#ifndef STG_FOURIER_FT_HPP
#define STG_FOURIER_FT_HPP

#include "../kriging/space.hpp"
#include "fftw3.h"
#include "mesh_builders/mesh_builders.hpp"
#include "rtable/cube_relation_table.hpp"
#include <array>
#include <concepts>
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

    template<std::floating_point T>
    std::vector<T> fourier3(const mesh::CubeRelationTable<T>& ps, const std::vector<T>& f) {
        const std::size_t n = ps.n();
        const std::size_t n3 = ps.n_vertices();

    std::vector<double> fourier3(const mesh::CubeRelationTable<double>& ps, const std::vector<double>& f);
    std::vector<double> inverse_fourier3(const mesh::CubeRelationTable<double>& fs, const std::vector<double>& fk);
}// namespace stg::kriging

#endif
