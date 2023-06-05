#ifndef STG_RECTILINEAR_GRID_PARSER_HPP
#define STG_RECTILINEAR_GRID_PARSER_HPP

#include "geometry/geometry.hpp"
#include "stg_tensor/tensor.hpp"
#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>
#include <concepts>
#include <filesystem>
#include <fmt/format.h>
#include <fstream>
#include <iostream>
#include <iterator>
#include <mesh_builders/cube_fe_mesh.hpp>
#include <mesh_builders/mesh_builders.hpp>
#include <rtable/cube_relation_table.hpp>
#include <string>
#include <string_view>
#include <vector>

namespace stg::mesh {

    /*
   * Пока только кубическая сетка
   */
    class RectilinearGridParser final {
    public:
        explicit RectilinearGridParser(std::string_view filename)
            : file_{open_file(filename)} {}

        RectilinearGridParser(const RectilinearGridParser&) = delete;
        RectilinearGridParser& operator=(const RectilinearGridParser&) = delete;

        // template<std::floating_point T>
        // std::shared_ptr<CubeRelationTable<T>> mesh() {
        //     boost::char_separator<char> sep{" "};
        //     std::string line;
        //     bool coordinates_begin = false;
        //     bool has_n_points = false;
        //     bool has_length = false;
        //     std::size_t n;
        //     double l;
        //     std::vector<T> vertices;

        //     while (std::getline(file_, line)) {
        //         boost::tokenizer tok(line, sep);
        //         if (*tok.begin() == "DIMENSIONS") {
        //             n = std::stoul(*(++tok.begin()));
        //             cached_vert_number_ = n * n * n;
        //             has_n_points = true;
        //         }

        //         if (*tok.begin() == "Y_COORDINATES") {
        //             coordinates_begin = false;
        //             has_length = true;
        //         }

        //         if (coordinates_begin) {
        //             const auto coordinate_left = *tok.begin();
        //             const auto coord = std::stod(coordinate_left);
        //             vertices.push_back(coord);
        //         }

        //         if (*tok.begin() == "X_COORDINATES") {
        //             coordinates_begin = true;
        //         }

        //         if (has_n_points && has_length) { break; }
        //     }

        //     l = 2 * vertices.back();

        //     CubeMeshBuilder<T> builder{l, n};
        //     return builder.build_relation_table(std::move(vertices));
        // }

        // template<std::floating_point T>
        // std::shared_ptr<CubeFiniteElementsMesh<T>> fe_mesh() {
        //     boost::char_separator<char> sep{" "};
        //     std::string line;
        //     bool coordinates_begin = false;
        //     bool has_n_points = false;
        //     bool has_length = false;
        //     std::size_t n;
        //     double l;
        //     std::vector<T> vertices;

        //     while (std::getline(file_, line)) {
        //         boost::tokenizer tok(line, sep);
        //         if (*tok.begin() == "DIMENSIONS") {
        //             n = std::stoul(*(++tok.begin()));
        //             cached_vert_number_ = n * n * n;
        //             has_n_points = true;
        //             has_cached_vert_number_ = true;
        //         }

        //         if (*tok.begin() == "Y_COORDINATES") {
        //             coordinates_begin = false;
        //             has_length = true;
        //         }

        //         if (coordinates_begin) {
        //             const auto coordinate_left = *tok.begin();
        //             const auto coord = std::stod(coordinate_left);
        //             vertices.push_back(coord);
        //         }

        //         if (*tok.begin() == "X_COORDINATES") {
        //             coordinates_begin = true;
        //         }

        //         if (has_n_points && has_length) { break; }
        //     }

        //     l = 2 * vertices.back();

        //     CubeMeshBuilder<T> builder{l, n};
        //     return builder.build(std::move(vertices));
        // }

        template<std::floating_point T>
        std::vector<T> scalar_data(std::string_view starting_expr = "LOOKUP_TABLE default",
                                   std::string_view end_expr = " ") {
            bool scalar_data_start = false;
            boost::char_separator<char> sep{" "};
            std::string line;
            std::vector<T> result;
            if (has_cached_vert_number_) result.reserve(cached_vert_number_);
            while (std::getline(file_, line)) {
                if (scalar_data_start) {
                    result.push_back(std::stod(line));
                }
                if (line == starting_expr) {
                    scalar_data_start = true;
                }
            }
            return result;
        }

        template<std::floating_point T>
        std::vector<Vector<T>> vector_data(std::string_view starting_expr = "LOOKUP_TABLE default",
                                           std::string_view end_expr = " ") {
            bool vector_data_start = false;
            boost::char_separator<char> sep{" "};
            std::string line;
            std::vector<Vector<T>> result;
            if (has_cached_vert_number_) result.reserve(cached_vert_number_);
            while (std::getline(file_, line)) {
                boost::tokenizer tok{line, sep};
                if (vector_data_start) {
                    std::vector<T> vel(3);
                    std::transform(tok.begin(), tok.end(), vel.begin(), [](const std::string& value) {
                        return std::stod(value);
                    });
                    result.push_back({vel[0], vel[1], vel[2]});
                }
                if (line == starting_expr) {
                    vector_data_start = true;
                }
            }
            return result;
        }

        template<std::floating_point T>
        std::vector<tensor::Tensor<T>> tensor_data(std::string_view starting_expr = "LOOKUP_TABLE default",
                                                   std::string_view end_expr = " ") {
            bool tensor_data_start = false;
            boost::char_separator<char> sep{" "};
            std::string line;
            std::vector<tensor::Tensor<T>> result;
            if (has_cached_vert_number_) result.reserve(cached_vert_number_);
            while (std::getline(file_, line)) {
                boost::tokenizer tok{line, sep};
                if (tensor_data_start) {
                    tensor::Tensor<T> tensor;
                    std::transform(tok.begin(), tok.end(), tensor.begin(), [](const std::string& value) {
                        return std::stod(value);
                    });
                    result.push_back(tensor);
                }
                if (line == starting_expr) {
                    tensor_data_start = true;
                }
            }
            return result;
        }

        template<std::floating_point T>
        std::shared_ptr<CubeRelationTable<T>> mesh() {
            boost::char_separator<char> sep{" "};
            std::string line;
            bool coordinates_begin = false;
            bool has_n_points = false;
            bool has_length = false;
            std::size_t n;
            double l;

            while (std::getline(file_, line)) {
                boost::tokenizer tok(line, sep);
                if (*tok.begin() == "DIMENSIONS") {
                    n = std::stoul(*(++tok.begin()));
                    cached_vert_number_ = n * n * n;
                    has_n_points = true;
                }

                if (coordinates_begin) {
                    const auto coordinate_left = *tok.begin();
                    l = 2 * std::fabs(std::stod(coordinate_left));
                    const auto dl = l / n;
                    l += dl;
                    has_length = true;
                }

                if (*tok.begin() == "X_COORDINATES") {
                    coordinates_begin = true;
                }

                if (has_n_points && has_length) { break; }
            }

            CubeMeshBuilder<T> builder{l, n};
            return builder.build_relation_table();
        }

        template<std::floating_point T>
        std::shared_ptr<CubeFiniteElementsMesh<T>> fe_mesh() {
            boost::char_separator<char> sep{" "};
            std::string line;
            bool coordinates_begin = false;
            bool has_n_points = false;
            bool has_length = false;
            std::size_t n;
            double l;

            while (std::getline(file_, line)) {
                boost::tokenizer tok(line, sep);
                if (*tok.begin() == "DIMENSIONS") {
                    n = std::stoul(*(++tok.begin()));
                    cached_vert_number_ = n * n * n;
                    has_n_points = true;
                    has_cached_vert_number_ = true;
                }

                if (coordinates_begin) {
                    const auto coordinate_left = *tok.begin();
                    l = 2 * std::fabs(std::stod(coordinate_left));
                    const auto dl = l / n;
                    l += dl;
                    has_length = true;
                }

                if (*tok.begin() == "X_COORDINATES") {
                    coordinates_begin = true;
                }

                if (has_n_points && has_length) { break; }
            }

            CubeMeshBuilder<T> builder{l, n};
            return builder.build();
        }

    private:
        std::ifstream file_;
        std::size_t cached_vert_number_;
        bool has_cached_vert_number_ = false;

        static std::ifstream open_file(std::string_view filename) {
            std::filesystem::path path = std::filesystem::absolute(filename);
            if (!std::filesystem::exists(path)) {
                fmt::print(std::cerr, "File doesn't exists {}", filename);
            }
            return std::ifstream{path, std::ios_base::in};
        }
    };
}// namespace stg::mesh

#endif//STG_RECTILINEAR_GRID_PARSER_HPP
