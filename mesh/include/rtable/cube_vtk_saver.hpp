#ifndef STG_CUBE_VTK_SAVER_HPP
#define STG_CUBE_VTK_SAVER_HPP

#include <string_view>
#include <execution>
#include <string_view>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <iostream>
#include <fmt/core.h>
#include <fmt/ostream.h>
#include <fmt/format.h>
#include <geometry/geometry.hpp>
#include <stg_tensor/tensor.hpp>
#include "cube_relation_table.hpp"


namespace stg::mesh {

  class VtkRectilinearGridSaver final {
  public:
    VtkRectilinearGridSaver(std::string_view filename)
      : file_{open_file(filename)} { }

    template<std::floating_point T>
    void save_mesh(const std::shared_ptr<CubeRelationTable<T>>& rtable) {
      const auto& vertices = rtable->vertices();
      const auto& bound_indices = rtable->elements_global_indices();
      const auto& elements_types = rtable->elements_types();

      write_header(rtable->n());
      write_coordinates(vertices.cbegin(), vertices.cend(), vertices.size());
/*      write_cell_data(bound_indices.cbegin(), bound_indices.cend(), rtable->n_elements());
      write_cell_types_data(elements_types.cbegin(), elements_types.cend(), rtable->n_elements());*/
    }

    template<std::forward_iterator Iter>
    void save_scalar_data(Iter begin, Iter end,
                          std::string_view table_name = "DefaultTable") {
      using IterValueType = typename std::iterator_traits<Iter>::value_type;
      const std::size_t size = std::distance(begin, end);
      write_point_data_header(size);
      fmt::print(file_, "SCALARS {} double\n", table_name);
      fmt::print(file_, "LOOKUP_TABLE default\n");
      std::copy(begin, end, std::ostream_iterator<IterValueType>{file_, "\n"});
      file_ << std::endl;
    }

    template<std::forward_iterator Iter>
    void save_vector_data(Iter begin, Iter end,
                          std::string_view table_name = "VectorField") {
      const std::size_t size = std::distance(begin, end);
      write_point_data_header(size);
      fmt::print(file_, "VECTORS {} double\n", table_name);
      std::for_each(begin, end, [&](const auto& vector) {
        fmt::print(file_, "{} {} {}\n",
                   vector.template get<0>(),
                   vector.template get<1>(),
                   vector.template get<2>());
      });
      file_ << std::endl;
    }

    template<ranges::viewable_range Range>
    void save_velocity_data(Range&& range,
                            std::string_view table_name = "VectorField") {
      const std::size_t size = ranges::distance(range);
      write_point_data_header(size);
      fmt::print(file_, "VECTORS {} double\n", table_name);
      ranges::for_each(range, [&](const auto& vec_zip) {
                                const auto [x, y, z] = vec_zip;
                                fmt::print(file_, "{} {} {}\n", x, y, z);
                              });
      file_ << std::endl;
    }

    template<std::forward_iterator Iter>
    void save_tensor_data(Iter begin, Iter end,
                          std::string_view table_name = "TensorData") {
      using IterValueType = typename std::iterator_traits<Iter>::value_type;
      using TensorValueType = typename IterValueType::value_type;
      const std::size_t size = std::distance(begin, end);
      write_point_data_header(size);
      fmt::print(file_, "TENSORS {} double\n", table_name);
      std::for_each(begin, end, [&](const auto& tensor) {
        std::copy(tensor.cbegin(), tensor.cend(), std::ostream_iterator<TensorValueType>{file_, " "});
        file_ << std::endl;
      });
      file_ << std::endl;
    }

    ~VtkRectilinearGridSaver() {
      // Уточнить насчёт RAII
      if (file_.is_open()) {
        file_.close();
      }
    }

  private:
    std::ofstream file_;
    mutable bool has_point_data_flag_ = false;

    static std::ofstream open_file(std::string_view filename) {
      std::filesystem::path path = std::filesystem::absolute(filename);
      if (!std::filesystem::exists(path)) {
        std::filesystem::create_directories(path.parent_path());
      }
      return std::ofstream{path, std::ios_base::out};
    }

    void write_header(std::size_t dimensions) {
      fmt::print(file_, "# vtk DataFile Version 2.0\n"
                        "Function\n"
                        "ASCII\n"
                        "DATASET RECTILINEAR_GRID\n"
                        "DIMENSIONS {} {} {}\n",
                        dimensions, dimensions, dimensions
      );
    }

    void write_point_data_header(std::size_t size) {
      if (!has_point_data_flag_) {
        fmt::print(file_, "POINT_DATA {}\n", size);
        has_point_data_flag_ = true;
      }
    }

    template<std::forward_iterator Iter>
    void write_coordinates(Iter begin, Iter end, std::size_t points_num) {
      write_coordinate_component(begin, end, points_num, "X");
      write_coordinate_component(begin, end, points_num, "Y");
      write_coordinate_component(begin, end, points_num, "Z");
    }

    template<std::forward_iterator Iter>
    void write_coordinate_component(Iter begin, Iter end,
                                    std::size_t points_num,
                                    std::string_view component) {
      using IterValueType = typename std::iterator_traits<Iter>::value_type;
      fmt::print(file_, "{}_COORDINATES {} double\n", component, points_num);
      std::copy(begin, end, std::ostream_iterator<IterValueType>{file_, "\n"});
      file_ << std::endl;
    }

    template<std::forward_iterator Iter>
    void write_cell_data(Iter begin, Iter end, size_t elements_amount) {
      fmt::print(file_, "CELLS {} {}\n", elements_amount, bounded_vertices(begin, end));
      std::for_each(begin, end, [&](const auto& value) {
        file_ << value.size() << " ";
        std::copy(value.cbegin(), value.cend(), std::ostream_iterator<size_t>{file_, " "});
        file_ << std::endl;
      });
      file_ << std::endl;
    }

    template<std::forward_iterator Iter>
    void write_cell_types_data(Iter begin, Iter end, size_t elements_amount) {
      fmt::print(file_, "CELL_TYPES {}\n", elements_amount);
      std::copy(begin, end, std::ostream_iterator<size_t>{file_, "\n"});
      file_ << std::endl;
    }

    template<std::forward_iterator Iter>
    size_t bounded_vertices(Iter begin, Iter end) const {
      auto result = std::transform_reduce(std::execution::par_unseq,
                                          begin, end, 0ul,
                                          [](auto left, auto right) {
                                            return left + right + 1ul;
                                          },
                                          [](const auto& elem) {
                                            return elem.size();
                                          });
      return result;
    }
  };
}
#endif //STG_CUBE_VTK_SAVER_HPP
