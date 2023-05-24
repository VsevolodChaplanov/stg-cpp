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
#include <fmt/os.h>
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

      write_header(rtable->n());
      write_coordinates(vertices.cbegin(), vertices.cend(), vertices.size());
    }

    template<std::forward_iterator Iter>
    void save_scalar_data(Iter begin, Iter end,
                          std::string_view table_name = "DefaultTable") {
      const std::size_t size = std::distance(begin, end);
      write_point_data_header(size);
      file_.print("SCALARS {} double\n", table_name);
      file_.print("LOOKUP_TABLE default\n");
      file_.print("{}", fmt::join(begin, end, "\n"));
      file_.print("\n");
    }

    template<std::forward_iterator Iter>
    void save_vector_data(Iter begin, Iter end,
                          std::string_view table_name = "VectorField") {
      const std::size_t size = std::distance(begin, end);
      write_point_data_header(size);
      file_.print("VECTORS {} double\n", table_name);
      std::for_each(begin, end, [&](const auto& vector) {
        file_.print("{} {} {}\n",
                   vector.template get<0>(),
                   vector.template get<1>(),
                   vector.template get<2>());
      });
      file_.print("\n");
    }

    template<std::ranges::viewable_range Range>
    void save_vector_data(Range&& range, std::string_view table_name = "VectorField") {
      const std::size_t size = std::ranges::distance(range);
      write_point_data_header(size);
      file_.print("VECTORS {} double\n", table_name);
      std::ranges::for_each(std::forward<Range>(range),
        [&](const auto& vector) {
          file_.print("{} {} {}\n", std::get<0>(vector), std::get<1>(vector), std::get<2>(vector));
      });
      file_.print("\n");
    }

    template<ranges::viewable_range Range>
    void save_velocity_data(Range&& range,
                            std::string_view table_name = "VectorField") {
      const std::size_t size = ranges::distance(range);
      write_point_data_header(size);
      file_.print("VECTORS {} double\n", table_name);
      ranges::for_each(range, [&](const auto& vec_zip) {
                                const auto [x, y, z] = vec_zip;
                                file_.print("{} {} {}\n", x, y, z);
                              });
      file_.print("\n");
    }

    template<std::forward_iterator Iter>
    void save_tensor_data(Iter begin, Iter end,
                          std::string_view table_name = "TensorData") {
      const std::size_t size = std::distance(begin, end);
      write_point_data_header(size);
      file_.print("TENSORS {} double\n", table_name);
      std::for_each(begin, end, [&](const auto& tensor) {
        file_.print("{}", fmt::join(tensor.cbegin(), tensor.cend(), " "));
        file_.print("\n");
      });
      file_.print("\n");
    }

    ~VtkRectilinearGridSaver() { }

  private:
    fmt::ostream file_;
    mutable bool has_point_data_flag_ = false;

    fmt::ostream open_file(std::string_view filename) {
      return fmt::output_file(filename.data(),
        fmt::file::CREATE | fmt::file::WRONLY | fmt::file::APPEND);
    }

    void write_header(std::size_t dimensions) {
      file_.print("# vtk DataFile Version 2.0\n"
                        "Function\n"
                        "ASCII\n"
                        "DATASET RECTILINEAR_GRID\n"
                        "DIMENSIONS {} {} {}\n",
                        dimensions, dimensions, dimensions
      );
    }

    void write_point_data_header(std::size_t size) {
      if (!has_point_data_flag_) {
        file_.print("POINT_DATA {}\n", size);
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
      file_.print("{}_COORDINATES {} double\n", component, points_num);
      file_.print("{}", fmt::join(begin, end, "\n"));
      file_.print("\n");
    }

    template<std::forward_iterator Iter>
    void write_cell_data(Iter begin, Iter end, size_t elements_amount) {
      fmt::print(file_, "CELLS {} {}\n", elements_amount, bounded_vertices(begin, end));
      std::for_each(begin, end, [&](const auto& value) {
        file_.print("{} {}", value.size(), fmt::join(value, " "));
        file_.print("\n");
      });
      file_.print("\n");
    }

    template<std::forward_iterator Iter>
    void write_cell_types_data(Iter begin, Iter end, size_t elements_amount) {
      file_.print("CELL_TYPES {}\n", elements_amount);
      file_.print("{}", fmt::join(begin, end, "\n"));
      file_.print("\n");
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
