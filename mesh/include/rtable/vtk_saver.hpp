#ifndef STG_VTK_SAVER_HPP
#define STG_VTK_SAVER_HPP

#include <execution>
#include <string_view>
#include <filesystem>
#include <fstream>
#include <type_traits>
#include <iterator>
#include <iostream>
#include <fmt/core.h>
#include <fmt/ostream.h>
#include <geometry/geometry.hpp>
#include <stg_tensor/tensor.hpp>
#include "i_relation_table.hpp"
#include "cube_relation_table.hpp"


namespace stg::mesh {

  namespace detail {
    template<typename T, template<typename...> typename ReferenceType>
    struct is_specialization : std::false_type {};

    template<template<typename...> typename ReferenceType, typename... Args>
    struct is_specialization<ReferenceType<Args...>, ReferenceType> : std::true_type {};

    template<typename T, template<typename...> typename ReferenceType>
    constexpr bool is_specialization_v = is_specialization<T, ReferenceType>::value;
  }


  class VtkSaver {
  public:
    VtkSaver() = default;

    template<std::floating_point T>
    void save_mesh(const std::shared_ptr<IRelationTable<T>>& rtable, std::string_view filename) const {
      const auto& vertices = rtable->vertices();
      const auto& bound_indices = rtable->elements_global_indices();
      const auto& elements_types = rtable->elements_types();
      std::filesystem::path path = std::filesystem::absolute(filename);
      if (!std::filesystem::exists(path)) {
        std::filesystem::create_directories(path.parent_path());
      }
      std::ofstream file{path, std::ios_base::out};
      write_header(file);
      write_point_data(file, vertices.cbegin(), vertices.cend(), rtable->n_vertices());
      write_cell_data(file, bound_indices.cbegin(), bound_indices.cend(), rtable->n_elements());
      write_cell_types_data(file, elements_types.cbegin(), elements_types.cend(), rtable->n_elements());
    }

    template<std::forward_iterator Iter>
    void save_scalar_data(std::string_view filename,
                          Iter begin, Iter end,
                          std::string_view table_name = "DefaultTable") const {
      using IterValueType = typename std::iterator_traits<Iter>::value_type;
      std::filesystem::path path = std::filesystem::absolute(filename);
      if (!std::filesystem::exists(path)) {
        std::filesystem::create_directories(path.parent_path());
      }
      std::ofstream file{path, std::ios_base::app};
      if (!has_point_data_flag_) {
        fmt::print(file, "POINT_DATA {}\n", std::distance(begin, end));
        has_point_data_flag_ = true;
      }
      fmt::print(file, "SCALARS {} double\n", table_name);
      fmt::print(file, "LOOKUP_TABLE default\n");
      std::copy(begin, end, std::ostream_iterator<IterValueType>{file, "\n"});
      file << std::endl;
    }

    template<std::forward_iterator Iter>
    void save_vector_data(std::string_view filename,
                          Iter begin, Iter end,
                          std::string_view table_name = "VectorField") const {
      std::filesystem::path path = std::filesystem::absolute(filename);
      if (!std::filesystem::exists(path)) {
        std::filesystem::create_directories(path.parent_path());
      }
      std::ofstream file{path, std::ios_base::app};
      if (!has_point_data_flag_) {
        fmt::print(file, "POINT_DATA {}\n", std::distance(begin, end));
        has_point_data_flag_ = true;
      }
      fmt::print(file, "VECTORS {} double\n", table_name);
      std::for_each(begin, end, [&file](const auto& vector) {
        fmt::print(file, "{} {} {}\n",
                   vector.template get<0>(),
                   vector.template get<1>(),
                   vector.template get<2>());
        file << std::endl;
      });
      file << std::endl;
    }

    template<std::forward_iterator Iter>
    void save_tensor_data(std::string_view filename,
                          Iter begin, Iter end,
                          std::string_view table_name = "TensorData") const {
      using IterValueType = typename std::iterator_traits<Iter>::value_type;
      using TensorValueType = typename IterValueType::value_type;
      std::filesystem::path path = std::filesystem::absolute(filename);
      if (!std::filesystem::exists(path)) {
        std::filesystem::create_directories(path.parent_path());
      }
      std::ofstream file{path, std::ios_base::app};
      if (!has_point_data_flag_) {
        fmt::print(file, "POINT_DATA {}\n", std::distance(begin, end));
        has_point_data_flag_ = true;
      }
      fmt::print(file, "TENSORS {} double\n", table_name);
      std::for_each(begin, end, [&file](const auto& tensor) {
        std::copy(tensor.cbegin(), tensor.cend(), std::ostream_iterator<TensorValueType>{file, " "});
        file << std::endl;
      });
      file << std::endl;
    }

  protected:
    void write_header(std::ofstream& file) const {
      fmt::print(file,
        "# vtk DataFile Version 3.0\n"
        "RFG method\n"
        "ASCII\n"
        "DATASET UNSTRUCTURED_GRID\n");
      file << std::endl;
    }

    template<std::forward_iterator Iter>
    void write_point_data(std::ofstream& file, Iter begin, Iter end, size_t points_number) const {
      fmt::print(file, "POINTS {} double\n", points_number);
      std::for_each(begin, end, [&file](const auto& value) {
        static size_t counter = 0;
        file << value << " ";
        if (++counter % 3 == 0) {
          file << std::endl;
        }
      });
      file << std::endl;
    }

    template<std::forward_iterator Iter>
    void write_cell_data(std::ofstream& file, Iter begin, Iter end, size_t elements_amount) const {
      fmt::print(file, "CELLS {} {}\n", elements_amount, bounded_vertices(begin, end));
      std::for_each(begin, end, [&file](const auto& value) {
        file << value.size() << " ";
        std::copy(value.cbegin(), value.cend(), std::ostream_iterator<size_t>{file, " "});
        file << std::endl;
      });
      file << std::endl;
    }

    template<std::forward_iterator Iter>
    void write_cell_types_data(std::ofstream& file, Iter begin, Iter end, size_t elements_amount) const {
      fmt::print(file, "CELL_TYPES {}\n", elements_amount);
      std::copy(begin, end, std::ostream_iterator<size_t>{file, "\n"});
      file << std::endl;
    }

    template<std::forward_iterator Iter>
    size_t bounded_vertices(Iter begin, Iter end) const {
      auto result = std::transform_reduce(std::execution::par_unseq, begin, end, 0ul,
        [](auto left, auto right) {
          return left + right + 1ul;
      },
        [](const auto& elem) {
          return elem.size();
      });
      return result;
    }

    mutable bool has_point_data_flag_ = false;
  };
}

#endif //STG_VTK_SAVER_HPP
