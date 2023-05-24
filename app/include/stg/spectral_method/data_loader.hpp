#ifndef STG_SPECTRAL_DATA_LOADER_HPP
#define STG_SPECTRAL_DATA_LOADER_HPP

#include <filesystem>
#include <vtk_parser/rectilinear_grid_parser.hpp>
#include <velocity_field/velocity_field.hpp>
#include <velocity_field/velocity_samples.hpp>
#include <velocity_field/velocity_samples_1d.hpp>
#include <mesh_builders.hpp>
#include <statistics.hpp>

namespace stg::spectral {
namespace fs = std::filesystem;
namespace rv = ranges::views;
using namespace stg::mesh;
using namespace stg::tensor;
using namespace stg::field;

  class DataLoader final {
  public:
    explicit DataLoader(fs::path directory_with_files)
      : work_dir_{directory_with_files} { }

    template<std::floating_point T>
    std::shared_ptr<CubeFiniteElementsMesh<T>> load_mesh(std::string_view file_with_grid) const {
      std::stringstream sstream;
      auto wd = work_dir_;
      wd.append(file_with_grid);
      const std::string grid_filename{sstream.str()};
      RectilinearGridParser parser{wd.string()};

      return parser.fe_mesh<T>();
    }

    template<std::floating_point T>
    std::vector<T> load_scalar_data(std::string_view filename) const {
      auto path_to_file = work_dir_;
      path_to_file.append(filename);
      RectilinearGridParser parser{path_to_file.string()};
      return parser.scalar_data<T>();
    }

    template<std::floating_point T>
    std::vector<Tensor<T>> load_covariation_data() const {
      RectilinearGridParser parser_xx{work_dir_.string() + "r_11.vtk"};
      auto xx = parser_xx.scalar_data<T>();
      RectilinearGridParser parser_xy{work_dir_.string() + "r_12.vtk"};
      auto xy = parser_xy.scalar_data<T>();
      RectilinearGridParser parser_xz{work_dir_.string() + "r_13.vtk"};
      auto xz = parser_xz.scalar_data<T>();
      RectilinearGridParser parser_yy{work_dir_.string() + "r_22.vtk"};
      auto yy = parser_yy.scalar_data<T>();
      RectilinearGridParser parser_yz{work_dir_.string() + "r_23.vtk"};
      auto yz = parser_yz.scalar_data<T>();
      RectilinearGridParser parser_zz{work_dir_.string() + "r_33.vtk"};
      auto zz = parser_zz.scalar_data<T>();

      const std::size_t data_size = ranges::size(xx);
      std::vector<Tensor<T>> cov_data(data_size);
      for (const std::size_t ivert: rv::iota(0ull, data_size)) {
        cov_data[ivert] = Tensor{ xx[ivert], xy[ivert], xz[ivert],
                                  xy[ivert], yy[ivert], yz[ivert],
                                  xz[ivert], yz[ivert], zz[ivert]
        };
      }

      return cov_data;
    }

    template<std::floating_point T>
    VelocitySamples<T> load_velocity_samples() const {
      std::vector<VelocityField<T>> samples;
      samples.reserve(100);
      for (const auto & entry : fs::directory_iterator(work_dir_)) {
        const bool is_velocity_field_file = entry.path().filename().string().starts_with("sg3");
        if (is_velocity_field_file) {
          RectilinearGridParser parser{entry.path().string()};
          const auto file_field = parser.vector_data<T>();
          std::vector<T> vx, vy, vz;
          for (const auto [vx_v, vy_v, vz_v] : rv::transform(file_field, [](const auto& val) {
            return std::make_tuple(val.template get<0>(), val.template get<1>(), val.template get<2>());
          })) {
            vx.emplace_back(vx_v);
            vy.emplace_back(vy_v);
            vz.emplace_back(vz_v);
          }
          samples.emplace_back(std::move(vx), std::move(vy), std::move(vz));
        }
      }
      return VelocitySamples<T>{std::move(samples)};
    }

    template<std::floating_point T>
    VelocitySamples1D<T> load_velocity_samples_1sg() const {
      std::vector<VelocityField1D<T>> samples;
      samples.reserve(100);
      for (const auto & entry : fs::directory_iterator(work_dir_)) {
        const bool is_velocity_field_file = entry.path().filename().string().starts_with("sg1");
        if (is_velocity_field_file) {
          RectilinearGridParser parser{entry.path().string()};
          auto file_field = parser.scalar_data<T>();
          samples.emplace_back(std::move(file_field));
        }
      }
      return VelocitySamples1D<T>{std::move(samples)};
    }

    template<std::floating_point T>
    std::vector<Tensor<T>> load_fert() const {
      RectilinearGridParser parser_xx{work_dir_.string() + "phi_11.vtk"};
      auto xx = parser_xx.scalar_data<T>();
      RectilinearGridParser parser_xy{work_dir_.string() + "phi_12.vtk"};
      auto xy = parser_xy.scalar_data<T>();
      RectilinearGridParser parser_xz{work_dir_.string() + "phi_13.vtk"};
      auto xz = parser_xz.scalar_data<T>();
      RectilinearGridParser parser_yy{work_dir_.string() + "phi_22.vtk"};
      auto yy = parser_yy.scalar_data<T>();
      RectilinearGridParser parser_yz{work_dir_.string() + "phi_23.vtk"};
      auto yz = parser_yz.scalar_data<T>();
      RectilinearGridParser parser_zz{work_dir_.string() + "phi_33.vtk"};
      auto zz = parser_zz.scalar_data<T>();

      const std::size_t data_size = ranges::size(xx);
      std::vector<Tensor<T>> cov_data(data_size);
      for (const std::size_t ivert: rv::iota(0ull, data_size)) {
        cov_data[ivert] = Tensor{ xx[ivert], xy[ivert], xz[ivert],
                                  xy[ivert], yy[ivert], yz[ivert],
                                  xz[ivert], yz[ivert], zz[ivert]
        };
      }

      return cov_data;
    }

  private:
    const fs::path work_dir_;
  };
}

#endif //STG_DATA_LOADER_HPP
