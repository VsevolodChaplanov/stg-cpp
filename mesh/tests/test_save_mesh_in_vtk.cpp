#include <filesystem>
#include <string>
#include <array>
#include <fmt/core.h>
#include <fmt/format.h>
#include "common.hpp"


struct VtkSaverFixture {

  static const inline double eps = 1.e-6;

  const double x_l = -1.;
  const double x_r = 1.;
  const double y_l = -1.;
  const double y_r = 1.;
  const double z_l = -1.;
  const double z_r = 1.;
  const std::size_t nx = 10;
  const std::size_t ny = 10;
  const std::size_t nz = 10;

  const std::size_t nvert = nx * ny * nz;

  const std::string filename = fmt::format("tests_mesh/simple_cube_mesh_{}_{}_{}.vtk", nx, ny, nz);
};


SCENARIO_METHOD(VtkSaverFixture, "Generate cube mesh and save as vtk with values") {
  GIVEN("Cubic mesh") {
    CubePrizmMeshBuilder<double> builder { x_l, x_r, y_l, y_r, z_l, z_r, nx, ny, nz };

    const auto test_mesh = builder.build();

    THEN("Check assembled cube mesh") {
      CHECK(test_mesh->n_elements() == (nx - 1) * (ny - 1) * (nz - 1));
      CHECK(test_mesh->n_vertices() == nx * ny * nz);
      CHECK(test_mesh->nx() == nx);
      CHECK(test_mesh->ny() == ny);
      CHECK(test_mesh->nz() == nz);
    }

    VtkSaver saver{};
    saver.save_mesh<double>(test_mesh->relation_table(), filename);
    size_t file_size = std::filesystem::file_size(filename);
    CHECK(file_size > 0);

    std::vector<double> iota_scalars(nvert);
    ranges::iota(iota_scalars, 0);
    saver.save_scalar_data(filename, iota_scalars.cbegin(), iota_scalars.cend());
    size_t new_file_size = std::filesystem::file_size(filename);
    CHECK(new_file_size > file_size);
    file_size = new_file_size;

    std::vector<stg::Vector<double>> iota_vectors;
    iota_vectors.resize(nvert);
    std::for_each(iota_vectors.begin(), iota_vectors.end(), [](auto& value) { value = {1, 1, 1}; });
    saver.save_vector_data(filename, iota_vectors.cbegin(), iota_vectors.cend());
    new_file_size = std::filesystem::file_size(filename);
    CHECK(new_file_size > file_size);
    file_size = new_file_size;

    std::vector<Tensor<double>> iota_tensors;
    iota_tensors.resize(nvert);
    std::for_each(iota_tensors.begin(), iota_tensors.end(), [](auto& value) {
      std::array<double, 9> tensor{1, 0, 0, 0, 1, 0, 0, 0, 1};
      value = Tensor<double>{tensor};
    });
    saver.save_tensor_data(filename, iota_tensors.cbegin(), iota_tensors.cend());
    new_file_size = std::filesystem::file_size(filename);
    CHECK(new_file_size > file_size);
    file_size = new_file_size;
  }
}