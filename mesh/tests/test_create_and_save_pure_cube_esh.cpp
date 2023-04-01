#include "common.hpp"

struct VtkSaverFixture {

  static const inline double eps = 1.e-6;

  const double l = 2;
  const std::size_t n = 20;

  const std::size_t nvert = n * n * n;

  const std::string filename = fmt::format("tests_mesh/simple_cube_mesh_{}_{}_{}.vtk", n, n, n);
};


SCENARIO_METHOD(VtkSaverFixture, "Generate pure cube mesh and save as vtk with values") {
  GIVEN("Cubic mesh") {
    CubeMeshBuilder<double> builder{l, n};
    const auto test_mesh = builder.build();

    THEN("Check assembled cube mesh") {
      CHECK(test_mesh->n_elements() == (n - 1) * (n - 1) * (n - 1));
      CHECK(test_mesh->n_vertices() == n * n * n);
      CHECK(test_mesh->relation_table()->n() == n);
    }

    VtkRectilinearGridSaver saver{filename};
    saver.save_mesh<double>(test_mesh->relation_table());
    size_t file_size = std::filesystem::file_size(filename);
    CHECK(file_size > 0);

    std::vector<double> iota_scalars(nvert);
    ranges::iota(iota_scalars, 0);
    saver.save_scalar_data(iota_scalars.cbegin(), iota_scalars.cend());
    size_t new_file_size = std::filesystem::file_size(filename);
    CHECK(new_file_size > file_size);
    file_size = new_file_size;

    std::vector<stg::Vector<double>> iota_vectors;
    iota_vectors.resize(nvert);
    std::for_each(iota_vectors.begin(), iota_vectors.end(), [](auto& value) { value = {1, 1, 1}; });
    saver.save_vector_data(iota_vectors.cbegin(), iota_vectors.cend());
    new_file_size = std::filesystem::file_size(filename);
    CHECK(new_file_size > file_size);
    file_size = new_file_size;

    std::vector<Tensor<double>> iota_tensors;
    iota_tensors.resize(nvert);
    std::for_each(iota_tensors.begin(), iota_tensors.end(), [](auto& value) {
    std::array<double, 9> tensor{1, 0, 0, 0, 1, 0, 0, 0, 1};
    value = Tensor<double>{tensor};
    });
    saver.save_tensor_data(iota_tensors.cbegin(), iota_tensors.cend());
    new_file_size = std::filesystem::file_size(filename);
    CHECK(new_file_size > file_size);
    file_size = new_file_size;
  }
}