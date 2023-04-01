#include "common.hpp"
#include <fmt/format.h>
#include <range/v3/experimental/utility/generator.hpp>

struct CubePrizmFEMeshTestFixture {
  static const inline double eps = 1.e-6;

  const double x_l = -1.5;
  const double x_r = 1.5;
  const double y_l = -1.5;
  const double y_r = 1.5;
  const double z_l = -1.5;
  const double z_r = 1.5;
  const std::size_t nx = 50;
  const std::size_t ny = 50;
  const std::size_t nz = 50;

  const CubePrizmMeshBuilder<double> builder { x_l, x_r, y_l, y_r, z_l, z_r, nx, ny, nz };
  const std::shared_ptr<CubePrizmFEMesh<double>> test_mesh = builder.build();
  const VtkSaver saver;
};

SCENARIO_METHOD(CubePrizmFEMeshTestFixture, "Create cube mesh") {
  GIVEN("Assembled cube finite elements mesh") {
    THEN("Check assembled cube mesh") {
      CHECK(test_mesh->n_elements() == (nx - 1) * (ny - 1) * (nz - 1));
      CHECK(test_mesh->n_vertices() == nx * ny * nz);
      CHECK(test_mesh->nx() == nx);
      CHECK(test_mesh->ny() == ny);
      CHECK(test_mesh->nz() == nz);
    }

    THEN("Integrate one to get volume of domain") {
      std::vector<double> ones(test_mesh->n_vertices(), 1.);

      const auto integral = test_mesh->integrate(ones);
      CHECK_THAT(integral, Catch::Matchers::WithinRel((x_r - x_l) * (y_r - y_l) * (z_r - z_l), eps));
    }

    THEN("Integrate real function values") {
      std::vector<double> function(test_mesh->n_vertices());
      for (const size_t index: ranges::views::iota(0ul, test_mesh->n_vertices())) {
        const auto vertex = test_mesh->relation_table()->vertex(index);
        function[index] = vertex.get<0>() * vertex.get<1>() * vertex.get<2>();
      }

      const auto integral = test_mesh->integrate(function);
      CHECK_THAT(integral, Catch::Matchers::WithinAbs(0., eps));
    }

    THEN("Integrate real function values 2") {
      std::vector<double> function(test_mesh->n_vertices());
      for (const size_t index: ranges::views::iota(0ul, test_mesh->n_vertices())) {
        const auto vertex = test_mesh->relation_table()->vertex(index);
        function[index] = vertex.get<0>() * vertex.get<0>() * vertex.get<1>() * vertex.get<2>();  // x**2 * y * z
      }

      const auto integral = test_mesh->integrate(function);
      CHECK_THAT(integral, Catch::Matchers::WithinAbs(0, eps));

      const std::string filename = fmt::format("tests_mesh/cube_mesh_tests_scenario_1_{}_{}_{}.vtk", nx, ny, nz);

      saver.save_mesh(test_mesh->relation_table(), filename);
      saver.save_scalar_data(filename,
                             function.cbegin(), function.cend());
    }
  }
}


struct NonSymmetricCubePrizmFEMeshTestFixture {
  static const inline double eps = 1.e-6;

  const double x_l = -0.5;
  const double x_r = 0.5;
  const double y_l = -1.;
  const double y_r = 1.;
  const double z_l = -1.5;
  const double z_r = 1.5;
  const std::size_t nx = 30;
  const std::size_t ny = 40;
  const std::size_t nz = 50;

  const CubePrizmMeshBuilder<double> builder { x_l, x_r, y_l, y_r, z_l, z_r, nx, ny, nz };
  const std::shared_ptr<CubePrizmFEMesh<double>> test_mesh = builder.build();
  const VtkSaver saver;
};

SCENARIO_METHOD(NonSymmetricCubePrizmFEMeshTestFixture, "Create cube mesh with non symm parameters") {
  GIVEN("Assembled cube finite elements mesh") {
    THEN("Check assembled cube mesh") {
      CHECK(test_mesh->n_elements() == (nx - 1) * (ny - 1) * (nz - 1));
      CHECK(test_mesh->n_vertices() == nx * ny * nz);
      CHECK(test_mesh->nx() == nx);
      CHECK(test_mesh->ny() == ny);
      CHECK(test_mesh->nz() == nz);
    }

    THEN("Integrate one to get volume of domain") {
      std::vector<double> ones(test_mesh->n_vertices(), 1.);

      const auto integral = test_mesh->integrate(ones);
      CHECK_THAT(integral, Catch::Matchers::WithinRel((x_r - x_l) * (y_r - y_l) * (z_r - z_l), eps));
    }

    THEN("Integrate real function values") {
      std::vector<double> function(test_mesh->n_vertices());
      for (const size_t index: ranges::views::iota(0ul, test_mesh->n_vertices())) {
        const auto vertex = test_mesh->relation_table()->vertex(index);
        function[index] = vertex.get<0>() * vertex.get<1>() * vertex.get<2>();
      }

      const auto integral = test_mesh->integrate(function);
      CHECK_THAT(integral, Catch::Matchers::WithinAbs(0., eps));
    }

    THEN("Integrate real function values 2") {
      std::vector<double> function(test_mesh->n_vertices());
      for (const size_t index: ranges::views::iota(0ul, test_mesh->n_vertices())) {
        const auto vertex = test_mesh->relation_table()->vertex(index);
        function[index] = vertex.get<0>() * vertex.get<0>() * vertex.get<1>() * vertex.get<2>();  // x**2 * y * z
      }

      const auto integral = test_mesh->integrate(function);
      CHECK_THAT(integral, Catch::Matchers::WithinAbs(0, eps));

      const std::string filename = fmt::format("tests_mesh/cube_mesh_tests_scenario_2_{}_{}_{}.vtk", nx, ny, nz);

      saver.save_mesh(test_mesh->relation_table(), filename);
      saver.save_scalar_data(filename,
                             function.cbegin(), function.cend());
    }
  }
}



struct NonSymmetricAndNonSymmAlongDomainCubePrizmFEMeshTestFixture {
  static const inline double eps = 1.e-4;

  const double x_l = 0.5;
  const double x_r = 1.;
  const double y_l = -0.5;
  const double y_r = 0.1;
  const double z_l = 0.5;
  const double z_r = 1.5;
  const std::size_t nx = 50;
  const std::size_t ny = 60;
  const std::size_t nz = 70;

  const CubePrizmMeshBuilder<double> builder { x_l, x_r, y_l, y_r, z_l, z_r, nx, ny, nz };
  const std::shared_ptr<CubePrizmFEMesh<double>> test_mesh = builder.build();
  const VtkSaver saver;
};

SCENARIO_METHOD(NonSymmetricAndNonSymmAlongDomainCubePrizmFEMeshTestFixture,
                "Create cube mesh with non symm parameters with asymm in domain params") {
  GIVEN("Assembled cube finite elements mesh") {
    THEN("Check assembled cube mesh") {
      CHECK(test_mesh->n_elements() == (nx - 1) * (ny - 1) * (nz - 1));
      CHECK(test_mesh->n_vertices() == nx * ny * nz);
      CHECK(test_mesh->nx() == nx);
      CHECK(test_mesh->ny() == ny);
      CHECK(test_mesh->nz() == nz);
    }

    THEN("Integrate one to get volume of domain") {
      std::vector<double> ones(test_mesh->n_vertices(), 1.);

      const auto integral = test_mesh->integrate(ones);
      CHECK_THAT(integral, Catch::Matchers::WithinRel((x_r - x_l) * (y_r - y_l) * (z_r - z_l), eps));
    }

    THEN("Integrate real function values") {
      std::vector<double> function(test_mesh->n_vertices());
      for (const size_t index: ranges::views::iota(0ul, test_mesh->n_vertices())) {
        const auto vertex = test_mesh->relation_table()->vertex(index);
        function[index] = vertex.get<0>() * vertex.get<1>() * vertex.get<2>();
      }

      const auto integral = test_mesh->integrate(function);
      CHECK_THAT(integral, Catch::Matchers::WithinRel(-0.045, eps));
    }

    THEN("Integrate real function values 2") {
      std::vector<double> function(test_mesh->n_vertices());
      for (const size_t index: ranges::views::iota(0ul, test_mesh->n_vertices())) {
        const auto vertex = test_mesh->relation_table()->vertex(index);
        function[index] = vertex.get<0>() * vertex.get<0>() * vertex.get<1>() * vertex.get<2>();  // x**2 * y * z
      }

      const auto integral = test_mesh->integrate(function);
      CHECK_THAT(integral, Catch::Matchers::WithinRel(-0.035, eps));

      const std::string filename = fmt::format("tests_mesh/cube_mesh_tests_scenario_3_{}_{}_{}.vtk", nx, ny, nz);

      saver.save_mesh(test_mesh->relation_table(), filename);
      saver.save_scalar_data(filename,
                             function.cbegin(), function.cend());
    }
  }
}


struct CubePrizmMeshIndicesChecksFixture {
  static const inline double eps = 1.e-6;

  const double x_l = -1.5;
  const double x_r = 1.5;
  const double y_l = -1.5;
  const double y_r = 1.5;
  const double z_l = -1.5;
  const double z_r = 1.5;
  const std::size_t nx = 5;
  const std::size_t ny = 6;
  const std::size_t nz = 7;

  const CubePrizmMeshBuilder<double> builder { x_l, x_r, y_l, y_r, z_l, z_r, nx, ny, nz };
  const std::shared_ptr<CubePrizmFEMesh<double>> test_mesh = builder.build();
  const VtkSaver saver;
};

SCENARIO_METHOD(CubePrizmMeshIndicesChecksFixture, "Indices conversion check") {
  GIVEN("Assembled cube finite elements mesh") {
    THEN("Check assembled cube mesh") {
      CHECK(test_mesh->n_elements() == (nx - 1) * (ny - 1) * (nz - 1));
      CHECK(test_mesh->n_vertices() == nx * ny * nz);
      CHECK(test_mesh->nx() == nx);
      CHECK(test_mesh->ny() == ny);
      CHECK(test_mesh->nz() == nz);

      const auto base_index = test_mesh->tri_index(0);
      CHECK(base_index[0] == 0);
      CHECK(base_index[1] == 0);
      CHECK(base_index[2] == 0);
      const auto on_x = test_mesh->tri_index(2);
      CHECK(on_x[0] == 2);
      CHECK(on_x[1] == 0);
      CHECK(on_x[2] == 0);
      const auto on_xy = test_mesh->tri_index(6);
      CHECK(on_xy[0] == 1);
      CHECK(on_xy[1] == 1);
      CHECK(on_xy[2] == 0);
      const auto on_z = test_mesh->tri_index(30);
      CHECK(on_z[0] == 0);
      CHECK(on_z[1] == 0);
      CHECK(on_z[2] == 1);
      const auto on_xyz = test_mesh->tri_index(36);
      CHECK(on_xyz[0] == 1);
      CHECK(on_xyz[1] == 1);
      CHECK(on_xyz[2] == 1);

      const auto base_lin_index = test_mesh->lin_index(0, 0, 0);
      const auto on_x_lin_index = test_mesh->lin_index(2, 0, 0);
      const auto on_xy_lin_index = test_mesh->lin_index(1, 1, 0);
      const auto on_z_lin_index = test_mesh->lin_index(0, 0, 1);
      const auto on_xyz_lin_index = test_mesh->lin_index(1, 1, 1);
      CHECK(base_lin_index == 0);
      CHECK(on_x_lin_index == 2);
      CHECK(on_xy_lin_index == 6);
      CHECK(on_z_lin_index == 30);
      CHECK(on_xyz_lin_index == 36);
    }
  }
}