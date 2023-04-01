#include "common.hpp"

struct ParsingFixture {
  static inline const double eps = 1.e-6;
  const std::string scalar_file_filename
    = "test_resources/phi_11.vtk";
};

SCENARIO_METHOD(ParsingFixture, "Parse file with scalar data") {
  RectilinearGridParser parser{scalar_file_filename};

  const auto result_mesh = parser.mesh<double>();

  THEN("Check assembled cube mesh") {
    CHECK(result_mesh->n_elements() == 50 * 50 * 50);
    CHECK(result_mesh->n_vertices() == 51 * 51 * 51);
  }

  const auto scalar_data = parser.scalar_data<double>();

  REQUIRE(scalar_data.size() == 132651);
  CHECK_THAT(scalar_data[0], WithinRel(3.69811e-07, eps));
  CHECK_THAT(scalar_data[2], WithinRel(4.29192e-07, eps));
  CHECK_THAT(scalar_data.back(), WithinRel(3.69811e-07, eps));
}

struct ParsingVectorDataFixture {
  static inline const double eps = 1.e-6;
  const std::string scalar_file_filename
    = "test_resources/sg3_L10_N21_eigcut200_try0.vtk";
};

SCENARIO_METHOD(ParsingVectorDataFixture, "Parse file with vector data") {
  RectilinearGridParser parser{scalar_file_filename};

  const auto result_mesh = parser.mesh<double>();

  THEN("Check assembled cube mesh") {
    CHECK(result_mesh->n_elements() == 20 * 20 * 20);
    CHECK(result_mesh->n_vertices() == 21 * 21 * 21);
    CHECK_THAT(result_mesh->vertex(0).get<0>(), WithinRel(-4.7619, eps));
    CHECK_THAT(result_mesh->vertex(0).get<1>(), WithinRel(-4.7619, eps));
    CHECK_THAT(result_mesh->vertex(0).get<2>(), WithinRel(-4.7619, eps));
  }

  const auto vector_data = parser.vector_data<double>();

  REQUIRE(vector_data.size() == 9261);
  CHECK_THAT(vector_data[0].get<0>(), WithinRel(-0.195596, eps));
  CHECK_THAT(vector_data[0].get<1>(), WithinRel(-0.00458737, eps));
  CHECK_THAT(vector_data[0].get<2>(), WithinRel(-0.00960132, eps));

  CHECK_THAT(vector_data[2].get<0>(), WithinRel(-0.218068, eps));
  CHECK_THAT(vector_data[2].get<1>(), WithinRel(0.0717094, eps));
  CHECK_THAT(vector_data[2].get<2>(), WithinRel(0.0634021, eps));

  CHECK_THAT(vector_data.back().get<0>(), WithinRel(0.0268544, eps));
  CHECK_THAT(vector_data.back().get<1>(), WithinRel(-0.214446, eps));
  CHECK_THAT(vector_data.back().get<2>(), WithinRel(-0.0763451, eps));
}