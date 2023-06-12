#include "common.hpp"

struct AnalysisKrigingDataFixture {
    const std::size_t n = 21;
    stg::spectral::DataLoader data_loader{"test_resources/"};
};

SCENARIO_METHOD(AnalysisKrigingDataFixture, "Do analysis of data obtained by kriging method") {
    GIVEN("Analysis class") {
        KrigingAnalysis<double> kriging_test{data_loader};

        THEN("Check assembled meshes") {
            const auto& real_space = kriging_test.real_space_mesh();
            const auto& fourier_space = kriging_test.fourier_space_mesh();


            CHECK(real_space->n_elements() == (n - 1) * (n - 1) * (n - 1));
            CHECK(real_space->n_vertices() == n * n * n);
            CHECK(real_space->relation_table()->n() == n);

            CHECK(fourier_space->n_elements() == (n - 1) * (n - 1) * (n - 1));
            CHECK(fourier_space->n_vertices() == n * n * n);
            CHECK(fourier_space->relation_table()->n() == n);

            AND_THEN("Calculate covariations") {
                const auto& calc_cov_test = kriging_test.calculate_covariations();
                const auto& cov_data = kriging_test.covariations_from_data();

                VtkRectilinearGridSaver saver{"tests_mesh/kriging_calculated_covariations.vtk"};
                saver.save_mesh(kriging_test.real_space_mesh()->relation_table());
                saver.save_tensor_data(calc_cov_test.cbegin(), calc_cov_test.cend());
            }
        }
    }
}