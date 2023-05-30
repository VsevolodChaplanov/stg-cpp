#include "common.hpp"
#include "stg/kriging/data_loader.hpp"
#include "stg/kriging/sequential_kriging_1d_analysis.hpp"
#include "stg/spectral_method/sequential_analysis.hpp"
#include <string>

struct SequentialCorrelationsTestsFixture {
    const std::string dir = "./spectral_result/";
    const std::string mesh_file = "velocity_field_0.vtk";
    const std::string save = "./spectral_result/covariations.vtk";
};

SCENARIO_METHOD(SequentialCorrelationsTestsFixture, "Calculate squentially covariation function") {

    stg::spectral::DataLoader loader{dir};
    SequentialCorrelations<double> analyser{std::move(loader), mesh_file};

    analyser.calc_covariations_for_amount(1000);
    analyser.save_calculated_covariations(save);
}