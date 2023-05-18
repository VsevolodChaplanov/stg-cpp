#include "common.hpp"
#include <fmt/format.h>
#include <fmt/printf.h>
#include <fmt/ostream.h>
#include <fmt/core.h>
#include <fmt/os.h>
#include <fmt/std.h>
#include <iomanip>
#include <fstream>

using namespace stg::spectral;

struct SpectralMethodApplicationTestFixture {
  const SpectralParameters<double> parameters {
    .seed = 42,
    .ampl_mean = 0, .ampl_std = 1,
    .wv_mean = 0, .wv_std = 1,
    .freq_mean = 0, .freq_std = 1. / 2.,
    .cube_edge_len = 10.,
    .k_min = 0.1, .k_max = 100,
    .length_scale = 1.e-3, .time_scale = 1.e-3,
    .n_spectra = 1000, .n_fourier = 1000,
    .edge_points = 21,
    .save_data_dir_path = "./spectral_result/"
  };

  const std::string directory_with_files = "./spectral_result/";

  stg::spectral::DataLoader loader{directory_with_files};
};

SCENARIO_METHOD(SpectralMethodApplicationTestFixture, "Generate samples along time");