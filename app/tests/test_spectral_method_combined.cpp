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
  const SpectralParameters<double> parameters { };

  const std::string directory_with_files = "./spectral_result/";
  static constexpr inline auto path_template = "./spectral_result/field_t{}.vtk";

  stg::spectral::DataLoader loader{directory_with_files};

  const double t_start = 0;
  const std::size_t n_time_shots = 1000;
};

SCENARIO_METHOD(SpectralMethodApplicationTestFixture, "Generate samples along time") {
  SpectralMethodApplication<double> application{loader, parameters, std::make_shared<KolmogorovSpectra<double>>()};

  const auto t_end = 100 * parameters.time_scale;
  const auto tau = (t_end - t_start) / (n_time_shots - 1ull);
  double time_moment = t_start;
  for (const std::size_t time_index : std::views::iota(0ull, n_time_shots)) {
    const std::string table_name = fmt::format("VelocityField{}", time_index);
    application.generate_velocity_field(time_moment);
    application.save_data_to(fmt::format(path_template, time_index), table_name);
    time_moment += tau;
  }
}