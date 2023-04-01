#include <iostream>
#include <vector>
#include <range/v3/all.hpp>
#include <statistics.hpp>

using namespace stg::statistics;

class VelocityField final {
public:

  auto zip_values_view() const {
    return ranges::view::zip(vx, vy, vz);
  }

  auto vx_view() const {
    return ranges::view::all(vx);
  }

  auto vy_view() const {
    return ranges::view::all(vy);
  }

  auto vz_view() const {
    return vz | ranges::views::all;
//    return ranges::view::all(vz);
  }

protected:
  std::vector<double> vx{0, 1, 2};
  std::vector<double> vy{3, 4, 5};
  std::vector<double> vz{6, 7, 8};
};

class VelocitySamples final {
public:
  VelocitySamples(std::size_t samples_amount)
    : velocity_samples_(samples_amount) {
  }

  VelocitySamples(std::vector<VelocityField> samples)
    : velocity_samples_{std::move(samples)} { }

  auto vx_component_for_vertex(std::size_t ivert) const {
    return velocity_samples_ | ranges::views::transform(
        [ivert](const VelocityField& field) {
          return field.vx_view() | ranges::views::drop(ivert) | ranges::views::take(1);
        }) | ranges::views::join;
  }

  auto vy_component_for_vertex(std::size_t ivert) const {
    return velocity_samples_ | ranges::views::transform(
      [ivert](const VelocityField& field) {
        return field.vy_view() | ranges::views::drop(ivert) | ranges::views::take(1);
      }) | ranges::views::join;
  }

  auto vz_component_for_vertex(std::size_t ivert) const {
    return velocity_samples_ | ranges::views::transform(
      [ivert](const VelocityField& field) {
        return field.vz_view() | ranges::views::drop(ivert) | ranges::views::take(1);
      }) | ranges::views::join;
  }

private:
  std::vector<VelocityField> velocity_samples_;
};

int main() {
  std::vector<VelocityField> sample(5);
  VelocitySamples samples(5);

  for (const auto& val : samples.vx_component_for_vertex(1)) {
    std::cout << val << std::endl;
  }

  auto f_x_range = samples.vx_component_for_vertex(0);
  auto f_y_range = samples.vy_component_for_vertex(0);
  auto f_z_range = samples.vz_component_for_vertex(0);
  auto s_x_range = samples.vx_component_for_vertex(1);
  auto s_y_range = samples.vy_component_for_vertex(1);
  auto s_z_range = samples.vz_component_for_vertex(1);

  const auto res = SpaceCovariance::covariance_tensor(f_x_range, s_x_range,
                                                      f_y_range, s_y_range,
                                                      f_z_range, s_z_range);

  const auto res_2 = SpaceCovariance::covariance_tensor(f_x_range, s_x_range,
                                                      f_y_range, s_y_range,
                                                      f_z_range, s_z_range,
                                                      0., 0., 0., 0., 0., 0.);

  const auto res_3 = SpaceCorrelation::correlation_tensor(f_x_range, s_x_range,
                                                          f_y_range, s_y_range,
                                                          f_z_range, s_z_range);

  return EXIT_SUCCESS;
}