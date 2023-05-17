#ifndef STG_SPECTRA_APPROXIMATOR_HPP
#define STG_SPECTRA_APPROXIMATOR_HPP

#include <concepts>
#include <type_traits>
#include <vector>
#include <range/v3/all.hpp>

namespace stg::generators {
namespace rv = ranges::views;

  template<std::floating_point T>
  class SpectraApproximator {
  public:
    using value_type = T;

    SpectraApproximator(T k_min, T k_max, std::size_t N);
  protected:
    const value_type k_min_;
    const value_type k_max_;
    const std::size_t n_;
  };

  template<std::floating_point T>
  class HuangSpectraApproximator final : SpectraApproximator<T> {
  public:
    using value_type = SpectraApproximator<T>::val;
    HuangSpectraApproximator(double k_min, double k_max, std::size_t N);

  protected:
    std::vector<value_type> k_m_;
  };

  template<std::floating_point T>
  HuangSpectraApproximator<T>::HuangSpectraApproximator(double k_min, double k_max, std::size_t N)
    : SpectraApproximator<T>(k_min, k_max, N) {}
}

#endif //STG_SPECTRA_APPROXIMATOR_HPP
