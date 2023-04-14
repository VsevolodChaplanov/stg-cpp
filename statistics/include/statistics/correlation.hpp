#ifndef STG_CORRELATION_HPP
#define STG_CORRELATION_HPP

#include <iterator>
#include <stg_tensor/tensor.hpp>
#include "concepts.hpp"
#include "covariance.hpp"
#include "standard_deviation.hpp"

namespace stg::statistics {
using namespace stg::tensor;

  class Correlation {
    template<typename Iter>
    using IteratorValueType = typename std::iterator_traits<Iter>::value_type;
  public:

    /*
     * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
     * * * * * * * * Ranges c++20 stl approach * * * * * * * *
     * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
     */

    template<NumericViewable Range>
    static Tensor<RangeValueType<Range>> correlation_tensor(Range&& f_range, Range&& s_range, Range&& t_range) {
      auto f_mean_f = std::async(std::launch::async, [&] { return Mean::mean(f_range); });
      auto s_mean_f = std::async(std::launch::async, [&] { return Mean::mean(s_range); });
      auto t_mean_f = std::async(std::launch::async, [&] { return Mean::mean(t_range); });

      const auto f_mean = f_mean_f.get();
      const auto s_mean = s_mean_f.get();
      const auto t_mean = t_mean_f.get();

      auto f_std_f = std::async(std::launch::async, [&] { return StandardDeviation::std(f_range, f_mean); });
      auto s_std_f = std::async(std::launch::async, [&] { return StandardDeviation::std(s_range, s_mean); });
      auto t_std_f = std::async(std::launch::async, [&] { return StandardDeviation::std(t_range, t_mean); });

      return correlation_tensor(std::forward<Range>(f_range),
                                std::forward<Range>(s_range),
                                std::forward<Range>(t_range),
                                f_mean, s_mean, t_mean,
                                f_std_f.get(), s_std_f.get(), t_std_f.get());
    }

    template<NumericViewable Range, std::floating_point MeanValueType>
    static auto correlation_tensor(Range&& f_range, Range&& s_range, Range&& t_range,
                                   MeanValueType f_mean, MeanValueType s_mean, MeanValueType t_mean) {
      auto f_std_f = std::async(std::launch::async, [&] { return StandardDeviation::std(f_range, f_mean); });
      auto s_std_f = std::async(std::launch::async, [&] { return StandardDeviation::std(s_range, s_mean); });
      auto t_std_f = std::async(std::launch::async, [&] { return StandardDeviation::std(t_range, t_mean); });

      return correlation_tensor(std::forward<Range>(f_range),
                                std::forward<Range>(s_range),
                                std::forward<Range>(t_range),
                                f_mean, s_mean, t_mean,
                                f_std_f.get(), s_std_f.get(), t_std_f.get());
    }

    template<NumericViewable Range, std::floating_point MeanValueType>
    static auto correlation_tensor(Range&& f_range, Range&& s_range, Range&& t_range,
                                   MeanValueType first_mean, MeanValueType second_mean, MeanValueType third_mean,
                                   MeanValueType first_std, MeanValueType second_std, MeanValueType third_std) {
      auto c_11_f = std::async(std::launch::async, [&] {
        return Correlation::correlation(f_range, f_range,
                                        first_mean, first_mean,
                                        first_std, first_std); });
      auto c_12_f = std::async(std::launch::async, [&] {
        return Correlation::correlation(f_range, s_range,
                                        first_mean, second_mean,
                                        first_std, second_std); });
      auto c_13_f = std::async(std::launch::async, [&] {
        return Correlation::correlation(f_range, t_range,
                                        first_mean, third_mean,
                                        first_std, third_std); });
      auto c_22_f = std::async(std::launch::async, [&] {
        return Correlation::correlation(s_range, s_range,
                                        second_mean, second_mean,
                                        second_std, second_std); });
      auto c_23_f = std::async(std::launch::async, [&] {
        return Correlation::correlation(s_range, t_range,
                                        second_mean, third_mean,
                                        second_std, third_std); });
      auto c_33_f = std::async(std::launch::async, [&] {
        return Correlation::correlation(t_range, t_range,
                                        third_mean, third_mean,
                                        third_std, third_std); });

      auto c_11 = c_11_f.get();
      auto c_12 = c_12_f.get();
      auto c_13 = c_13_f.get();
      auto c_22 = c_22_f.get();
      auto c_23 = c_23_f.get();
      auto c_33 = c_33_f.get();

      return Tensor{ c_11, c_12, c_13,
                     c_12, c_22, c_23,
                     c_13, c_23, c_33};
    }

    template<NumericViewable FirstRange, NumericViewable SecondRange>
    static auto correlation(FirstRange&& f_range, SecondRange&& s_range) {
      const auto f_mean = Mean::mean(f_range);
      const auto s_mean = Mean::mean(s_range);
      const auto f_std = StandardDeviation::std(f_range, f_mean);
      const auto s_std = StandardDeviation::std(s_range, s_mean);

      return correlation(std::forward<FirstRange>(f_range),
                         std::forward<FirstRange>(s_range),
                         f_mean, s_mean, f_std, s_std);
    }

    template<NumericViewable FirstRange, NumericViewable SecondRange, std::floating_point MeanValueType>
    static auto correlation(FirstRange&& f_range, SecondRange&& s_range,
                            MeanValueType f_mean, MeanValueType s_mean) {
      const auto f_std = StandardDeviation::std(f_range, f_mean);
      const auto s_std = StandardDeviation::std(s_range, s_mean);

      return correlation(std::forward<FirstRange>(f_range),
                         std::forward<SecondRange>(s_range),
                         f_mean, s_mean, f_std, s_std);
    }


    template<NumericViewable FirstRange, NumericViewable SecondRange, std::floating_point MeanValueType>
    static auto correlation(FirstRange&& f_range, SecondRange&& s_range,
                            MeanValueType first_mean, MeanValueType second_mean,
                            MeanValueType first_std, MeanValueType second_std) {
      const std::size_t f_size = ranges::distance(f_range);
      const std::size_t s_size = ranges::distance(s_range);
      if (f_size != s_size) {
        throw std::logic_error("Ranges sizes are not equal");
      }

      const auto covariance = Covariance::covariance(std::forward<FirstRange>(f_range),
                                                     std::forward<SecondRange>(s_range),
                                                     first_mean, second_mean);
      auto result = covariance / (first_std * second_std);
      return result;
    }

    /*
     * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
     * * * * * * * * Iterator old stl approach * * * * * * * *
     * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
     */

    template<std::forward_iterator Iter>
    static Tensor<IteratorValueType<Iter>> correlation_tensor(Iter first_begin, Iter first_end,
                                                              Iter second_begin, Iter second_end,
                                                              Iter third_begin, Iter third_end,
                                                              IteratorValueType<Iter> first_mean,
                                                              IteratorValueType<Iter> second_mean,
                                                              IteratorValueType<Iter> third_mean,
                                                              IteratorValueType<Iter> first_std,
                                                              IteratorValueType<Iter> second_std,
                                                              IteratorValueType<Iter> third_std) {
      const auto c_11 = correlation(first_begin, first_end,
                                    first_begin, first_end,
                                    first_mean, first_mean,
                                    first_std, first_std);
      const auto c_12 = correlation(first_begin, first_end,
                                    second_begin, second_end,
                                    first_mean, second_mean,
                                    first_std, second_std);
      const auto c_13 = correlation(first_begin, first_end,
                                    third_begin, third_end,
                                    first_mean, third_mean,
                                    first_std, third_std);
      const auto c_22 = correlation(second_begin, second_end,
                                    second_begin, second_end,
                                    second_mean, second_mean,
                                    second_std, second_std);
      const auto c_23 = correlation(second_begin, second_end,
                                    third_begin, third_end,
                                    second_mean, third_mean,
                                    second_std, third_std);
      const auto c_33 = correlation(third_begin, third_end,
                                    third_begin, third_end,
                                    third_mean, third_mean,
                                    third_std, third_std);

      return {c_11, c_12, c_13, c_12, c_22, c_23, c_13, c_23, c_33};
    }

    template<std::forward_iterator Iter>
    static Tensor<IteratorValueType<Iter>> correlation_tensor(Iter first_begin, Iter first_end,
                                                              Iter second_begin, Iter second_end,
                                                              Iter third_begin, Iter third_end) {
      const auto first_mean = Mean::mean(first_begin, first_end);
      const auto second_mean = Mean::mean(second_begin, second_end);
      const auto third_mean = Mean::mean(third_begin, third_end);
      const auto first_std = StandardDeviation::std(first_begin, first_end, first_mean);
      const auto second_std = StandardDeviation::std(second_begin, second_end, second_mean);
      const auto third_std = StandardDeviation::std(third_begin, third_end, third_mean);

      auto result = correlation_tensor(first_begin, first_end,
                                       second_begin, second_end,
                                       third_begin, third_end,
                                       first_mean, second_mean,
                                       third_mean, first_std,
                                       second_std, third_std);
      return result;
    }

    template<std::forward_iterator Iter>
    static IteratorValueType<Iter> correlation(Iter first_begin, Iter first_end,
                                               Iter second_begin, Iter second_end) {
      const auto first_mean = Mean::mean(first_begin, first_end);
      const auto second_mean = Mean::mean(second_begin, second_end);
      const auto first_std = StandardDeviation::std(first_begin, first_end, first_mean);
      const auto second_std = StandardDeviation::std(second_begin, second_end, second_mean);
      const auto covariance = Covariance::covariance(first_begin, first_end,
                                                     second_begin, second_end,
                                                     first_mean, second_mean);
      auto result = covariance / (first_std * second_std);
      return result;
    }

    template<std::forward_iterator Iter>
    static IteratorValueType<Iter> correlation(Iter first_begin, Iter first_end,
                                               Iter second_begin, Iter second_end,
                                               IteratorValueType<Iter> first_mean,
                                               IteratorValueType<Iter> second_mean,
                                               IteratorValueType<Iter> first_std,
                                               IteratorValueType<Iter> second_std) {
      const auto covariance = Covariance::covariance(first_begin, first_end,
                                                     second_begin, second_end,
                                                     first_mean, second_mean);
      auto result = covariance / (first_std * second_std);
      return result;
    }
  };
}

#endif //STG_CORRELATION_HPP
