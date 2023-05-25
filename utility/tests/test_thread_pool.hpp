#include "common.hpp"
#include "stg_thread_pool.hpp"
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <future>
#include <range/v3/all.hpp>
#include <vector>

struct ThreadPoolTestsFixture {
    std::size_t amount_;
    stg::utility::ThreadPool pool_{4};
    std::vector<std::size_t> values_{amount_};
};


SCENARIO_METHOD(ThreadPoolTestsFixture, "Use thread pool to manage tasks") {

    for (const std::size_t i: ranges::views::iota(0ull, amount_)) {
        auto function = [i, this]() {
            values_[i] = i;
        };
        pool_.post(function);
    }

    for (const std::size_t i: ranges::views::iota(0ull, amount_)) {
        CHECK(i == values_[i]);
    }
}