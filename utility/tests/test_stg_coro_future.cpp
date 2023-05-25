#include "common.hpp"
#include <catch2/catch_test_macros.hpp>
#include <chrono>
#include <fmt/core.h>
#include <fmt/printf.h>
#include <future>
#include <range/v3/numeric/accumulate.hpp>
#include <range/v3/view/iota.hpp>
#include <sstream>
#include <thread>

using namespace std::chrono_literals;

unsigned long long no_async() {
    auto function = []() {
        auto val = ranges::accumulate(ranges::views::iota(0ull, 10'000'000ull), 0);
        return val;
    };
    auto first = function();
    auto second = function();

    return first + second;
}

unsigned long long no_co_await() {
    auto function = []() {
        auto val = ranges::accumulate(ranges::views::iota(0ull, 10'000'000ull), 0);
        return val;
    };
    auto first = std::async(std::launch::async, function);
    auto second = std::async(std::launch::async, function);

    return first.get() + second.get();
}

std::future<unsigned long long> first_test(as_coroutine) {
    auto first_result = co_await std::async(std::launch::async, [] {
        auto val = ranges::accumulate(ranges::views::iota(0ull, 10'000'000ull), 0);
        return val;
    });
    auto second_result = co_await std::async(std::launch::async, [] {
        auto val = ranges::accumulate(ranges::views::iota(0ull, 10'000'000ull), 0);
        return val;
    });
    co_return first_result + second_result;
}

std::future<unsigned long long> second_test(as_coroutine) {
    auto first_result = co_await std::async(std::launch::async, [] {
        auto val = ranges::accumulate(ranges::views::iota(0ull, 10'000'000ull), 0);
        return val;
    });
    auto second_result = co_await std::async(
            std::launch::async, [](auto prev) {
        auto val = ranges::accumulate(ranges::views::iota(0ull, 10'000'000ull), 0);;
        return val + prev; }, first_result);
    co_return second_result;
}

std::future<unsigned long long> third_test(as_coroutine) {
    auto first_result = std::async(std::launch::async, [] {
        auto val = ranges::accumulate(ranges::views::iota(0ull, 10'000'000ull), 0);
        return val;
    });
    auto second_result = co_await std::async(
            std::launch::async, [](auto prev) {
        auto val = ranges::accumulate(ranges::views::iota(0ull, 10'000'000ull), 0);;
        return val + prev; }, first_result.get());
    co_return second_result;
}

SCENARIO("Coroitunes utilites tests") {
    auto exact_result = 2 * ranges::accumulate(ranges::views::iota(0ull, 10'000'000ull), 0);

    {
        auto start = std::chrono::high_resolution_clock::now();
        auto no_async_v = no_async();
        CHECK(no_async_v == exact_result);
        fmt::print("no_async {}\n", (std::chrono::high_resolution_clock::now() - start).count());
    }

    {
        auto start = std::chrono::high_resolution_clock::now();
        auto no_co_await_v = no_co_await();
        CHECK(no_co_await_v == exact_result);
        fmt::print("no co await {}\n", (std::chrono::high_resolution_clock::now() - start).count());
    }

    {
        auto start = std::chrono::high_resolution_clock::now();
        auto first_test_v = first_test({});
        CHECK(first_test_v.get() == exact_result);
        fmt::print("first_test {}\n", (std::chrono::high_resolution_clock::now() - start).count());
    }

    {
        auto start = std::chrono::high_resolution_clock::now();
        auto second_test_v = second_test({});
        CHECK(second_test_v.get() == exact_result);
        fmt::print("second_test {}\n", (std::chrono::high_resolution_clock::now() - start).count());
    }

    {
        auto start = std::chrono::high_resolution_clock::now();
        auto third_test_v = third_test({});
        CHECK(third_test_v.get() == exact_result);
        fmt::print("third_test {}\n", (std::chrono::high_resolution_clock::now() - start).count());
    }
}

// 1000389620
// 499354
// 2000870535