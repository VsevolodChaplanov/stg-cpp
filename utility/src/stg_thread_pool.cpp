#include "stg_thread_pool.hpp"
#include <algorithm>
#include <bits/ranges_algo.h>
#include <boost/asio/io_context.hpp>
#include <cstddef>
#include <ranges>
#include <thread>
#include <tuple>
#include <vector>

namespace stg::utility {

    ThreadPool::ThreadPool(std::size_t threads)
        : thread_amount_{threads} {
        std::ranges::for_each(std::views::iota(0ull, threads),
                              [this](std::size_t) {
                                  workers_.emplace_back([this] { io_.run(); });
                              });
    }

    net::io_context& ThreadPool::get_context() {
        return io_;
    }
}// namespace stg::utility