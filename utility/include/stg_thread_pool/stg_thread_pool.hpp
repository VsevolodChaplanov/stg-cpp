#ifndef STG_UTILITY_STG_THREAD_POOL_HPP
#define STG_UTILITY_STG_THREAD_POOL_HPP

#include <algorithm>
#include <boost/asio/execution_context.hpp>
#include <boost/asio/io_context.hpp>
#include <boost/asio/post.hpp>
#include <boost/asio/thread_pool.hpp>
#include <concepts>
#include <ranges>
#include <stop_token>
#include <thread>
#include <type_traits>
#include <vector>

namespace stg::utility {
    namespace net = boost::asio;
    class ThreadPool final {
    public:
        explicit ThreadPool(std::size_t threads = std::max(1u, std::thread::hardware_concurrency()))
            : thread_amount_{threads} {
            // std::ranges::for_each(std::views::iota(0ull, threads),
            //                       [this](std::size_t) {
            //                           workers_.emplace_back([this](std::stop_token) { io_.run(); });
            //                       });
        }

        ThreadPool(const ThreadPool&) = delete;
        ThreadPool& operator=(const ThreadPool&) = delete;

        ThreadPool(ThreadPool&&) = delete;
        ThreadPool& operator=(ThreadPool&&) = delete;

        net::thread_pool& get_context() {
            return thread_pool_;
        }

        template<typename Task>
        void post(Task&& task) {
            net::post(thread_pool_, std::forward<Task>(task));
        }

    private:
        const std::size_t thread_amount_;
        net::thread_pool thread_pool_{thread_amount_};
        // std::vector<std::jthread> workers_{thread_amount_};
    };
}// namespace stg::utility

#endif