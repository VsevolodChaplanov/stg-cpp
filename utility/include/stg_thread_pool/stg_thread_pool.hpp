#ifndef STG_UTILITY_STG_THREAD_POOL_HPP
#define STG_UTILITY_STG_THREAD_POOL_HPP

#include <algorithm>
#include <boost/asio/execution_context.hpp>
#include <boost/asio/io_context.hpp>
#include <boost/asio/post.hpp>
#include <boost/asio/thread_pool.hpp>
#include <concepts>
#include <cstddef>
#include <thread>
#include <type_traits>
#include <vector>

namespace stg::utility {
    namespace net = boost::asio;
    class ThreadPool final {
    public:
        explicit ThreadPool(std::size_t threads = std::max(1u, std::thread::hardware_concurrency()));

        ThreadPool(const ThreadPool&) = delete;
        ThreadPool& operator=(const ThreadPool&) = delete;

        ThreadPool(ThreadPool&&) = delete;
        ThreadPool& operator=(ThreadPool&&) = delete;

        net::io_context& get_context();

        template<typename Task>
        void post(Task&& task) {
            net::post(io_, std::forward<Task>(task));
        }

    private:
        const std::size_t thread_amount_;
        net::io_context io_{static_cast<int>(thread_amount_)};
        std::vector<std::jthread> workers_{thread_amount_};
    };
}// namespace stg::utility

#endif