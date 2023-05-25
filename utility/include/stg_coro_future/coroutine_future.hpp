#ifndef STG_UTILITY_FUTURE_COROUTINE_HPP
#define STG_UTILITY_FUTURE_COROUTINE_HPP

#include <coroutine>
#include <exception>
#include <future>
#include <thread>
#include <type_traits>

/*

    Here provides mechanism of lazy async evalutaions

    usage:

        std::future<EvalResult> do_evalutation() {
            auto first_value = co_await std::async([] () { return SomeLongConstuctObject; })
            auto first_value = co_await std::async([] () { return SomeOtherLongConstuctObject; })
            auto main_result = co_await std::async([] (SomeLongConstuctObject a, SomeOtherLongConstuctObject b) {
                return a + b;
            })
            co_return main_result
        }

*/

struct as_coroutine {};

namespace {

    template<typename T>
    concept FutureTypeConcept = !std::is_void_v<T> && !std::is_reference_v<T>;
}// namespace

namespace std {
    template<FutureTypeConcept T, typename... Args>
    struct std::coroutine_traits<std::future<T>, as_coroutine, Args...> {
        struct promise_type : std::promise<T> {
            std::future<T> get_return_object() { return this->get_future(); }
            std::suspend_never initial_suspend() noexcept { return {}; }
            std::suspend_never final_suspend() noexcept { return {}; }
            void return_value(T&& value) noexcept(std::is_nothrow_move_constructible_v<T>) {
                this->set_value(std::forward<T>(value));
            }
            void return_value(const T& value) noexcept(std::is_nothrow_copy_constructible_v<T>) {
                this->set_value(std::forward<T>(value));
            }
            void unhandled_exception() noexcept(std::is_nothrow_copy_constructible_v<decltype(std::current_exception())>) {
                this->set_exception(std::current_exception());
            }
        };
    };

    // std::future<void>
    template<typename... Args>
    struct std::coroutine_traits<std::future<void>, as_coroutine, Args...> {
        struct promise_type : std::promise<void> {
            std::future<void> get_return_object() { return this->get_future(); }
            std::suspend_never initial_suspend() noexcept { return {}; }
            std::suspend_never final_suspend() noexcept { return {}; }
            void unhandled_exception() noexcept(std::is_nothrow_copy_constructible_v<decltype(std::current_exception())>) {
                this->set_exception(std::current_exception());
            }
        };
    };
}// namespace std

namespace stg::utils {
    // future awaiter
    template<typename T>
    struct future_awaiter final : public std::future<T> {
        bool await_ready() noexcept {
            using namespace std::chrono_literals;
            return this->wait_for(0s) != std::future_status::timeout;
        }
        T await_resume() noexcept(std::is_nothrow_move_constructible_v<T> || std::is_nothrow_copy_constructible_v<T>) {
            return this->get();
        }
        void await_suspend(std::coroutine_handle<> coro_handle) {
            std::thread([this, coro_handle] {
                this->wait();
                coro_handle();
            }).detach();
        }
    };
}// namespace stg::utils

template<typename T>
    requires(!std::is_reference_v<T>)
auto operator co_await(std::future<T> future) { return stg::utils::future_awaiter<T>{std::move(future)}; }

#endif