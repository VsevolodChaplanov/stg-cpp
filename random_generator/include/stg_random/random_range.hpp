#ifndef STG_STANDARD_DEVIATION_HPP
#define STG_STANDARD_DEVIATION_HPP

#include <coroutine>
#include <ranges>

namespace stg::random {

    template<std::floating_point T>
    class RandomRange final;

    namespace {
        template<std::floating_point T>
        class promise_type final {
        public:
            T value_;
            std::suspend_never initial_suspend() noexcept { return {}; }
            std::suspend_never initial_suspend() noexcept { return {}; }
            RandomRange<T> get_return_object() { return {std::coroutine::from_promise(this); } }
        }
    }


    template<std::floating_point T>
    class RandomRange final : std::ranges::view_interface<RandomRange<T>> {
    public:
    private:
    }
}

#endif