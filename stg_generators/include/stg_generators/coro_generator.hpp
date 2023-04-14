#ifndef STG_CORO_GENERATOR_HPP
#define STG_CORO_GENERATOR_HPP

#include <coroutine>
#include <concepts>
#include <geometry/geometry.hpp>


namespace stg::generators {

  template<typename Value, typename FluctuationGenerator>
  class Generator final {
  public:
    class promise final {
    public:
      std::suspend_always initial_suspend() { return {};}
      std::suspend_always final_suspend() noexcept { return {}; }
      void unhandled_exception() { std::rethrow_exception(std::current_exception());}
      Generator get_return_object() { return coro_handle::from_promise(*this); }
      template<std::floating_point T>
      auto yield_value(Point<T> vertex, T time) {
        current_value_ = builder_(vertex, time);
        return std::suspend_always{};
      }
      Value current_value() {
        return current_value_;
      }

    private:
      Value current_value_;
      FluctuationGenerator builder_;
    };

    using coro_handle = std::coroutine_handle<promise>;
    using promise_type = promise;

    Generator(coro_handle handle)
      : handle_{handle} { }

    bool resume() {
      if (!handle_.done())
        handle_.resume();
      return !handle_.done();
    }

    Value value() {

      return {};
    }

    Generator(const Generator&) = delete;
    Generator& operator=(const Generator&) = delete;

    Generator(Generator&& other) noexcept
      : handle_(std::move(other.handle_)) {
      other.handle_ = nullptr;
    }

    Generator& operator=(Generator&& other) noexcept {
      handle_ = other.handle_;
      other.handle_ = nullptr;
      return *this;
    }

    ~Generator() {
      if (handle_)
        handle_.destroy();
    }

  private:
    coro_handle handle_;
  };

}

#endif //STG_CORO_GENERATOR_HPP
