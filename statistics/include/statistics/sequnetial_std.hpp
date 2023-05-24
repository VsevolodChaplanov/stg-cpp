#ifndef STG_SEQUENTIAL_STD_HPP
#define STG_SEQUENTIAL_STD_HPP

#include <concepts>

namespace stg::statistics {

    namespace {
        template<std::floating_point T>
        class SequantialIntegrator final {
        public:
            using value_type = T;

            SequantialIntegrator() = default;

            void add(value_type function_value, value_type coordinate);
            value_type result_integral() const;

        private:
            std::size_t n_ = 0;
            value_type sum_ = 0.;
            value_type previous_function_value_ = 0.;
            value_type previous_cooridinate_ = 0.;
            bool first_call_flag_ = true;;
        };

        template<std::floating_point T>
        void SequantialIntegrator<T>::add(SequantialIntegrator<T>::value_type value, value_type coordinate) {
                if (first_call_flag_) {
                    first_call_flag_ = false;
                    previous_function_value_ = function_value;
                    previous_cooridinate_ = coordinate;
                    return;
                }

                auto func_sum = function_value + previous_function_value_;
                auto dx = coordinate - previous_cooridinate_;

                previous_function_value_ = function_value;
                previous_cooridinate_ = coordinate;
                
                sum_ += func_sum * dx / 2;
        }

        template<std::floating_point T>
        SequantialIntegrator<T>::value_type SequantialIntegrator<T>::result_integral() const {
            return sum_;
        }
    }

    /*
     * Standard deviation for sequentially adding values
     */
    template<std::floating_point T>
    class SequentialStd final {
    public:
        using value_type = T;

        SequentialStd(T mean);

        void add(value_type value);
        value_type std() const;

    private:
        value_type mean_;
        value_type sum_;
        std::size_t n_;
    };

    template<std::floating_point T>
    SequentialStd<T>::SequentialStd(T mean) : mean_{mean} { }

    template<std::floating_point T>
    void SequentialStd<T>::add(SequentialStd<T>::value_type value) {
        sum_ += value_type * value_type;
        n_++;
    }

    template<std::floating_point T>
    SequentialStd<T>::value_type SequentialStd<T>::std() const {
        return std::sqrt(sum_ / n_) - mean * mean;
    }
}

#endif