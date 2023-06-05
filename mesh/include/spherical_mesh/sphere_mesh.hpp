#ifndef STG_SPHERE_MESH_HPP
#define STG_SPHERE_MESH_HPP

#include <cmath>
#include <concepts>
#include <execution>
#include <geometry/geometry.hpp>
#include <iterator>
#include <numbers>
#include <numeric>
#include <range/v3/numeric/accumulate.hpp>
#include <range/v3/range/concepts.hpp>
#include <range/v3/view/iota.hpp>
#include <range/v3/view/transform.hpp>
#include <range/v3/view/zip.hpp>
#include <thread>
#include <vector>

namespace stg::mesh {
    namespace rv = ranges::views;


    template<std::floating_point T>
    class SphericalMesh final {
    public:
        using value_type = T;

        SphericalMesh(T r, std::size_t n_theta, std::size_t n_phi)
            : r_{r}, n_theta_{n_theta}, n_phi_{n_phi},
              elements_square_routine_thread_{
                      new std::thread(&SphericalMesh::elements_squares, this, r_, n_theta_, n_phi_)} {}

        template<std::forward_iterator Iter>
        value_type integrate(Iter begin, Iter end) const {
            if (elements_square_routine_thread_->joinable()) {
                elements_square_routine_thread_->join();
            }
            auto result = std::inner_product(
                    begin, end, elements_squares_.cbegin(),
                    static_cast<value_type>(0));
            return result;
        }

        template<ranges::viewable_range Range>
        value_type integrate(Range&& values) const {
            if (elements_square_routine_thread_->joinable()) {
                elements_square_routine_thread_->join();
            }
            auto paired = ranges::views::zip(std::forward<Range>(values), elements_squares_) | ranges::views::transform([](const auto& pair) {
                              const auto [value, square] = pair;
                              return value * square;
                          });
            auto result = ranges::accumulate(paired, 0.);
            return result;
        }

        [[nodiscard]] constexpr size_t n_elements() const { return (n_theta_ - 1) * (n_phi_ - 1); }

        std::vector<Point<value_type>> elements_centers() const {
            std::vector<Point<value_type>> k_vectors;

            for (const std::size_t itheta: rv::iota(0ul, n_theta_ - 1)) {
                const double theta_val = itheta * dtheta_;
                const double theta_center_val = theta_val + dtheta_ / 2.;

                for (const std::size_t iphi: rv::iota(0ul, n_phi_ - 1)) {
                    const double phi_val = iphi * dphi_;
                    const double phi_center_val = phi_val + dphi_ / 2.;

                    k_vectors.emplace_back(
                            r_ * std::sin(theta_center_val) * std::cos(phi_center_val),
                            r_ * std::sin(theta_center_val) * std::sin(phi_center_val),
                            r_ * std::cos(theta_center_val));
                }
            }

            return k_vectors;
        }

        ~SphericalMesh() {
            if (elements_square_routine_thread_->joinable()) {
                elements_square_routine_thread_->join();
            }
        }

    private:
        const value_type r_;
        const std::size_t n_theta_;
        const std::size_t n_phi_;
        const value_type dtheta_ = std::numbers::pi_v<value_type> / n_theta_;
        const value_type dphi_ = std::numbers::pi_v<value_type> / n_phi_;
        std::unique_ptr<std::thread> elements_square_routine_thread_;
        std::vector<double> elements_squares_;

        void elements_squares(value_type r, std::size_t n_theta, std::size_t n_phi) {
            elements_squares_.reserve((n_theta_ - 1) * (n_phi_ - 1));

            for (const std::size_t itheta: rv::iota(0ul, n_theta_ - 1)) {
                const double theta_val = itheta * dtheta_;
                const double theta_center_val = theta_val + dtheta_ / 2.;

                for (const std::size_t iphi: rv::iota(0ul, n_phi_ - 1)) {
                    const double phi_val = iphi * dphi_;
                    const double phi_center_val = phi_val + dphi_ / 2.;

                    if (itheta == 0 || itheta == n_theta - 1) {
                        const Point<value_type> point_a = {
                                r_ * std::sin(theta_center_val + dtheta_) * std::cos(phi_center_val),
                                r_ * std::sin(theta_center_val + dtheta_) * std::sin(phi_center_val),
                                r_ * std::cos(theta_center_val + dtheta_)};
                        const Point<value_type> point_b = {
                                r_ * std::sin(theta_center_val + dtheta_) * std::cos(phi_center_val + dphi_),
                                r_ * std::sin(theta_center_val + dtheta_) * std::sin(phi_center_val + dphi_),
                                r_ * std::cos(theta_center_val + dtheta_)};
                        const Point<value_type> point_c = {
                                r_ * std::sin(theta_center_val) * std::cos(phi_center_val),
                                r_ * std::sin(theta_center_val) * std::sin(phi_center_val),
                                r_ * std::cos(theta_center_val)};
                        const auto triangle_sq = triangle_square({point_a, point_b, point_c});
                        elements_squares_.push_back(triangle_sq);
                    } else {
                        const Point<value_type> point_a = {
                                r_ * std::sin(theta_center_val + dtheta_) * std::cos(phi_center_val),
                                r_ * std::sin(theta_center_val + dtheta_) * std::sin(phi_center_val),
                                r_ * std::cos(theta_center_val + dtheta_)};
                        const Point<value_type> point_b = {
                                r_ * std::sin(theta_center_val + dtheta_) * std::cos(phi_center_val + dphi_),
                                r_ * std::sin(theta_center_val + dtheta_) * std::sin(phi_center_val + dphi_),
                                r_ * std::cos(theta_center_val + dtheta_)};
                        const Point<value_type> point_c = {
                                r_ * std::sin(theta_center_val) * std::cos(phi_center_val + dphi_),
                                r_ * std::sin(theta_center_val) * std::sin(phi_center_val + dphi_),
                                r_ * std::cos(theta_center_val)};
                        const Point<value_type> point_d = {
                                r_ * std::sin(theta_center_val) * std::cos(phi_center_val),
                                r_ * std::sin(theta_center_val) * std::sin(phi_center_val),
                                r_ * std::cos(theta_center_val)};
                        const auto element_sq = element_square({point_a, point_b, point_c, point_d});
                        elements_squares_.push_back(element_sq);
                    }
                }
            }
        }

        value_type element_square(const std::vector<Point<T>>& vertices) const {
            const auto& a = vertices[0];
            const auto& b = vertices[1];
            const auto& c = vertices[2];
            const auto& d = vertices[3];

            const auto result_1 = triangle_square({a, b, c});
            const auto result_2 = triangle_square({b, c, d});

            return result_1 + result_2;
        }

        value_type triangle_square(const std::vector<Point<T>>& vertices) const {
            const auto& a = vertices[0];
            const auto& b = vertices[1];
            const auto& c = vertices[2];

            const auto cross_prod_res = cross_product(c - a, c - b);
            auto square = std::sqrt(dot_product(cross_prod_res, cross_prod_res));

            return square;
        }
    };
}// namespace stg::mesh

#endif//STG_SPHERE_MESH_HPP
