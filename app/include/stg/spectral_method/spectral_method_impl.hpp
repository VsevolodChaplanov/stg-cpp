#ifndef STG_SPECTRAL_METHOD_IMPL_HPP
#define STG_SPECTRAL_METHOD_IMPL_HPP

#include "data_loader.hpp"
#include <boost/asio/post.hpp>
#include <boost/asio/thread_pool.hpp>
#include <concepts>
#include <fem.hpp>
#include <memory>
#include <mesh_builders.hpp>
#include <range/v3/view/iota.hpp>
#include <statistics.hpp>
#include <stg_generators.hpp>
#include <stg_tensor/tensor.hpp>
#include <velocity_field.hpp>

namespace stg::spectral {
    using namespace stg::mesh;
    using namespace stg::field;
    using namespace stg::generators;
    using namespace stg::statistics;
    namespace net = boost::asio;
    namespace rv = ranges::views;

    template<std::floating_point T, std::size_t seed = 42>
    class SpectralMethodApplicationImpl final {
    public:
        struct Statistics {
            /*
       * Assuming that the mean of fluctuations is 0 for physical reason
       */
            double mean_x = 0.;
            double mean_y = 0.;
            double mean_z = 0.;
            double std_x;
            double std_y;
            double std_z;
        };

        using value_type = T;

        SpectralMethodApplicationImpl(T cube_edge_len, std::size_t edge_points,
                                      SpectralGeneratorConfig<T> config,
                                      std::size_t samples_amount = 100,
                                      std::size_t concurrency_hint = std::thread::hardware_concurrency() - 1)
            : SpectralMethodApplicationImpl(CubeMeshBuilder<value_type>{cube_edge_len, edge_points},
                                            std::forward<SpectralGeneratorConfig<value_type>>(config),
                                            samples_amount,
                                            concurrency_hint) {}

        SpectralMethodApplicationImpl(CubeMeshBuilder<T> builder,
                                      SpectralGeneratorConfig<T> config,
                                      std::size_t samples_amount = 100,
                                      std::size_t concurrency_hint = std::thread::hardware_concurrency() - 1)
            : SpectralMethodApplicationImpl(
                      builder.build(),
                      std::forward<SpectralGeneratorConfig<value_type>>(config),
                      samples_amount,
                      concurrency_hint) {}

        SpectralMethodApplicationImpl(std::shared_ptr<const CubeFiniteElementsMesh<T>> mesh,
                                      SpectralGeneratorConfig<T> config,
                                      std::size_t samples_amount = 100,
                                      std::size_t concurrency_hint = std::thread::hardware_concurrency() - 1)
            : fe_mesh_{std::move(mesh)}, generator_config_(std::move(config)), spectral_generator_{generator_config_}, velocity_samples_{100, fe_mesh_->n_vertices()}, thread_pool_{new net::thread_pool{std::max(1ul, concurrency_hint)}} {}

        /*
     * Use generator to directly generate value
     */
        Vector<value_type> generate_at_point(const Point<value_type>& point, value_type time) const {
            auto result = spectral_generator_(point, time);
            return result;
        }

        /*
     * Generate velocity fluctuations in mesh vertices
     */
        void generate_on_mesh(value_type time) {
            auto executor_ctx = thread_pool_->get_executor();
            for (const std::size_t index: rv::iota(0ul, fe_mesh_->n_vertices())) {
                net::post(*thread_pool_.get(),
                          [this, index, time] { generate_at_vertex(index, time); });
            }
            thread_pool_->join();
        }

        /*
     * Generate samples for statistics in mesh vertices
     */
        void generate_samples_on_mesh(value_type time) {
            for (const std::size_t index: rv::iota(0ul, velocity_samples_.size())) {
                generate_sample(index, time);
                /*net::post(*thread_pool_.get(),
                  [this, index, time]{ generate_sample(index, time); });*/
            }
            thread_pool_->join();
        }

        Statistics collect_ansamble_statistics() const {
            if (cache_.is_ansamble_cache_) { return cache_.ansamble_cache_; }

            const auto vx_view = velocity_field_.vx_view();
            const auto vy_view = velocity_field_.vy_view();
            const auto vz_view = velocity_field_.vz_view();

            const auto std_x_f = std::async(std::launch::async, &Mean::mean, vx_view, 0.);
            const auto std_y_f = std::async(std::launch::async, &Mean::mean, vy_view, 0.);
            const auto std_z_f = std::async(std::launch::async, &Mean::mean, vz_view, 0.);

            cache_.ansamble_cache_ = Statistics{
                    .std_x = std_x_f.get(),
                    .std_y = std_y_f.get(),
                    .std_z = std_z_f.get()};
            cache_.is_ansamble_cache_ = true;

            return cache_.ansamble_cache_;
        }

        /*
     * Save velocity, tensors and all scalar data to vtk file
     */
        void save_to_vtk(std::string_view filename,
                         std::string_view divergences_table_name = "Divergences",
                         std::string_view velocity_table_name = "VelocityField",
                         std::string_view corr_tensors_table_name = "CorrelationTensors") const {
            VtkRectilinearGridSaver saver{filename};

            saver.save_mesh(fe_mesh_->relation_table());
            if (!divergences_.empty())
                saver.save_scalar_data(divergences_.cbegin(),
                                       divergences_.cend(),
                                       divergences_table_name);
            if (!velocity_field_.empty()) {
                auto val_view = velocity_field_.values_view();
                saver.save_velocity_data(val_view, velocity_table_name);
            }
            if (!corr_tensor_data_.empty())
                saver.save_tensor_data(corr_tensor_data_.cbegin(),
                                       corr_tensor_data_.cend(),
                                       corr_tensors_table_name);
        }

        /*
     * Generate correlations tensors relative to vertex with index i,j,k
     */
        void generate_correlations(std::size_t ix, std::size_t jy, std::size_t kz) {
            corr_tensor_data_.resize(fe_mesh_->n_vertices());
            const std::size_t base_vert_index = fe_mesh_->relation_table()->lin_index(ix, jy, kz);

            for (const std::size_t ivert: rv::iota(0ul, fe_mesh_->n_vertices())) {
                auto space_corr = [&]() mutable {
                    auto base_x_sample = velocity_samples_.vx_component_for_vertex(base_vert_index);
                    auto base_y_sample = velocity_samples_.vy_component_for_vertex(base_vert_index);
                    auto base_z_sample = velocity_samples_.vz_component_for_vertex(base_vert_index);

                    auto vert_x_sample = velocity_samples_.vx_component_for_vertex(ivert);
                    auto vert_y_sample = velocity_samples_.vy_component_for_vertex(ivert);
                    auto vert_z_sample = velocity_samples_.vz_component_for_vertex(ivert);

                    corr_tensor_data_[ivert] = SpaceCorrelation::correlation_tensor(base_x_sample, vert_x_sample,
                                                                                    base_y_sample, vert_y_sample,
                                                                                    base_z_sample, vert_z_sample,
                                                                                    0., 0., 0., 0., 0., 0.);
                };
                space_corr();
                /*net::post(*thread_pool_.get(), space_corr);*/
            }

            thread_pool_->join();
        }

        void generate_correlations_relative_to_space_center() {
            const auto center_ind = fe_mesh_->center_tri_index();
            generate_correlations(center_ind[0], center_ind[1], center_ind[2]);
        }

        ~SpectralMethodApplicationImpl() {
            thread_pool_->join();
        }

    private:
        const std::shared_ptr<const CubeFiniteElementsMesh<value_type>> fe_mesh_;
        const SpectralGeneratorConfig<value_type> generator_config_;
        const SpectralGenerator<value_type, seed> spectral_generator_;
        VelocityField<value_type> velocity_field_{fe_mesh_->n_vertices()};
        VelocitySamples<value_type> velocity_samples_{100, fe_mesh_->n_vertices()};
        std::vector<tensor::Tensor<value_type>> corr_tensor_data_;
        std::vector<value_type> divergences_;
        std::unique_ptr<net::thread_pool> thread_pool_;

        struct {
            Statistics ansamble_cache_;
            bool is_ansamble_cache_ = false;
        } cache_;

        void generate_at_vertex(std::size_t ivert, value_type time) {
            const auto vert_point = fe_mesh_->relation_table()->vertex(ivert);
            const auto velocity = spectral_generator_(vert_point, time);
            velocity_field_.set_value(velocity, ivert);
        }

        void generate_sample(std::size_t isample, value_type time) {
            VelocityField<value_type> sample{fe_mesh_->n_vertices()};
            for (const std::size_t index: rv::iota(0ul, fe_mesh_->n_vertices())) {
                auto&& vertex = fe_mesh_->relation_table()->vertex(index);
                auto&& value = spectral_generator_(std::move(vertex), time);
                sample.set_value(std::move(value), index);
            }
            velocity_samples_.set_sample(std::move(sample), isample);
        }
    };


    namespace fs = std::filesystem;

    template<std::floating_point T>
    class SpectralMethodApplication final {
    public:
        using value_type = T;

        SpectralMethodApplication(DataLoader loader, SpectralParameters<value_type> params, std::shared_ptr<ISpectra<value_type>> spectra)
            : SpectralMethodApplication(std::move(loader), params,
                                        std::make_shared<SpectralGeneratorV2<value_type>>(params),
                                        CubeMeshBuilder<value_type>{params.cube_edge_len, params.edge_points},
                                        std::move(spectra)) {}

        SpectralMethodApplication(DataLoader loader, SpectralParameters<value_type> params,
                                  std::shared_ptr<SpectralGeneratorV2<value_type>> generator,
                                  const CubeMeshBuilder<value_type>& builder,
                                  std::shared_ptr<ISpectra<value_type>> spectra)
            : loader_{std::move(loader)}, parameters_{std::move(params)}, fe_mesh_{builder.build()}, spectral_generator_{std::move(generator)} {
            spectral_generator_->initialize_spectra(spectra);
            spectral_generator_->initialize_wave_vector_amplitudes(parameters_.k_min, parameters_.k_max, parameters_.n_spectra);
            spectral_generator_->initialize_random_coefficients();
            spectral_generator_->initialize_inner_generators();

            velocity_field_.resize(fe_mesh_->n_vertices());
        }

        void generate_velocity_field(value_type time) {
            auto func = [time, this](std::size_t g_index) {
                const auto vertex = fe_mesh_->relation_table()->vertex(g_index);
                velocity_field_.set_value(spectral_generator_->operator()(vertex, time), g_index);
            };
            auto executor = pool_.get_executor();
            for (const std::size_t g_index: rv::iota(0ull, fe_mesh_->n_vertices())) {
                // func(g_index);
                net::post(executor, std::bind(func, g_index));
            }
            pool_.join();
        }

        value_type get_max_period() const {
            auto max_period = spectral_generator_->max_period();
            return max_period;
        }

        void save_data_to(std::filesystem::path path, std::string_view table_name = "Vector field") {
            VtkRectilinearGridSaver saver{path.string()};
            saver.template save_mesh<value_type>(fe_mesh_->relation_table());
            saver.save_vector_data(velocity_field_.values_view(), table_name);
        }

        void set_spectra(std::shared_ptr<ISpectra<value_type>> spectra) {
            spectral_generator_->initialize_spectra(std::move(spectra));
        }

    private:
        DataLoader loader_;
        SpectralParameters<value_type> parameters_;
        const std::shared_ptr<const CubeFiniteElementsMesh<value_type>> fe_mesh_ = CubeMeshBuilder<value_type>{parameters_.cube_edge_len, parameters_.edge_points}.build();
        const std::shared_ptr<SpectralGeneratorV2<value_type>> spectral_generator_ = std::make_shared<SpectralGeneratorV2<value_type>>(parameters_);
        VelocityField<value_type> velocity_field_;


        net::thread_pool pool_{std::thread::hardware_concurrency()};
    };
}// namespace stg::spectral

#endif//STG_SPECTRAL_METHOD_IMPL_HPP
