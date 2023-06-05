#ifndef _SPACE_HPP
#define _SPACE_HPP

#include <array>
#include <fstream>
#include <string>
#include <utility>
#include <vector>

namespace stg::kriging {
    using point_t = std::array<double, 3>;

    // 3d symmetric cube [-L/2, L/2]
    // divided by N segments in each direction
    struct Space {
        Space(size_t N, double L);

        [[nodiscard]] point_t point(size_t i, size_t j, size_t k) const;
        [[nodiscard]] point_t point(size_t ind) const;

        [[nodiscard]] size_t lin_index(size_t i, size_t j, size_t k) const;
        [[nodiscard]] std::array<size_t, 3> tri_index(size_t ind) const;

        [[nodiscard]] bool point_within(const point_t& p) const;
        [[nodiscard]] double interpolate_at(const point_t& p, const std::vector<double>& f) const;
        [[nodiscard]] std::array<double, 3> interpolate_at3(const point_t& p,
                                                            const std::vector<double>& f1,
                                                            const std::vector<double>& f2,
                                                            const std::vector<double>& f3) const;
        [[nodiscard]] std::array<double, 6> interpolate_at6(const point_t& p,
                                                            const std::vector<double>& f1,
                                                            const std::vector<double>& f2,
                                                            const std::vector<double>& f3,
                                                            const std::vector<double>& f4,
                                                            const std::vector<double>& f5,
                                                            const std::vector<double>& f6) const;

        const size_t N;
        const size_t N3;
        const double L;
        const double A;
        const double B;
        const double h;
        const std::vector<double> coo;

        void tovtk(std::string fn, const std::vector<double>& f) const {
            return tovtk(std::move(fn), f.begin(), f.end());
        }

        template<typename DataIter>
        void tovtk(std::string fn, DataIter begin, DataIter end) const {
            std::ofstream fs(fn);
            tovtk_init(fs);
            for (DataIter it = begin; it != end; ++it) {
                fs << *it << std::endl;
            }
        }

        void tovtk(std::string fn,
                   const std::vector<double>& u,
                   const std::vector<double>& v) const {
            return tovtk(std::move(fn), u.begin(), u.end(), v.begin(), v.end());
        }

        template<typename DataIter>
        void tovtk(std::string fn,
                   DataIter ubegin, DataIter uend,
                   DataIter vbegin, DataIter vend) const {
            std::ofstream fs(fn);
            tovtk_init(fs, 2);
            DataIter uit = ubegin;
            DataIter vit = vbegin;
            while (uit != uend && vit != vend) {
                fs << *uit++ << " " << *vit++ << std::endl;
            }
        }

        void tovtk(std::string fn,
                   const std::vector<double>& u,
                   const std::vector<double>& v,
                   const std::vector<double>& w) const {
            return tovtk(std::move(fn), u.begin(), u.end(), v.begin(), v.end(), w.begin(), w.end());
        }

        template<typename DataIter>
        void tovtk(std::string fn,
                   DataIter ubegin, DataIter uend,
                   DataIter vbegin, DataIter vend,
                   DataIter wbegin, DataIter wend) const {
            std::ofstream fs(fn);
            tovtk_init(fs, 3);
            DataIter uit = ubegin;
            DataIter vit = vbegin;
            DataIter wit = wbegin;
            while (uit != uend && vit != vend && wit != wend) {
                fs << *uit++ << " " << *vit++ << " " << *wit++ << std::endl;
            }
        }

    private:
        void tovtk_init(std::ostream& fs, int datadim = 1) const;
        [[nodiscard]] double use_interpolation_polynom(const std::array<size_t, 8>& indices, const std::array<double, 8>& bases, const std::vector<double>& v) const;
        void interpolation_polynom(const point_t& p, std::array<size_t, 8>& indices, std::array<double, 8>& bases) const;
    };

    struct PhisicalSpace;
    struct FourierSpace;

    struct PhysicalSpace : public Space {
        PhysicalSpace(size_t N, double L);
        [[nodiscard]] FourierSpace fourier_space() const;
    };

    struct FourierSpace : public Space {
        FourierSpace(size_t N, double Lk);
        [[nodiscard]] PhysicalSpace physical_space() const;
    };
}// namespace stg::kriging
#endif
