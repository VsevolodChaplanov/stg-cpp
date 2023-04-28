#ifndef STG_FE_ELEMENT_CRTP_HPP
#define STG_FE_ELEMENT_CRTP_HPP

#include <array>
#include <range/v3/all.hpp>
#include <range/v3/view/all.hpp>
#include <range/v3/numeric/inner_product.hpp>
#include "concepts.hpp"

namespace std::mesh {

  namespace detail {
      
    template<typename Derived, template<typename> typename CRTPType>
    class crtp_typed {
    public:
      Derived& underlying() { return static_cast<Derived&>(*this); }
      const Derived& underlying() const { return static_cast<Derived const&>(*this); }
    protected:
      crtp_typed() = default;
      friend CRTPType<Derived>;
    };

    template<typename Derived>
    class crtp {
    public:
      Derived& underlying() { return static_cast<Derived&>(*this); }
      const Derived& underlying() const { return static_cast<Derived const&>(*this); }
    protected:
      crtp() = default;
    };
  } // namespace detail

  namespace matrix {
    template<std::floating_point T, std::size_t np>
    class LumpedMatrix {
    public: 
      using value_type = T;

      constexpr value_type lumped(std::size_t index) const { return this->lumped_matrix_[index]; }
      auto lumped_view() const { return this->lumped_matrix_ | ranges::views::all; }
      
      template<ranges::viewable_range Matrix>
      void fill_lumped(Matrix&& values) {
        ranges::copy(std::forward<Matrix>(values), this->lumped_matrix_);
      }

    protected:
      std::array<value_type, np> lumped_matrix_; 
    };
  }

  namespace mixin {
    template<typename Element>
    class Integrate : public detail::crtp_typed<Element, Integrate> {
    public:
      using value_type = typename Element::value_type;
      
      template<ranges::viewable_range Range>
      constexpr value_type integrate(Range&& values) const {
        static_assert(std::is_floating_point_v<ranges::range_value_t<Range>>, 
                      "Not floating point type under range");
        assert(std::ranges::distance(values) == this->basis_functions_n());
        auto lumped_und = this->underlying().lumped_view();
        auto values_begin = ranges::begin(values);
        auto lumped_begin = std::begin(lumped_und);
        auto lumped_end = std::begin(lumped_und);
        auto result = std::inner_product(lumped_begin, lumped_end, values_begin, 0.);
        return result;
      }
    };
  }

  namespace elements {
    
    template<typename Impl, std::size_t np>
    class LagrangianElement : detail::crtp<Impl> {
    public:
      constexpr std::size_t basis_functions_n() const { return np; }
      constexpr double lumped(std::size_t index) const { return this->underlying().lumped(index); }
    };

    class VoxelElement : public LagrangianElement<VoxelElement, 8> { };
  }

  template<typename ElementType>
  class FiniteElement final : public detail::crtp<ElementType> {

  };

} // namespace std::mesh


#endif