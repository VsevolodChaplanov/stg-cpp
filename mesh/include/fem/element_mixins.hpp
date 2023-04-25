#ifndef STG_ELEMENT_MIXINS_HPP
#define STG_ELEMENT_MIXINS_HPP

#include <ranges>
#include <range/v3/all.hpp>
#include "concepts.hpp"

namespace stg::mesh::mixin {
using namespace concepts;

  namespace detail {
    template<typename Derived, template<typename> typename CRTPType>
    class CRTP {
    public:
      Derived& down_cast() { return static_cast<Derived&>(*this); }
      const Derived& down_cast() const { return static_cast<Derived const&>(*this); }
    private:
      CRTP() = default;
      friend CRTPType<Derived>;
    };
  }

  template<HasLumped FiniteElement>
  class IntegrateMixin : public detail::CRTP<FiniteElement, IntegrateMixin> {
  public:
    template<std::ranges::viewable_range Range>
    auto integrate(Range&& values) const {
      static_assert(std::is_floating_point_v<std::ranges::range_value_t<Range>>,
                    "Not floating point value undex range");
      assert(std::ranges::distance(values) == this->basis_functions_n());
      return ranges::inner_product(std::forward<Range>(values),
                                   this->down_cast().lumped_, 0.);
    }
  };


}

#endif //STG_ELEMENT_MIXINS_HPP
