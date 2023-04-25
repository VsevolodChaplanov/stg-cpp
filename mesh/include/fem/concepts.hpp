#ifndef STG_CONCEPTS_HPP
#define STG_CONCEPTS_HPP

#include <concepts>

namespace stg::mesh::concepts {

  template<typename Element>
  concept HasBasisFunctionAmount = requires (Element element) { element.basis_functions_n(); };

  template<typename Element>
  concept HasLumped = HasBasisFunctionAmount<Element>
    && requires (Element element) { element.lumped(0), element.lumped_; };
}

#endif //STG_CONCEPTS_HPP
