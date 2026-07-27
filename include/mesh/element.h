#ifndef HUMMINGBIRD_MESH_ELEMENT_H_
#define HUMMINGBIRD_MESH_ELEMENT_H_

#include <memory>
#include <vector>

#include "material.h"
#include "node.h"
#include "quadrature/gauss_lobatto_legendre.h"

namespace hummingbird::mesh {
/**
 * @brief Defines a subset of the domain (an "element")
 *
 */
class Element {
 public:
  /**
   * @brief Construct a new Element object
   *
   * @param material_id Material ID
   * @param gll_quadrature GaussLobattoLegendre quadrature set
   */
  Element(const int material_id,
          const quadrature::GaussLobattoLegendre& gll_quadrature);

  /**
   * @brief Create nodes by mapping the Gauss-Lobatto-Legendre quadrature set
   * from the reference to the real domain
   *
   * @param gll_quadrature GaussLobattoLegendre quadrature set
   */
  virtual std::vector<Node> CreateNodes(
      const quadrature::GaussLobattoLegendre& gll_quadrature) = 0;

 protected:
  /// @brief Material ID
  const int material_id_;

  /// @brief Gauss-Lobatto-Legendre quadrature set
  const quadrature::GaussLobattoLegendre& gll_quadrature_;

  /// @brief Nodes defining the element
  std::vector<Node> nodes_;
};
}  // namespace hummingbird::mesh

#endif  // HUMMINGBIRD_MESH_ELEMENT_H_
