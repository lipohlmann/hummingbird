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
   */
  Element(const int material_id);

  /**
   * @brief Create nodes by mapping the Gauss-Lobatto-Legendre quadrature set
   * from the reference to the real domain
   *
   * @param gll_quadrature GaussLobattoLegendre quadrature set
   */
  virtual std::vector<Node> CreateInteriorNodes(
      const quadrature::GaussLobattoLegendre& gll_quadrature) = 0;

  /**
   * @brief Get node IDs
   *
   * @return std::vector<size_t>
   */
  std::vector<size_t> node_ids() const { return node_ids_; }

  /**
   * @brief Set the new Node ID
   *
   * @param prev_id Previous ID (currently stored in object)
   * @param new_id New ID
   */
  void SetNewNodeID(const size_t prev_id, const size_t new_id);

 protected:
  /// @brief Material ID
  const int material_id_;

  /// @brief IDs of nodes defining the element
  std::vector<size_t> node_ids_;
};
}  // namespace hummingbird::mesh

#endif  // HUMMINGBIRD_MESH_ELEMENT_H_
