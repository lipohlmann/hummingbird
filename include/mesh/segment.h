// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, Liam Pohlmann

#ifndef HUMMINGBIRD_MESH_SEGMENT_H_
#define HUMMINGBIRD_MESH_SEGMENT_H_

#include <array>

#include "mesh/element.h"
#include "mesh/node.h"

namespace hummingbird {
class Segment : public Element {
 public:
  /**
   * @brief Construct a new Segment object
   *
   * @param boundary_node_ids Node IDs defining the Segment bounds
   * @param material_id Material ID
   * @param source_id Source ID
   */
  Segment(const std::array<size_t, 2> boundary_node_ids, const int material_id,
          const int source_id);

  /**
   * @brief Create interior Node objects using the Gauss-Lobatto-Legendre
   * quadrature set. The mapping from the reference to the real domain was based
   * on the equations provided by the book Computational Seismology by Heiner
   * Igel
   * (https://www.google.com/url?sa=t&source=web&rct=j&opi=89978449&url=https://www.geokniga.org/bookfiles/geokniga-computationalseismologyapracticalintroductionbyigelheinerz-liborg.pdf&ved=2ahUKEwiMuanq7_OVAxXVN4YAHVgAD1QQFnoECBUQAQ&usg=AOvVaw35g1Tl3ZSFUMusz8LKQxHk).
   * For a line segment defined by two endpoints, the mapping is:
   *
   * \f[
   * \mathbf{r}(\xi) = \mathbf{r}_1+\left( \frac{\xi+1}{2} \right)\left(
   * \mathbf{r}_2-\mathbf{r}_1 \right) \f]
   *
   * where \f$x\in[-1,1]\f$.
   *
   * @param existing_nodes Nodes already existing in Mesh
   * @param gll_quadrature Gauss-Lobatto-Legendre quadrature set
   * @return Additional nodes to be defined along the GLL quadrature set. Will
   * be of length N-2 (the two endpoint account for the remaining two nodes)
   */
  std::vector<Node> CreateInteriorNodes(
      const std::vector<Node>& existing_nodes,
      const GaussLobattoLegendre& gll_quadrature) override;

 private:
  std::array<size_t, 2> boundary_node_ids_;
};
}  // namespace hummingbird

#endif  // HUMMINGBIRD_MESH_SEGMENT_H_
