// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, Liam Pohlmann

#include <gtest/gtest.h>

#include <stdexcept>
#include <vector>

#include "banks/material_bank.h"
#include "mesh/element.h"
#include "mesh/node.h"
#include "quadrature/angular/ordinate.h"
#include "quadrature/gauss_lobatto_legendre.h"

namespace hummingbird {

namespace {

// Element is abstract, so its base-class behavior (node/material/source ID
// bookkeeping) is exercised through this minimal concrete subclass, which
// trivially implements the pure-virtual, formulation-specific methods.
class TestElement : public Element {
 public:
  using Element::Element;

  std::vector<Node> CreateInteriorNodes(
      const std::vector<Node>& existing_nodes,
      const GaussLobattoLegendre& gll_quadrature) override {
    return {};
  }

  arma::Col<double> LocalForcingVector(const GaussLobattoLegendre& gll_quad,
                                       const Mesh& mesh,
                                       const size_t ordinate_index) override {
    return arma::Col<double>();
  }

  arma::Mat<double> LocalStiffnessMatrix(const GaussLobattoLegendre& gll_quad,
                                         const MaterialBank& material_bank,
                                         const Ordinate& ordinate) override {
    return arma::Mat<double>();
  }

  arma::SpMat<double> LocalMassMatrix(
      const GaussLobattoLegendre& gll_quad,
      const MaterialBank& material_bank) override {
    return arma::SpMat<double>();
  }

  unsigned int dimension() const override { return 1; }
};

}  // namespace

TEST(ElementTest, ConstructorSetsMaterialAndSourceID) {
  TestElement element(3, 7);
  EXPECT_EQ(element.material_id(), 3);
  EXPECT_EQ(element.source_id(), 7);
}

TEST(ElementTest, NodeIdsStartsEmpty) {
  TestElement element(0, 0);
  EXPECT_TRUE(element.node_ids().empty());
}

TEST(ElementTest, SetMaterialIDUpdatesMaterialID) {
  TestElement element(1, 2);
  element.SetMaterialID(9);
  EXPECT_EQ(element.material_id(), 9);
}

TEST(ElementTest, SetSourceIDUpdatesSourceID) {
  TestElement element(1, 2);
  element.SetSourceID(9);
  EXPECT_EQ(element.source_id(), 9);
}

TEST(ElementTest, SetNodeIDsReplacesAllNodeIDs) {
  TestElement element(0, 0);
  element.SetNodeIDs({0, 1, 2});
  EXPECT_EQ(element.node_ids(), (std::vector<size_t>{0, 1, 2}));

  element.SetNodeIDs({5, 6});
  EXPECT_EQ(element.node_ids(), (std::vector<size_t>{5, 6}));
}

TEST(ElementTest, SetNewNodeIDReplacesMatchingID) {
  TestElement element(0, 0);
  element.SetNodeIDs({0, 1, 2});

  element.SetNewNodeID(1, 42);

  EXPECT_EQ(element.node_ids(), (std::vector<size_t>{0, 42, 2}));
}

TEST(ElementTest, SetNewNodeIDThrowsWhenPreviousIDNotFound) {
  TestElement element(0, 0);
  element.SetNodeIDs({0, 1, 2});

  EXPECT_THROW(element.SetNewNodeID(99, 42), std::runtime_error);
}

}  // namespace hummingbird
