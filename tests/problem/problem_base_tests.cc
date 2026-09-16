#include "problem/problem_base.h"

#include <gtest/gtest.h>

#include <stdexcept>
#include <vector>

namespace hummingbird {

namespace {

// ProblemBase is abstract (AssembleGlobalMatrixData/AssembleGlobalForcingData
// are pure virtual), so this stub exists purely to get a concrete instance
// whose CheckGlobalMatrixData/CheckGlobalForcingData (protected, unmodified
// here) can be exercised directly. CGProblem isn't used for this because its
// AssembleGlobalForcingData override isn't implemented yet.
class TestableProblem : public ProblemBase {
 public:
  // n_dofs/n_ordinates are irrelevant here since these tests only exercise
  // CheckGlobalMatrixData/CheckGlobalForcingData directly.
  TestableProblem() : ProblemBase(0, 0) {}

  using ProblemBase::CheckGlobalForcingData;
  using ProblemBase::CheckGlobalMatrixData;

  std::vector<GlobalMatrixData> AssembleGlobalMatrixData(
      const Mesh&, const GaussLobattoLegendre&, const MaterialBank&,
      const Ordinate&) override {
    return {};
  }

  std::vector<GlobalForcingData> AssembleGlobalForcingData(
      const GaussLobattoLegendre&, const Mesh&, const size_t) override {
    return {};
  }

  void Apply1DBCs(const Mesh&, const QuadratureBase<Ordinate>,
                  const BCBank&) override {}
};

}  // namespace

// ---------------------------------------------------------------------------
// CheckGlobalForcingData validation tests
// ---------------------------------------------------------------------------

TEST(CheckGlobalForcingDataTest, AcceptsContiguousUniqueRowIDs) {
  TestableProblem problem;
  std::vector<GlobalForcingData> gfd = {{0, 1.0}, {1, 2.0}, {2, 3.0}};
  EXPECT_NO_THROW(problem.CheckGlobalForcingData(gfd));
}

TEST(CheckGlobalForcingDataTest, AcceptsEmptyData) {
  TestableProblem problem;
  std::vector<GlobalForcingData> gfd;
  EXPECT_NO_THROW(problem.CheckGlobalForcingData(gfd));
}

TEST(CheckGlobalForcingDataTest, ThrowsOnDuplicateRowID) {
  TestableProblem problem;
  std::vector<GlobalForcingData> gfd = {{0, 1.0}, {0, 2.0}, {1, 3.0}};
  EXPECT_THROW(problem.CheckGlobalForcingData(gfd), std::runtime_error);
}

TEST(CheckGlobalForcingDataTest, ThrowsOnGapInRowIDs) {
  TestableProblem problem;
  std::vector<GlobalForcingData> gfd = {{0, 1.0}, {1, 2.0}, {3, 3.0}};
  EXPECT_THROW(problem.CheckGlobalForcingData(gfd), std::runtime_error);
}

// ---------------------------------------------------------------------------
// CheckGlobalMatrixData validation tests
// ---------------------------------------------------------------------------

TEST(CheckGlobalMatrixDataTest, AcceptsContiguousDenseIDs) {
  TestableProblem problem;
  std::vector<GlobalMatrixData> gmd = {
      {0, 0, 1.0}, {0, 1, 2.0}, {1, 0, 2.0}, {1, 1, 3.0}};
  EXPECT_NO_THROW(problem.CheckGlobalMatrixData(gmd));
}

TEST(CheckGlobalMatrixDataTest, AcceptsEmptyData) {
  TestableProblem problem;
  std::vector<GlobalMatrixData> gmd;
  EXPECT_NO_THROW(problem.CheckGlobalMatrixData(gmd));
}

TEST(CheckGlobalMatrixDataTest, AcceptsDuplicateRowColPairsFromSharedNodes) {
  // Two elements both contributing to node 1 (a shared node), producing
  // duplicate (row_id, col_id) pairs -- this is expected/correct FE assembly
  // behavior (Armadillo's batch SpMat constructor sums the duplicates) and
  // must not throw.
  TestableProblem problem;
  std::vector<GlobalMatrixData> gmd = {
      {0, 0, 1.0}, {0, 1, 1.0}, {1, 0, 1.0}, {1, 1, 1.0},  // element A: 0, 1
      {1, 1, 2.0}, {1, 2, 2.0}, {2, 1, 2.0}, {2, 2, 2.0},  // element B: 1, 2
  };
  EXPECT_NO_THROW(problem.CheckGlobalMatrixData(gmd));
}

TEST(CheckGlobalMatrixDataTest, ThrowsWhenReferencedIDsHaveGap) {
  // Only node IDs 0 and 3 are referenced; 1 and 2 are missing.
  TestableProblem problem;
  std::vector<GlobalMatrixData> gmd = {
      {0, 0, 1.0}, {0, 3, 2.0}, {3, 0, 2.0}, {3, 3, 3.0}};
  EXPECT_THROW(problem.CheckGlobalMatrixData(gmd), std::runtime_error);
}

}  // namespace hummingbird
