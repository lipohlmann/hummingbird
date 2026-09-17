// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, Liam Pohlmann

#include "utils/output.h"

#include <gtest/gtest.h>

#include <string>

#include "utils/enums.h"

namespace hummingbird {

TEST(OutputTest, PrintHeaderPrintsDescriptionAndCopyright) {
  testing::internal::CaptureStdout();
  print_header();
  const std::string output = testing::internal::GetCapturedStdout();

  EXPECT_NE(output.find("Spectral Elements for Radiation Transport"),
           std::string::npos);
  EXPECT_NE(output.find("Copyright"), std::string::npos);
}

TEST(OutputTest, PrintColumnsPrintsBothColumnHeaders) {
  testing::internal::CaptureStdout();
  print_columns();
  const std::string output = testing::internal::GetCapturedStdout();

  EXPECT_NE(output.find("k-eff Iter."), std::string::npos);
  EXPECT_NE(output.find("Scatter Iter."), std::string::npos);
}

TEST(OutputTest, PrintInputFilesPrintsBothFileNames) {
  testing::internal::CaptureStdout();
  print_input_files("input.json", "mesh.msh");
  const std::string output = testing::internal::GetCapturedStdout();

  EXPECT_NE(output.find("input.json"), std::string::npos);
  EXPECT_NE(output.find("mesh.msh"), std::string::npos);
}

TEST(OutputTest, PrintKStatusFirstIterationShowsPlaceholderInsteadOfError) {
  testing::internal::CaptureStdout();
  print_k_status(1.0, 0.5, 1);
  const std::string output = testing::internal::GetCapturedStdout();

  EXPECT_NE(output.find("-"), std::string::npos);
}

TEST(OutputTest, PrintKStatusLaterIterationShowsKEffAndError) {
  testing::internal::CaptureStdout();
  print_k_status(1.23456, 0.01, 3);
  const std::string output = testing::internal::GetCapturedStdout();

  EXPECT_NE(output.find("3"), std::string::npos);
  EXPECT_NE(output.find("1.23456e+00"), std::string::npos);
}

TEST(OutputTest, PrintScatterStatusFirstIterationShowsPlaceholderInsteadOfError) {
  testing::internal::CaptureStdout();
  print_scatter_status(1.0, 0.5, 1);
  const std::string output = testing::internal::GetCapturedStdout();

  EXPECT_NE(output.find("-"), std::string::npos);
}

TEST(OutputTest, PrintScatterStatusLaterIterationShowsScatterAndError) {
  testing::internal::CaptureStdout();
  print_scatter_status(2.5, 0.02, 4);
  const std::string output = testing::internal::GetCapturedStdout();

  EXPECT_NE(output.find("4"), std::string::npos);
  EXPECT_NE(output.find("2.50000e+00"), std::string::npos);
}

TEST(OutputTest, PrintScatterCompleteAnnouncesCompletion) {
  testing::internal::CaptureStdout();
  print_scatter_complete(1.0, RunMode::FIXED_SOURCE);
  const std::string output = testing::internal::GetCapturedStdout();

  EXPECT_NE(output.find("Source iterations complete."), std::string::npos);
}

}  // namespace hummingbird
