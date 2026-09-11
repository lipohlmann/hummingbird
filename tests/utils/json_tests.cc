// Unit tests for starling::utils::JSONFromFile.

#include <gtest/gtest.h>

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <string>

#include "utils/json.h"

namespace starling::utils {

// Helper: writes `contents` to a temp file and returns its path. File is
// removed by the returned guard going out of scope.
class TempFile {
 public:
  TempFile(const std::string& contents,
           const std::string& extension = ".json") {
    path_ = std::filesystem::temp_directory_path() /
            (std::string("starling_json_test_") +
             std::to_string(reinterpret_cast<uintptr_t>(this)) + extension);
    std::ofstream out(path_);
    out << contents;
  }

  ~TempFile() { std::filesystem::remove(path_); }

  const std::string path() const { return path_.string(); }

 private:
  std::filesystem::path path_;
};

TEST(JSONFromFile, ParsesValidObject) {
  TempFile file(R"({"a": 1, "b": "two", "c": [1, 2, 3]})");
  json result = JSONFromFile(file.path());
  EXPECT_EQ(result["a"], 1);
  EXPECT_EQ(result["b"], "two");
  EXPECT_EQ(result["c"], std::vector<int>({1, 2, 3}));
}

TEST(JSONFromFile, ParsesValidArrayAtTopLevel) {
  TempFile file(R"([1, 2, 3])");
  json result = JSONFromFile(file.path());
  EXPECT_TRUE(result.is_array());
  EXPECT_EQ(result.size(), 3u);
}

TEST(JSONFromFile, ParsesNestedStructures) {
  TempFile file(R"({"outer": {"inner": [true, false, null]}})");
  json result = JSONFromFile(file.path());
  ASSERT_TRUE(result.contains("outer"));
  EXPECT_TRUE(result["outer"]["inner"][0].get<bool>());
  EXPECT_FALSE(result["outer"]["inner"][1].get<bool>());
  EXPECT_TRUE(result["outer"]["inner"][2].is_null());
}

TEST(JSONFromFile, ThrowsOnNonexistentFile) {
  // No TempFile created -- path should not exist.
  EXPECT_THROW(
      { JSONFromFile("/this/path/almost_certainly_does_not_exist.json"); },
      std::runtime_error);
}

TEST(JSONFromFile, ThrowsOnNonexistentFileWithAbsolutePathInMessage) {
  // The implementation resolves and embeds an absolute path in the error
  // message even when given a relative filename. Confirms that codepath
  // doesn't itself throw/crash (e.g. std::filesystem::absolute on a bad
  // relative path) and that the message contains a path, not just a bare
  // string like the raw input.
  const std::string relative_missing = "definitely_missing_file.json";
  try {
    JSONFromFile(relative_missing);
    FAIL() << "Expected std::runtime_error to be thrown";
  } catch (const std::runtime_error& e) {
    std::string message = e.what();
    EXPECT_NE(message.find("Could not open mesh file:"), std::string::npos);
  }
}

TEST(JSONFromFile, ThrowsNlohmannParseErrorOnMalformedJson) {
  // json::parse throws nlohmann::json::parse_error (not caught/rethrown as
  // std::runtime_error here), so malformed-but-openable files surface a
  // different exception type than the "file not found" case. Documenting
  // that asymmetry via this test rather than assuming callers expect a
  // single exception type from this function.
  TempFile file("{ this is not valid json ");
  EXPECT_THROW({ JSONFromFile(file.path()); }, nlohmann::json::parse_error);
}

TEST(JSONFromFile, ThrowsOnEmptyFile) {
  // An empty file is openable but has no JSON content to parse.
  TempFile file("");
  EXPECT_THROW({ JSONFromFile(file.path()); }, nlohmann::json::parse_error);
}

TEST(JSONFromFile, HandlesEmptyJsonObject) {
  TempFile file("{}");
  json result = JSONFromFile(file.path());
  EXPECT_TRUE(result.is_object());
  EXPECT_TRUE(result.empty());
}

TEST(JSONFromFile, AcceptsPathWithDirectoryComponent) {
  // Confirms filename isn't assumed to be a bare name (e.g. no directory
  // separators) -- constructs a temp file under a nested directory.
  std::filesystem::path dir =
      std::filesystem::temp_directory_path() / "starling_json_test_subdir";
  std::filesystem::create_directories(dir);
  std::filesystem::path nested_file = dir / "nested.json";
  {
    std::ofstream out(nested_file);
    out << R"({"nested": true})";
  }

  json result = JSONFromFile(nested_file.string());
  EXPECT_TRUE(result["nested"].get<bool>());

  std::filesystem::remove_all(dir);
}

}  // namespace starling::utils