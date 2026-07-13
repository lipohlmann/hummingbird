# tests/CMakeLists.txt

include(FetchContent)

find_package(GTest QUIET)
if(NOT GTest_FOUND)
        FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/5a9c3f9e8d9b90bbbe8feb32902146cb8f7c1757.zip
        )
        FetchContent_MakeAvailable(googletest)
endif()
include(GoogleTest)

file(GLOB_RECURSE TEST_SOURCES CONFIGURE_DEPENDS "*.cc")
add_executable(test_hummingbird ${TEST_SOURCES})
target_include_directories(test_hummingbird PRIVATE ${INCLUDE_DIR})
target_compile_definitions(test_hummingbird PRIVATE TEST_DIR="${TEST_DIR}")
target_link_libraries(test_hummingbird PRIVATE GTest::gtest_main hummingbird_lib)

gtest_discover_tests(test_hummingbird)

