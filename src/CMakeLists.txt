file(GLOB_RECURSE SOURCES "*.cc")

add_library(hummingbird_lib)
target_compile_features(hummingbird_lib PRIVATE cxx_std_20)
target_include_directories(hummingbird_lib PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_sources(hummingbird_lib PRIVATE ${SOURCES})

add_executable(hummingbird main.cc)
target_link_libraries(hummingbird PRIVATE hummingbird_lib)

