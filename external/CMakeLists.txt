# external/CMakeLists

find_package(pugixml REQUIRED)
target_link_libraries(hummingbird_lib PRIVATE pugixml::pugixml)

find_package(armadillo REQUIRED)
target_link_libraries(hummingbird_lib PRIVATE armadillo)

