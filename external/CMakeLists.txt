# external/CMakeLists

find_package(pugixml REQUIRED)
target_link_libraries(hummingbird_lib PRIVATE pugixml::pugixml)

find_package(Armadillo REQUIRED)
target_link_libraries(hummingbird_lib PRIVATE Armadillo)

