# external/CMakeLists

find_package(pugixml REQUIRED)
target_link_libraries(hummingbird_lib PRIVATE pugixml::pugixml)

find_package(Armadillo CONFIG REQUIRED)

if(Armadillo_FOUND AND NOT TARGET Armadillo::Armadillo)
  add_library(Armadillo::Armadillo INTERFACE IMPORTED)
  set_target_properties(
    Armadillo::Armadillo
    PROPERTIES
      INTERFACE_LINK_LIBRARIES "${ARMADILLO_LIBRARIES}"
      INTERFACE_INCLUDE_DIRECTORIES "${ARMADILLO_INCLUDE_DIRS}"
  )
endif()
target_link_libraries(hummingbird_lib PRIVATE Armadillo::Armadillo)

