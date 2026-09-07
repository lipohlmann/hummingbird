find_package(Armadillo CONFIG REQUIRED)

if(Armadillo_FOUND AND NOT TARGET Armadillo::Armadillo)
  add_library(Armadillo::Armadillo INTERFACE IMPORTED)
  set_target_properties(
    Armadillo::Armadillo
    PROPERTIES INTERFACE_LINK_LIBRARIES "${ARMADILLO_LIBRARIES}"
               INTERFACE_INCLUDE_DIRECTORIES "${ARMADILLO_INCLUDE_DIRS}")
endif()

find_package(nlohmann_json REQUIRED)
find_package(fmt)

add_library(
  hummingbird_lib
  math/legendre_polynomials.cc
  mesh/element.cc
  mesh/mesh.cc
  mesh/segment.cc
  physics/parsed_volumetric_source.cc
  quadrature/gauss_chebyshev.cc
  quadrature/gauss_legendre.cc
  quadrature/gauss_lobatto_legendre.cc
  quadrature/angular/gauss_legendre_chebyshev.cc
  quadrature/angular/gauss_legendre_trapezoid.cc
  quadrature/angular/ordinate.cc
  utils/misc.cc)

target_compile_features(hummingbird_lib PRIVATE cxx_std_20)

target_link_libraries(hummingbird_lib PRIVATE nlohmann_json::nlohmann_json)
target_link_libraries(hummingbird_lib PRIVATE Armadillo::Armadillo)
target_link_libraries(hummingbird_lib PRIVATE fmt::fmt)

target_include_directories(hummingbird_lib
                           PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_sources(hummingbird_lib PRIVATE ${SOURCES})

add_executable(hummingbird main.cc)
target_link_libraries(hummingbird PRIVATE hummingbird_lib)
