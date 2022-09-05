
# copy geometry files for testing
configure_file(${CMAKE_SOURCE_DIR}/geo/sphere.dat ${CMAKE_CURRENT_BINARY_DIR}/ COPYONLY)
configure_file(${CMAKE_SOURCE_DIR}/geo/torus.dat ${CMAKE_CURRENT_BINARY_DIR}/ COPYONLY)

###############################################################################
# examples being run from CI
###############################################################################

set(CIFILES
    AnsatzSpace
    BlockClusterTree
    Geometry
    LaplaceSingleLayerFull
    LaplaceSingleLayerH2
    HelmholtzSingleLayerFull
    HelmholtzSingleLayerH2
    MaxwellSingleLayerFull
    MaxwellSingleLayerH2
    )

foreach(file IN LISTS CIFILES)
  add_executable(${file} ${file})
  target_link_libraries(${file} Eigen3::Eigen)
  add_test(${file}Example ${file})
  set_tests_properties(${file}Example PROPERTIES LABELS "E2E")
endforeach()

enable_testing()

add_custom_target(examples)
add_dependencies(examples ${CIFILES})

###############################################################################
# examples not being run from CI
###############################################################################

set(NOCIFILES
    Quadrature
    VTKSurfaceExport
    VTKDomainExport
    FullLaplaceWorkflow
    )

  foreach(file IN LISTS NOCIFILES)
  add_executable(${file}.out ${file})
  target_link_libraries(${file}.out Eigen3::Eigen)
endforeach()
