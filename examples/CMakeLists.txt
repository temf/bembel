
# copy geometry files for testing
configure_file(${CMAKE_SOURCE_DIR}/geo/sphere.dat ${CMAKE_CURRENT_BINARY_DIR}/ COPYONLY)
configure_file(${CMAKE_SOURCE_DIR}/geo/torus.dat ${CMAKE_CURRENT_BINARY_DIR}/ COPYONLY)
configure_file(${CMAKE_SOURCE_DIR}/geo/cube_small.dat ${CMAKE_CURRENT_BINARY_DIR}/ COPYONLY)

###############################################################################
# examples being run from CI
###############################################################################

set(CIFILES
    AnsatzSpace
    BlockClusterTree
    Geometry
    MassMatrix
    LaplaceBeltrami
    LaplaceSingleLayerFull
    LaplaceSingleLayerH2
    LaplaceHypersingularH2
    LaplaceDoubleLayerH2
    LaplaceAdjointDoubleLayerH2
    LazyEigenSum
    HelmholtzSingleLayerFull
    HelmholtzSingleLayerH2
    HelmholtzDoubleLayerH2
    HelmholtzAdjointDoubleLayerH2
    HelmholtzHypersingularH2
    HomogenisedLaplaceSingleLayerFull
    HomogenisedLaplaceSingleLayerH2
    MaxwellSingleLayerFull
    MaxwellSingleLayerH2
    )

foreach(file IN LISTS CIFILES)
  add_executable(${file} ${file}.cpp)
  target_link_libraries(${file} Eigen3::Eigen)
  add_test(${file}Example ${file})
  set_tests_properties(${file}Example PROPERTIES LABELS "E2Etests")
endforeach()

enable_testing()

add_custom_target(E2Etests)
add_dependencies(E2Etests ${CIFILES})

###############################################################################
# examples not being run from CI
###############################################################################

set(NOCIFILES
    Quadrature
    VTKSurfaceExport
    VTKDomainExport
    FullLaplaceWorkflow
    FullMaxwellWorkflow
    )

  foreach(file IN LISTS NOCIFILES)
  add_executable(${file} ${file}.cpp)
  target_link_libraries(${file} Eigen3::Eigen)
endforeach()
