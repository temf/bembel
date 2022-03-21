
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
    LaplaceDoubleLayerH2
    HelmholtzSingleLayerFull
    HelmholtzSingleLayerH2
    HelmholtzDoubleLayerH2
    HelmholtzBrakhageWernerH2
    MaxwellSingleLayerFull
    MaxwellSingleLayerH2
    )

foreach(file IN LISTS CIFILES)
  add_executable(${file}.out ${file})
  target_link_libraries(${file}.out Eigen3::Eigen)
  add_test(${file}Example ${file}.out)
endforeach()

enable_testing()

###############################################################################
# examples not being run from CI
###############################################################################

set(NOCIFILES
    Quadrature
    MassMatrix
    VTKSurfaceExport
    VTKDomainExport
    FullLaplaceWorkflow
    )

  foreach(file IN LISTS NOCIFILES)
  add_executable(${file}.out ${file})
  target_link_libraries(${file}.out Eigen3::Eigen)
endforeach()
