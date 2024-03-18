configure_file(${CMAKE_SOURCE_DIR}/geo/sphere.dat ${CMAKE_CURRENT_BINARY_DIR}/ COPYONLY)
configure_file(${CMAKE_SOURCE_DIR}/geo/cube_small.dat ${CMAKE_CURRENT_BINARY_DIR}/ COPYONLY)

set(EXAMPLES
    test_MassMatrix
    test_LaplaceBeltrami
    test_LaplaceSingleLayerH2
    test_HelmholtzSingleLayerH2
		test_MaxwellSingleLayerFull
    test_MaxwellSingleLayerH2
    test_HomogenisedLaplaceSingleLayerFull
    test_HomogenisedLaplaceSingleLayerH2
    CACHE INTERNAL "Examples as Unittests"
		)

foreach(file IN LISTS EXAMPLES)
  add_executable(${file} ${file}.cpp)
  target_link_libraries(${file} Eigen3::Eigen)
  add_test(${file}Unittest ${file})
  set_tests_properties(${file}Unittest PROPERTIES LABELS "unittests")
endforeach()

enable_testing()
