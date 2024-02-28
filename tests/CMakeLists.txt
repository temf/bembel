
# copy geometry files for testing
configure_file(${CMAKE_SOURCE_DIR}/geo/sphere.dat ${CMAKE_CURRENT_BINARY_DIR}/ COPYONLY)

set(UNITTESTS
		test_GeometryImportAndEval
		test_GeometryIGS
		test_GeometryIO
		test_Bezierextraction
		test_SurfacePointUpdate
		test_Spline
		test_Bernstein
		test_Projector
		test_Glue
		test_TangentialTrace
		test_FMMTransferMatrices
		test_FMMForwardTransformation
		test_FMMBackwardTransformation
		test_Laplace
		test_Helmholtz
		test_DuffyTrick
		test_LazyEigenSum
		test_HomogenisedCoefficients
		)

foreach(file IN LISTS UNITTESTS)
  add_executable(${file} ${file}.cpp)
  target_link_libraries(${file} Eigen3::Eigen)
  add_test(${file}Unittest ${file})
  set_tests_properties(${file}Unittest PROPERTIES LABELS "unittests")
endforeach()

enable_testing()

add_subdirectory(examples_as_unittests)
set(EXAMPLES_CACHE ${EXAMPLES} CACHE INTERNAL "My set")

add_custom_target(unittests)
add_dependencies(unittests ${UNITTESTS} ${EXAMPLES_CACHE})
