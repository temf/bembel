set(UNITTESTS
		test_GeometryImportAndEval
		test_SurfacePointUpdate
		test_Projector
		test_Glue
		test_FMMTransferMatrices
		test_FMMForwardTransformation
		test_FMMBackwardTransformation
		test_DuffyTrick
		)

foreach(file IN LISTS UNITTESTS)
  add_executable(${file} ${file}.cpp)
  target_link_libraries(${file} Eigen3::Eigen)
  add_test(${file}Unittest ${file})
  set_tests_properties(${file}Unittest PROPERTIES LABELS "unittests")
endforeach()

enable_testing()

add_custom_target(unittests)
add_dependencies(unittests ${UNITTESTS})