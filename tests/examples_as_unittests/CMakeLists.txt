configure_file(${CMAKE_SOURCE_DIR}/geo/sphere.dat ${CMAKE_CURRENT_BINARY_DIR}/ COPYONLY)

set(EXAMPLES
    test_MassMatrix
		test_MaxwellSingleLayerFull
		)

foreach(file IN LISTS EXAMPLES)
  add_executable(${file} ${file}.cpp)
  target_link_libraries(${file} Eigen3::Eigen)
  add_test(${file}Unittest ${file})
  set_tests_properties(${file}Unittest PROPERTIES LABELS "unittests")
endforeach()

enable_testing()
