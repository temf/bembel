
###############################################################################
# examples not being run from CI
###############################################################################

set(OTHER_NOCIFILES
    AugmentedEFIE
    AugmentedEFIEStabilized
    AugmentedEFIEBalun
    )

foreach(file IN LISTS OTHER_NOCIFILES)
  add_executable(${file} ${file}.cpp)
  target_link_libraries(${file} Eigen3::Eigen)
endforeach()
