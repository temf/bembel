// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.

#ifndef BEMBEL_TEST_TEST_H_
#define BEMBEL_TEST_TEST_H_

#include <array>
#include <fstream>
#include <iostream>
#include <string>

#include "TestGeometries.hpp"

#define BEMBEL_TEST_IF(in_bool)                                                \
  {                                                                            \
    if (!(in_bool)) {                                                        \
      assert(in_bool);                                                         \
      std::cout << "A test failed. Please recompile without -DNDEBUG to see "; \
      std::cout << "the corresponding assert.";                                \
      std::cout << std::endl;                                                  \
      return 1;                                                                \
    }                                                                          \
  }

namespace Test {
namespace Constants {
constexpr std::array<double, 11> eq_points = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5,
                                              0.6, 0.7, 0.8, 0.9, 1.0};
constexpr double test_tolerance_geometry = 1e-9;

constexpr double coefficient_accuracy = 1e-6;

// if this is defined, an anlytic testfunction is used for testing the
// Duffytrick, otherwise, a nearly singular function is used.
// #define BEMBEL_TEST_DUFFYTRICK_USE_ANALYTIC_FUNCTION_
constexpr unsigned int DuffyTrick_quadrature_degree = 16;
constexpr unsigned int DuffyTrick_refinement_level = 2;
constexpr double integrate0_tolerance = 1e-12;
constexpr double integrate1_tolerance = 1e-12;
constexpr double integrate2_tolerance = 1e-6;
constexpr double integrate3_tolerance = 1e-4;
constexpr double integrate4_tolerance = 1e-4;
}  // namespace Constants
}  // namespace Test

#endif
