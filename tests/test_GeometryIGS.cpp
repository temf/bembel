// This file is part of Bembel, the higher order C++ boundary element library.
//
// Copyright (C) 2024 see <http://www.bembel.eu>
//
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.

#include <Bembel/Geometry>

#include "tests/Test.hpp"
#include "tests/TestHelpers.hpp"

int main() {
  using namespace Bembel;

  Test::TestGeometryWriter::writeIGSScreen();

  Bembel::Geometry geometry("test_Screen.igs");
  assert(geometry.get_geometry().size() == 2);

  const PatchVector& patches = geometry.get_geometry();

  BEMBEL_TEST_IF(patches[0].polynomial_degree_x_ == 5);
  BEMBEL_TEST_IF(patches[0].polynomial_degree_y_ == 6);

  writeIGSFile(patches, "test_Screen_Export.igs");

  BEMBEL_TEST_IF(
      Test::compareFiles("test_Screen.igs", "test_Screen_Export.igs"));

  for (auto x : Test::Constants::eq_points) {
    for (auto y : Test::Constants::eq_points) {
      auto point = geometry.get_geometry()[0].eval(Eigen::Vector2d(x, y));
      auto normal = geometry.get_geometry()[0]
                        .evalNormal(Eigen::Vector2d(x, y))
                        .normalized();
      BEMBEL_TEST_IF(std::abs((point.norm() - 1.0)) <
                     Test::Constants::test_tolerance_geometry);
      BEMBEL_TEST_IF((point - normal).norm() <
                     Test::Constants::test_tolerance_geometry);
    }
  }

  return 0;
}
