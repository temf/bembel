// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.

#include <Bembel/Geometry>
#include "Test.hpp"

int main() {
  using namespace Bembel;

  constexpr std::array<double, 11> eq_points = {0,   0.1, 0.2, 0.3, 0.4, 0.5,
                                                0.6, 0.7, 0.8, 0.9, 1};
  Test::TestGeometryWriter::writeSpherePanel();

  Bembel::Geometry geometry("test_SpherePanel.dat");
  assert(geometry.get_geometry().size() == 1);

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
  Test::TestGeometryWriter::writeScaledSpherePanel();

  Bembel::Geometry geometry_fail("test_ScaledSpherePanel.dat");
  assert(geometry.get_geometry().size() == 1);

  for (auto x : Test::Constants::eq_points) {
    for (auto y : Test::Constants::eq_points) {
      auto point = geometry.get_geometry()[0].eval(Eigen::Vector2d(x, y));
      auto normal = geometry.get_geometry()[0]
                        .evalNormal(Eigen::Vector2d(x, y))
                        .normalized();
      BEMBEL_TEST_IF((std::abs((point.norm() - 1.0)) <
                      Test::Constants::test_tolerance_geometry));
      BEMBEL_TEST_IF((point - normal).norm() <
                     Test::Constants::test_tolerance_geometry);
    }
  }

  return 0;
}
