// This file is part of Bembel, the higher order C++ boundary element library.
//
// Copyright (C) 2022 see <http://www.bembel.eu>
//
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.

#include <Bembel/Geometry>

#include "tests/Test.hpp"

int main() {
  using namespace Bembel;

  // initialize geometry
  Test::TestGeometryWriter::writeScreen();
  Bembel::Geometry geometry("test_Screen.dat");
  assert(geometry.get_geometry().size() == 1);

  // initialize tolerance
  double tol = Test::Constants::test_tolerance_geometry;

  for (auto x : Test::Constants::eq_points) {
    for (auto y : Test::Constants::eq_points) {
      auto pt = Eigen::Vector2d(x, y);
      SurfacePoint srf_pt;
      geometry.get_geometry()[0].updateSurfacePoint(&srf_pt, pt, 3.1415, pt);
      assert((srf_pt.get_xi() - pt).norm() < tol);
      assert(std::abs(srf_pt.get_w() - 3.1415) < tol);
      assert((srf_pt.get_f() - geometry.get_geometry()[0].eval(pt)).norm() <
             tol);
      assert(
          (srf_pt.get_jacobian() - geometry.get_geometry()[0].evalJacobian(pt))
              .norm() < tol);
    }
  }
  return 0;
}
