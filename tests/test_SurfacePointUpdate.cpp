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

  Test::TestGeometryWriter::writeScreen();

  Bembel::Geometry geometry("test_Screen.dat");
  assert(geometry.get_geometry().size() == 1);

  for (auto x : Test::Constants::eq_points) {
    for (auto y : Test::Constants::eq_points) {
      auto pt = Eigen::Vector2d(x, y);
      Eigen::Matrix<double, 12, 1> srf_pt, srf_pt_ref;

      srf_pt_ref.head(2) = pt;
      srf_pt_ref(2) = 3.1415;
      srf_pt_ref.segment(3, 3) = geometry.get_geometry()[0].eval(pt);
      auto dummy = geometry.get_geometry()[0].evalJacobian(pt);
      srf_pt_ref.segment(6, 3) = dummy.col(0);
      srf_pt_ref.segment(9, 3) = dummy.col(1);

      auto point = geometry.get_geometry()[0].eval(pt);
      auto jacobian = geometry.get_geometry()[0].evalJacobian(pt);

      geometry.get_geometry()[0].updateSurfacePoint(&srf_pt, pt, 3.1415, pt);

      if ((srf_pt - srf_pt_ref).norm() >
          Test::Constants::test_tolerance_geometry)
        return 1;
    }
  }
  return 0;
}
