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

int main() {
  using namespace Bembel;

  Test::TestGeometryWriter::writePatchInternalKnots();

  Bembel::Geometry geometry("patch_internal_knots.dat");
  BEMBEL_TEST_IF(geometry.get_geometry().size() == 4);

  return 0;
}
