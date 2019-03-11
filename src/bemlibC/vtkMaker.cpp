// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#include <iostream>
#include "Spline.hpp"
#include "Visualize.hpp"

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cout << "You need to pass a .dat-geometry as input.\n";
    return 0;
  }
  auto geom = Spl::loadGeometryFile(argv[1]);
  Bembel::Vis::plot_geometry_only(geom, 5, "geometry.vtk");
  return 0;
}
