// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.

#include <chrono>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <string>
#include <utility>
#include <vector>
#include <Bembel/Geometry>
#include <Bembel/src/ClusterTree/ElementTreeNode.hpp>
#include <Bembel/src/ClusterTree/ElementTree.hpp>
#include <Bembel/src/IO/Stopwatch.hpp>
#include <iostream>

#include "writeVTK.hpp"

int main() {
  Bembel::Geometry geometry("sphere.dat");
  std::cout << "The geometry has " << geometry.get_number_of_patches()
            << " patches." << std::endl;
  Bembel::IO::Stopwatch T;
  T.tic();
  Bembel::ElementTree et(geometry, 10);
  std::cout << "time: " << T.toc() << " s.\n"; 
  std::cout << "got element tree\n";
  Eigen::MatrixXd P = et.generatePointList();
  Eigen::MatrixXi E = et.generateElementList();
  Eigen::VectorXd z = P.row(2);
  writeMesh2vtk("geometry.vtk", P, E, z);
  return 0;
}
