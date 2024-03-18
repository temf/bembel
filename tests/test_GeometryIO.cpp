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

bool compareDATFiles(const std::string& file_name1,
                     const std::string& file_name2) {
  std::ifstream file1(file_name1);
  std::ifstream file2(file_name2);

  if (!file1.is_open() || !file2.is_open()) {
    std::cerr << "Error opening files." << std::endl;
    return false;
  }

  std::string line1, line2;

  while (std::getline(file1, line1) && std::getline(file2, line2)) {
    if (line1.compare(line2) != 0) {
      return false;
    }
  }

  // Check if one file has more lines than the other
  if (std::getline(file1, line1) || std::getline(file2, line2)) {
    return false;
  }

  return true;
}

int main() {
  using namespace Bembel;

  Test::TestGeometryWriter::writePatchesWithDifferentDegree();

  Geometry geometry("patches_with_different_degree.dat");
  assert(geometry.get_geometry().size() == 2);

  const PatchVector& patches = geometry.get_geometry();

  BEMBEL_TEST_IF(patches[0].polynomial_degree_x_ == 5);
  BEMBEL_TEST_IF(patches[0].polynomial_degree_y_ == 6);

  WriteDATFile(patches, "export.dat");

  BEMBEL_TEST_IF(
      compareDATFiles("patches_with_different_degree.dat", "export.dat"));

  return 0;
}
