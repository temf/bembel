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

bool compareIGESFiles(const std::string& file_name1,
                      const std::string& file_name2) {
  std::ifstream file1(file_name1);
  std::ifstream file2(file_name2);

  if (!file1.is_open() || !file2.is_open()) {
    std::cerr << "Error opening files." << std::endl;
    return false;
  }

  std::string line1, line2;
  int line_count = 0;

  while (std::getline(file1, line1) && std::getline(file2, line2)) {
    // line 5 contains the date and need to be treated different
    if (line_count == 5) {
      if (line1.substr(0, 35).compare(line2.substr(0, 35)) != 0 ||
          line1.substr(54).compare(line2.substr(54)) != 0) {
        return false;
      }
      ++line_count;
      continue;
    }
    if (line1.compare(line2) != 0) {
      return false;
    }

    ++line_count;
  }

  // Check if one file has more lines than the other
  if (std::getline(file1, line1) || std::getline(file2, line2)) {
    return false;
  }

  return true;
}

int main() {
  using namespace Bembel;

  Test::TestGeometryWriter::writeIGSScreen();

  Bembel::Geometry geometry("test_Screen.igs");
  assert(geometry.get_geometry().size() == 2);

  const PatchVector& patches = geometry.get_geometry();

  BEMBEL_TEST_IF(patches[0].polynomial_degree_x_ == 5);
  BEMBEL_TEST_IF(patches[0].polynomial_degree_y_ == 6);

  const int precision = 15;
  writeIGSFile(patches, "test_Screen_Export.igs", precision);

  BEMBEL_TEST_IF(compareIGESFiles("test_Screen.igs", "test_Screen_Export.igs"));

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
