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

#ifndef BEMBEL_SRC_GEOMETRY_GEOMETRYIGS_HPP_
#define BEMBEL_SRC_GEOMETRY_GEOMETRYIGS_HPP_

namespace Bembel {

/**
 * \ingroup Geometry
 * \brief loads  geometry from IGES file. Note that the direction
 *        of the normals must be consistent.
 * \param name path/filename pointing to the geometry file
 * \return std::vector of NURBS::Patch describing geometry
 */
std::vector<Patch> LoadGeometryFileIGS(const std::string &file_name) noexcept {
  std::ifstream file;
  std::vector<int> patch_lines;
  file.open(file_name);
  if (!file) {
    std::cerr << "File " << file_name << " doesn't exist!";
    exit(1);
  }
  std::string row;
  getline(file, row);

  // character 72 denotes the section and the information starts in Section D
  while (row[72] != 'D') {
    getline(file, row);
  }

  // collect two lines info of each patch from Directory section
  while (row[72] != 'P') {
    std::stringstream iss;
    std::string word;
    std::vector<int> info1, info2;
    iss.str(row);
    while (iss >> word) {
      info1.push_back(std::stoi(word));
    }
    assert(info1[0] == 128 && "Entry type must be NURBS surface!");
    // this info1 in not going to be used further

    getline(file, row);
    iss.clear();
    iss.str(row);
    while (iss >> word) {
      info2.push_back(std::stoi(word));
    }
    getline(file, row);
    // This entry denotes how many lines correspond to a patch
    patch_lines.push_back(info2[3]);
  }

  // main loop over patches
  std::vector<Patch> out;
  out.reserve(patch_lines.size());
  const int number_of_patches = patch_lines.size();
  for (auto i = 0; i < number_of_patches; ++i) {
    std::stringstream iss;
    std::string word;
    std::vector<double> data;
    for (auto j = 0; j < patch_lines[i]; ++j) {
      // characters 0 to 63 contain data
      std::string raw_data = row.substr(0, 64);
      raw_data.erase(std::remove(raw_data.begin(), raw_data.end(), ' '),
                     raw_data.end());

      iss.str(raw_data);
      while (std::getline(iss, word, ',')) {
        data.push_back(atof(word.c_str()));
      }
      iss.clear();

      getline(file, row);
    }

    const int K1 = data[1];
    const int K2 = data[2];
    const int M1 = data[3];
    const int M2 = data[4];

    const int N1 = 1 + K1 - M1;
    const int N2 = 1 + K2 - M2;
    const int A = N1 + 2 * M1;
    const int B = N2 + 2 * M2;
    const int C = (1 + K1) * (1 + K2);

    // The data of the first entries is read and not needed any more
    data.erase(data.begin(), data.begin() + 10);

    // the + 1 is necessary because the construct excludes the last iterator
    std::vector<double> tempknt1(data.begin(), data.begin() + A + 1);
    std::vector<double> tempknt2(data.begin() + A + 1,
                                 data.begin() + A + B + 2);

    const int M = K2 + 1;
    const int N = K1 + 1;

    Eigen::MatrixXd weights(M, N);
    auto it_weights = data.begin() + A + B + 2;
    for (auto entry = it_weights; entry != it_weights + C; ++entry) {
      int index = entry - it_weights;
      weights((int)index / N, index % N) = *entry;
    }

    Eigen::MatrixXd x_coordinates(M, N);
    Eigen::MatrixXd y_coordinates(M, N);
    Eigen::MatrixXd z_coordinates(M, N);
    auto it_points = data.begin() + A + B + C + 2;
    for (auto entry = it_points; entry != it_points + 3 * C; entry += 3) {
      int index = (int)(entry - it_points) / 3;
      x_coordinates((int)index / N, index % N) = *entry;
      y_coordinates((int)index / N, index % N) = *(entry + 1);
      z_coordinates((int)index / N, index % N) = *(entry + 2);
    }

    // we need to transfer to homogeneous coordinates
    x_coordinates = x_coordinates.cwiseProduct(weights);
    y_coordinates = y_coordinates.cwiseProduct(weights);
    z_coordinates = z_coordinates.cwiseProduct(weights);

    std::vector<Eigen::MatrixXd> tmp = {x_coordinates, y_coordinates,
                                        z_coordinates, weights};

    Bembel::Patch tempPatch;
    tempPatch.init_Patch(tmp, tempknt1, tempknt2);
    out.push_back(tempPatch);
  }
  file.close();
  return out;
}

}  // namespace Bembel
#endif  // BEMBEL_SRC_GEOMETRY_GEOMETRYIGS_HPP_
