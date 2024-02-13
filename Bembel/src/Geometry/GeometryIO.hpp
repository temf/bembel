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

#ifndef BEMBEL_SRC_GEOMETRY_GEOMETRYIO_HPP_
#define BEMBEL_SRC_GEOMETRY_GEOMETRYIO_HPP_

namespace Bembel {

/**
 * \ingroup Geometry
 * \brief loads geometry from file with GEOPDE-format. Note that the direction
 *        of the normals must be consistent
 * \param name path/filename pointing to the geometry file
 * \return std::vector of NURBS::Patch describing geometry
 */
inline std::vector<Patch> LoadGeometryFileDAT(
    const std::string &file_name) noexcept {
  std::vector<Bembel::Patch> out;
  std::stringstream iss;
  std::string word;
  std::string row;
  int infoInt[5];
  std::ifstream file;
  file.open(file_name);
  if (!file) {
    std::cerr << "File " << file_name << " doesn't exist!";
    exit(1);
  }
  // first 4 rows are irrelevant
  for (int i = 0; i < 5; i++) {
    getline(file, row);
  }
  // main information
  iss.str(row);
  for (int i = 0; i < 5; i++) {
    iss >> word;

    infoInt[i] = stoi(word);
  }
  // main loop - patches
  for (int patchNr = 1; patchNr <= infoInt[2]; patchNr++) {
    Bembel::Patch tempPatch;
    std::vector<double> tempknt1;
    std::vector<double> tempknt2;
    std::vector<int> info;  // p and ncp / 0,1-p  2,3-ncp
    std::vector<Eigen::MatrixXd> tmp;

    getline(file, row);  // file_name

    getline(file, row);
    iss.str(row);
    while (iss >> word) {
      info.push_back(stoi(word));
    }
    getline(file, row);

    word = "";
    iss.clear();
    iss.str(row);
    while (iss >> word) {
      info.push_back(stoi(word));
    }

    word = "";
    iss.clear();

    // In textfiles are only 2 knotVectors and 4 Matrices

    getline(file, row);
    iss.str(row);

    // first knotVector

    for (int k = 0; k < info[0] + info[2] + 1; k++) {
      iss >> word;
      tempknt1.push_back(atof(word.c_str()));
    }
    iss.clear();
    getline(file, row);
    iss.str(row);
    word = "";
    // second knotVector
    for (int k = 0; k < info[1] + info[3] + 1; k++) {
      iss >> word;
      tempknt2.push_back(atof(word.c_str()));
    }
    iss.clear();
    // 4 Matrices
    int N = info[2];
    int M = info[3];
    for (int k = 0; k < 4; k++) {
      Eigen::Matrix<double, -1, -1> tempMatrix(
          M, N);  // == Eigen::MatrixXd tempMatrix(M,N);
      getline(file, row);
      iss.str(row);
      for (int i = 0; i < M; i++)
        for (int j = 0; j < N; j++) {
          iss >> word;

          tempMatrix(i, j) = atof(word.c_str());
        }

      tmp.push_back(tempMatrix);
      iss.clear();
    }
    // Important
    tempPatch.init_Patch(tmp, tempknt1, tempknt2);
    out.push_back(tempPatch);
  }

  file.close();
  return out;
}

/**
 * \ingroup Gemetry
 * \brief method to generate textfile for the geometry
 * \param name name of new textfile
 * \param patchnumber overall number of patches
 **/
inline void MakeFile(const std::string &file_name,
                     int number_of_patches) noexcept {
  std::ofstream file;
  file.open(file_name);
  file << "# nurbs mesh v.2.1"
       << "\r\n";
  file << "# " << file_name << "\r\n";
  file << "# Generated by BEMBEL, see www.bembel.eu"
       << "\r\n";
  file << "#"
       << "\r\n";
  file << "2 3 " << number_of_patches << " 0 0 "
       << "\r\n";
  file.close();
}
/**
 * \ingroup Gemetry
 * \brief method to write Patch information into textfile
 * \param knt1 knotVector1
 * \param knt2 knotVector2
 * \param tmp Vector with x,y,z,w Matices
 * \param name filename
 * \param patchnumberCurr current patch number
 **/
void WritePatch(const std::string &file_name, int current_patch_number,
                const std::vector<Eigen::MatrixXd> &xyzw,
                const std::vector<double> &knt1,
                const std::vector<double> &knt2) noexcept {
  std::ofstream file(file_name, std::ios::app);
  int N = xyzw[0].cols();  // ncp
  int M = xyzw[0].rows();  // ncp
  int p1 = knt1.size() - N - 1;
  int p2 = knt2.size() - M - 1;
  file << "PATCH " << current_patch_number << " \r\n";
  file << p1 << " " << p2 << " \r\n";
  file << N << " " << M << " \r\n";
  // knotVectors
  for (unsigned int i = 0; i < knt1.size() - 1; i++) {
    file << std::fixed << std::setprecision(15) << knt1[i] << "   ";
  }
  file << std::fixed << std::setprecision(15) << knt1[knt1.size() - 1]
       << "   \r\n";
  for (unsigned int i = 0; i < knt2.size() - 1; i++) {
    file << std::fixed << std::setprecision(15) << knt2[i] << "   ";
  }
  file << std::fixed << std::setprecision(15) << knt2[knt2.size() - 1]
       << "   \r\n";

  // Matrices
  for (unsigned int n = 0; n < xyzw.size(); n++) {
    for (int i = 0; i < M; i++)
      for (int j = 0; j < N; j++) {
        file << std::fixed << std::setprecision(15) << xyzw[n](i, j);
        if (N * i + j == M * N - 1) {
          file << "   \r\n";
          continue;
        }
        file << "   ";
      }
  }

  file.close();
}

}  // namespace Bembel
#endif  // BEMBEL_SRC_GEOMETRY_GEOMETRYIO_HPP_
