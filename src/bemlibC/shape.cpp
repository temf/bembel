// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#include "spline/shape.h"

namespace Spl {

/**
 * \brief Loads geometry from file with GEOPDE-format. Note that the direction
 * of the normals must be consistent
 *
 * \param[in]  name      Path to geometry file
 *
 * \return std::vector of NURBS::Patch describing geometry
 *
 */
std::vector<Spl::Patch> loadGeometryFile(std::string name) {
  std::vector<Spl::Patch> out;
  std::stringstream iss;
  std::string word;
  std::string row;
  int infoInt[5];
  std::ifstream file;
  file.open(name);
  if (!file) {
    std::cerr << "File " << name << " doesn't exist!";
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
  // std::cout << "number of patches: \n" << infoInt[2] << std::endl;

  // main loop - patches
  for (int patchNr = 1; patchNr <= infoInt[2]; patchNr++) {
    Spl::Patch tempPatch;
    std::vector<double> tempknt1;
    std::vector<double> tempknt2;
    std::vector<int> info;  // p and ncp / 0,1-p  2,3-ncp
    std::vector<Eigen::MatrixXd> tmp;

    getline(file, row);  // name
    // std::cout << "loading :  " << row << std::endl;

    getline(file, row);
    iss.str(row);
    // std::cout << "p :  " << std::endl;
    while (iss >> word) {
      // std::cout << word << std::endl;
      info.push_back(stoi(word));
    }
    getline(file, row);

    word = "";
    iss.clear();
    iss.str(row);
    // std::cout << "ncp :  " << std::endl;
    while (iss >> word) {
      // std::cout << word << std::endl;
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
    // std::cout << "first knotVector length :  " << tempknt1.size() <<
    // std::endl;
    iss.clear();
    getline(file, row);
    iss.str(row);
    word = "";
    // second knotVector
    for (int k = 0; k < info[1] + info[3] + 1; k++) {
      iss >> word;
      tempknt2.push_back(atof(word.c_str()));
    }
    // std::cout << "second knotVector length :  " << tempknt2.size() <<
    // std::endl;
    iss.clear();
    // 4 Matrices
    int N = info[2];
    int M = info[3];
    // std::cout << "Matrix size:  :  " << M << " x " << N << std::endl;
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
    tempPatch.initHom(tmp, tempknt1, tempknt2);
    out.push_back(tempPatch);
    // std::cout << "PATCH " << patchNr << " LOADED" << std::endl;
    // std::cout << "\n" << std::endl;
    //
  }

  file.close();
  return out;
}

/**
 * Method to write Patch information in Textfile
 *[input] knt1 -knotVector1
 *        knt2 -knotVector2
 *        tmp - Vector with x,y,z,w Matices
 *        name - filename
 *        patchNrCurr - current patch number
 *[output] - void but adds information in Textfile
 **/
void writePatch(std::string name, int patchNrCurr,
                std::vector<Eigen::MatrixXd> xyzw, std::vector<double> knt1,
                std::vector<double> knt2) {
  std::ofstream file(name, std::ios::app);
  int N = xyzw[0].cols();  // ncp
  int M = xyzw[0].rows();  // ncp
  int p1 = knt1.size() - N - 1;
  int p2 = knt2.size() - M - 1;
  file << "PATCH " << patchNrCurr << " \r\n";
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

}  // namespace Spl
