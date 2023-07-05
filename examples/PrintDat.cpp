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

#include "Bembel/src/IO/Stopwatch.hpp"
#include "Bembel/src/IO/print2dat.hpp"

int main() {
  Bembel::IO::Stopwatch S;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> bla =
      Eigen::MatrixXd::Random(10000, 100);
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> bla2;
  S.tic();
  Bembel::IO::print2bin("mat.dat", 2 * bla);
  Bembel::IO::bin2Mat("mat.dat", &bla2);
  std::cout << "error: " << (2 * bla - bla2).norm() << " " << bla.norm()
            << std::endl;
  std::cout << "time: " << S.toc() << std::endl;
  // std::cout << bla << std::endl;
  return 0;
}
