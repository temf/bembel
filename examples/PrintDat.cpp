#include "Bembel/src/IO/print2dat.hpp"
#include "Bembel/src/IO/Stopwatch.hpp"

int main() {
  Bembel::IO::Stopwatch S;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> bla = Eigen::MatrixXd::Random(10000, 100);
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>  bla2;
  S.tic();
  Bembel::IO::print2bin("mat.dat", 2 * bla);
  Bembel::IO::bin2Mat("mat.dat", &bla2);
  std::cout << "error: " << (2 * bla - bla2).norm() << " " << bla.norm() <<  std::endl;
  std::cout << "time: " << S.toc() << std::endl;
  //std::cout << bla << std::endl;
  return 0;
}
