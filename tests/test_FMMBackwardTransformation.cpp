// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.

#include <Bembel/AnsatzSpace>
#include <Bembel/DummyOperator>
#include <Bembel/Geometry>

#include "Test.hpp"
#include "TestGeometries.hpp"

int main() {
  using namespace Bembel;
  using namespace Eigen;

  int number_of_points = 9;
  int number_of_points2 = number_of_points * number_of_points;
  int refinement_level = 3;
  int cluster_level = 1;

  // generate geometry
  Test::TestGeometryWriter::writeScreen();
  Geometry geometry("test_Screen.dat");
  AnsatzSpace<DummyOperator> ansatz_space(geometry, refinement_level, 0);

  // generate matrices
  MatrixXd fmm_transfer_matrices =
      H2Multipole::computeTransferMatrices<H2Multipole::ChebychevRoots>(
          number_of_points);
  std::vector<MatrixXd> fmm_moment_matrix1 =
      H2Multipole::Moment2D<H2Multipole::ChebychevRoots, DummyOperator>::
          compute2DMoment(ansatz_space.get_superspace(), cluster_level,
                          refinement_level - cluster_level, number_of_points);
  std::vector<MatrixXd> fmm_moment_matrix2 = H2Multipole::
      Moment2D<H2Multipole::ChebychevRoots, DummyOperator>::compute2DMoment(
          ansatz_space.get_superspace(), cluster_level - 1,
          refinement_level - cluster_level + 1, number_of_points);

  fmm_moment_matrix1[0].conservativeResize(2 * number_of_points2,
                                           fmm_moment_matrix1[0].cols());
  fmm_moment_matrix1[0].block(number_of_points2, 0, number_of_points2,
                              fmm_moment_matrix1[0].cols()) =
      fmm_moment_matrix1[0].block(0, 0, number_of_points2,
                                  fmm_moment_matrix1[0].cols());
  fmm_moment_matrix2[0].conservativeResize(2 * number_of_points2,
                                           fmm_moment_matrix2[0].cols());
  fmm_moment_matrix2[0].block(number_of_points2, 0, number_of_points2,
                              fmm_moment_matrix2[0].cols()) =
      fmm_moment_matrix2[0].block(0, 0, number_of_points2,
                                  fmm_moment_matrix2[0].cols());

  // generate random vector
  VectorXd random_vector = VectorXd::Random(2 * number_of_points2);

  // apply moment matrices to these matrices
  MatrixXd testcase = fmm_moment_matrix2[0].transpose() * random_vector;

  // do backward transformation
  std::vector<MatrixXd> backward_dst;
  backward_dst.push_back(MatrixXd::Zero(2 * number_of_points2, 4));
  backward_dst.push_back(random_vector);
  MatrixXd backward = H2Multipole::backwardTransformation(
      fmm_moment_matrix1[0], fmm_transfer_matrices, 1, backward_dst);

  BEMBEL_TEST_IF((testcase - backward).norm() < 1e-12);
}
