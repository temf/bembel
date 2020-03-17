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
  int refinement_level = 6;
  int cluster_level = 2;

  // generate geometry
  Test::TestGeometryWriter::writeScreen();
  Geometry geometry("test_Screen.dat");
  AnsatzSpace<DummyOperator> ansatz_space(geometry, refinement_level, 0);

  // generate matrices
  MatrixXd fmm_transfer_matrices =
      H2Multipole::computeTransferMatrices<H2Multipole::ChebychevRoots>(
          number_of_points);
  std::vector<MatrixXd> fmm_moment_matrix0 = H2Multipole::Moment2D<
      Bembel::H2Multipole::ChebychevRoots,
      DummyOperator>::compute2DMoment(ansatz_space.get_superspace(),
                                      cluster_level,
                                      refinement_level - cluster_level,
                                      number_of_points);
  std::vector<MatrixXd> fmm_moment_matrix1 = H2Multipole::Moment2D<
      Bembel::H2Multipole::ChebychevRoots,
      DummyOperator>::compute2DMoment(ansatz_space.get_superspace(),
                                      cluster_level - 1,
                                      refinement_level - cluster_level + 1,
                                      number_of_points);
  std::vector<MatrixXd> fmm_moment_matrix2 = H2Multipole::Moment2D<
      Bembel::H2Multipole::ChebychevRoots,
      DummyOperator>::compute2DMoment(ansatz_space.get_superspace(),
                                      cluster_level - 2,
                                      refinement_level - cluster_level + 2,
                                      number_of_points);

  fmm_moment_matrix0[0].conservativeResize(2 * number_of_points2,
                                           fmm_moment_matrix0[0].cols());
  fmm_moment_matrix0[0].block(number_of_points2, 0, number_of_points2,
                              fmm_moment_matrix0[0].cols()) =
      fmm_moment_matrix0[0].block(0, 0, number_of_points2,
                                  fmm_moment_matrix0[0].cols());

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

  // generate random vector and reshape to matrices
  VectorXd random_vector = VectorXd::Random(fmm_moment_matrix2[0].cols());
  MatrixXd random_matrix0 =
      Map<MatrixXd>(random_vector.data(), fmm_moment_matrix0[0].cols(),
                    random_vector.rows() / fmm_moment_matrix0[0].cols());
  MatrixXd random_matrix1 =
      Map<MatrixXd>(random_vector.data(), fmm_moment_matrix1[0].cols(),
                    random_vector.rows() / fmm_moment_matrix1[0].cols());
  MatrixXd random_matrix2 =
      Map<MatrixXd>(random_vector.data(), fmm_moment_matrix2[0].cols(),
                    random_vector.rows() / fmm_moment_matrix2[0].cols());

  // apply moment matrices to these matrices
  MatrixXd testcase0 = fmm_moment_matrix0[0] * random_matrix0;
  MatrixXd testcase1 = fmm_moment_matrix1[0] * random_matrix1;
  MatrixXd testcase2 = fmm_moment_matrix2[0] * random_matrix2;

  // do forward transformation
  std::vector<MatrixXd> forward = H2Multipole::forwardTransformation(
      fmm_moment_matrix0[0], fmm_transfer_matrices, refinement_level - 1 - 1,
      random_matrix0);

  BEMBEL_TEST_IF(((testcase0 - forward[0]).norm() +
                  (testcase1 - forward[1]).norm() +
                  (testcase2 - forward[2]).norm()) < 1e-12);
}
