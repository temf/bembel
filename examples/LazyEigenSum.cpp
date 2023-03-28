// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.

#include <Bembel/AnsatzSpace>
#include <Bembel/Geometry>
#include <Bembel/H2Matrix>
#include <Bembel/IO>
#include <Bembel/Laplace>
#include <iostream>

int main() {
  using namespace Bembel;
  using namespace Eigen;

  Bembel::IO::Stopwatch sw;

  int polynomial_degree_max = 3;
  int refinement_level_max = 3;

  // Load geometry from file "sphere.dat", which must be placed in the same
  // directory as the executable
  Geometry geometry("sphere.dat");

  // Iterate over polynomial degree.
  for (int polynomial_degree = 0; polynomial_degree < polynomial_degree_max + 1;
       ++polynomial_degree) {
    // Iterate over refinement levels
    for (int refinement_level = 0; refinement_level < refinement_level_max + 1;
         ++refinement_level) {
      sw.tic();
      std::cout << "Degree " << polynomial_degree << " Level "
                << refinement_level;
      // Build ansatz space
      AnsatzSpace<LaplaceSingleLayerOperator> ansatz_space(
          geometry, refinement_level, polynomial_degree);

      // Set up and compute discrete operator
      DiscreteOperator<H2Matrix<double>, LaplaceSingleLayerOperator> disc_op(
          ansatz_space);
      disc_op.compute();
      const H2Matrix<double> &H2 = disc_op.get_discrete_operator();

      // reference random vector
      VectorXd a = VectorXd::Random(H2.cols());

      // H2matrix vector
      auto b1 = H2 * a;
      VectorXd c1 = a + b1;

      // set up dense identity matrix
      MatrixXd S(H2.rows(), H2.cols());
      S.setIdentity();

      // set up sum of matrices
      auto sum = S + H2;

      SparseMatrix<double> gaga(H2.rows(), H2.cols());
      auto sum2 = S + gaga;

      auto b2 = sum * a;
      VectorXd c2 = b2;
      std::cout << (c1 - c2).norm() << std::endl;
    }

    std::cout << std::endl;
  }

  return 0;
}
