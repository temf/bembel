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

#include <Bembel/AnsatzSpace>
#include <Bembel/Geometry>
#include <Bembel/H2Matrix>
#include <Bembel/Helmholtz>
#include <Bembel/IO>
#include <iostream>

template <typename Matrix1, typename Matrix2>
void testEqual(const Matrix1& m1, const Matrix2& m2, const double tol) {
  typedef typename Eigen::Matrix<typename Matrix1::Scalar, Eigen::Dynamic, 1>
      VectorType;
  assert(m1.rows() == m2.rows());
  assert(m1.cols() == m2.cols());
  VectorType a = VectorType::Random(m1.cols());
  VectorType b1 = m1 * a;
  VectorType b2 = m2 * a;
  std::cout << "Error is " << (b1 - b2).norm() << std::endl;
  assert((b1 - b2).norm() < tol);
}

int main() {
  using namespace Bembel;
  using namespace Eigen;

  int polynomial_degree = 1;
  int refinement_level = 2;

  // Load geometry from file "sphere.dat", which must be placed in the same
  // directory as the executable
  Geometry geometry("sphere.dat");

  // define types
  typedef HelmholtzSingleLayerOperator Operator;
  typedef typename LinearOperatorTraits<Operator>::Scalar Scalar;
  typedef typename LinearOperatorTraits<Operator>::EigenType VectorType;
  typedef Matrix<Scalar, Dynamic, Dynamic> DenseMatrixType;
  typedef H2Matrix<Scalar> H2MatrixType;
  typedef Matrix<Scalar, Dynamic, Dynamic> DenseMatrixType;
  typedef SparseMatrix<Scalar> SparseMatrixType;

  // Build ansatz space
  AnsatzSpace<Operator> ansatz_space(geometry, refinement_level,
                                     polynomial_degree);

  // Set up and compute discrete operator
  DiscreteOperator<DenseMatrixType, Operator> disc_op_dense(ansatz_space);
  disc_op_dense.compute();
  const DenseMatrixType& D = disc_op_dense.get_discrete_operator();
  DiscreteOperator<H2MatrixType, Operator> disc_op_H2(ansatz_space);
  disc_op_H2.compute();
  const H2MatrixType& H2 = disc_op_H2.get_discrete_operator();

  {  // interaction with itself
    double tol = 1e-8;

    testEqual(D + D, H2 + H2, tol);
    testEqual(D - D, H2 - H2, tol);
    testEqual(D + D + D, H2 + H2 + H2, tol);
    testEqual(D - D + D, H2 - H2 + H2, tol);
    testEqual(D + D - D, H2 + H2 - H2, tol);
    testEqual(D - D - D, H2 - H2 - H2, tol);
  }
  {  // interaction with dense matrix
    double tol = 1e-8;

    testEqual(D + D, D + H2, tol);
    testEqual(D + D, H2 + D, tol);
    testEqual(D - D, D - H2, tol);
    testEqual(D - D, H2 - D, tol);
    testEqual(D + D + D, D + H2 + H2, tol);
    testEqual(D + D + D, H2 + D + H2, tol);
    testEqual(D + D + D, H2 + H2 + D, tol);
    testEqual(D - D + D, D - H2 + H2, tol);
    testEqual(D - D + D, H2 - D + H2, tol);
    testEqual(D - D + D, H2 - H2 + D, tol);
    testEqual(D + D - D, D + H2 - H2, tol);
    testEqual(D + D - D, H2 + D - H2, tol);
    testEqual(D + D - D, H2 + H2 - D, tol);
    testEqual(D - D - D, D - H2 - H2, tol);
    testEqual(D - D - D, H2 - D - H2, tol);
    testEqual(D - D - D, H2 - H2 - D, tol);
  }
  {  // interaction with sparse matrix
    double tol = 1e-8;

    SparseMatrixType S(H2.rows(), H2.cols());
    S.setIdentity();

    testEqual(S + D, S + H2, tol);
    testEqual(D + S, H2 + S, tol);
    testEqual(S - D, S - H2, tol);
    testEqual(D - S, H2 - S, tol);
    testEqual(S + D + D, S + H2 + H2, tol);
    testEqual(D + S + D, H2 + S + H2, tol);
    testEqual(D + D + S, H2 + H2 + S, tol);
    testEqual(S - D + D, S - H2 + H2, tol);
    testEqual(D - S + D, H2 - S + H2, tol);
    testEqual(D - D + S, H2 - H2 + S, tol);
    testEqual(S + D - D, S + H2 - H2, tol);
    testEqual(D + S - D, H2 + S - H2, tol);
    testEqual(D + D - S, H2 + H2 - S, tol);
    testEqual(S - D - D, S - H2 - H2, tol);
    testEqual(D - S - D, H2 - S - H2, tol);
    testEqual(D - D - S, H2 - H2 - S, tol);
  }
  {  // intercation with scalar
    double tol = 1e-8;
    Scalar b = 3.;

    testEqual(D, H2, tol);
    testEqual(b * D, b * H2, tol);
    testEqual(b*D+D, b*H2 + H2, tol);
  }

  return 0;
}
