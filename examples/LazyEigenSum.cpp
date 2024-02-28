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

#include <Bembel/AnsatzSpace>
#include <Bembel/Geometry>
#include <Bembel/H2Matrix>
#include <Bembel/Helmholtz>
#include <Bembel/IO>
#include <iostream>

template <typename Matrix1, typename Matrix2>
void testOpPlus(const Matrix1& m1, const Matrix2& m2) {
  typedef typename Eigen::Matrix<typename Matrix1::Scalar, Eigen::Dynamic, 1>
      VectorType;
  assert(m1.rows() == m2.rows());
  assert(m1.cols() == m2.cols());
  VectorType a = VectorType::Random(m1.cols());
  VectorType b1 = m1 * a;
  VectorType b2 = m2 * a;
  VectorType ref = b1.eval() + b2.eval();
  auto matrix_op = m1 + m2;
  VectorType test = matrix_op * a;
  assert((ref - test).norm() < 1e-14);
}
template <typename Matrix1, typename Matrix2>
void testOpMinus(const Matrix1& m1, const Matrix2& m2) {
  typedef typename Eigen::Matrix<typename Matrix1::Scalar, Eigen::Dynamic, 1>
      VectorType;
  assert(m1.rows() == m2.rows());
  assert(m1.cols() == m2.cols());
  VectorType a = VectorType::Random(m1.cols());
  VectorType b1 = m1 * a;
  VectorType b2 = m2 * a;
  VectorType ref = b1.eval() - b2.eval();
  auto matrix_op = m1 - m2;
  VectorType test = matrix_op * a;
  assert((ref - test).norm() < 1e-14);
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
  typedef H2Matrix<Scalar> H2MatrixType;
  typedef Matrix<Scalar, Dynamic, Dynamic> DenseMatrixType;
  typedef SparseMatrix<Scalar> SparseMatrixType;

  // Build ansatz space
  AnsatzSpace<Operator> ansatz_space(geometry, refinement_level,
                                     polynomial_degree);

  // Set up and compute discrete operator
  DiscreteOperator<H2MatrixType, Operator> disc_op(ansatz_space);
  disc_op.compute();
  const H2MatrixType& H2 = disc_op.get_discrete_operator();

  {  // interaction with itself
    testOpPlus(H2, H2);
    testOpPlus(H2, H2 + H2);
    testOpPlus(H2 + H2, H2);
    testOpPlus(H2, H2 - H2);
    testOpPlus(H2 - H2, H2);
    testOpMinus(H2, H2);
    testOpMinus(H2, H2 + H2);
    testOpMinus(H2 + H2, H2);
    testOpMinus(H2, H2 - H2);
    testOpMinus(H2 - H2, H2);
  }
  {  // interaction with dense matrix
    DenseMatrixType S(H2.rows(), H2.cols());
    S.setIdentity();

    testOpPlus(S, H2);
    testOpPlus(S, H2 + H2);
    testOpPlus(S, H2 - H2);
    testOpPlus(H2, S);
    testOpPlus(H2 + H2, S);
    testOpPlus(H2 - H2, S);
    testOpMinus(S, H2);
    testOpMinus(S, H2 + H2);
    testOpMinus(S, H2 - H2);
    testOpMinus(H2, S);
    testOpMinus(H2 + H2, S);
    testOpMinus(H2 - H2, S);
  }
  {  // interaction with sparse matrix
    SparseMatrixType S(H2.rows(), H2.cols());
    S.setIdentity();

    testOpPlus(S, H2);
    testOpPlus(S, H2 + H2);
    testOpPlus(S, H2 - H2);
    testOpPlus(H2, S);
    testOpPlus(H2 + H2, S);
    testOpPlus(H2 - H2, S);
    testOpMinus(S, H2);
    testOpMinus(S, H2 + H2);
    testOpMinus(S, H2 - H2);
    testOpMinus(H2, S);
    testOpMinus(H2 + H2, S);
    testOpMinus(H2 - H2, S);
  }
  {  // interaction with scalar
    Scalar b = 3.;

    testOpPlus(H2, b * H2);
    testOpPlus(H2, H2 + b * H2);
    testOpPlus(H2, H2 - b * H2);
    testOpPlus(H2, b * H2 + H2);
    testOpPlus(H2, b * H2 - H2);
  }

  return 0;
}
