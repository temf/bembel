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

#ifndef BEMBEL_SRC_MULTIGRID_PROLONGATIONBSPLINES_HPP_
#define BEMBEL_SRC_MULTIGRID_PROLONGATIONBSPLINES_HPP_

namespace Bembel {
namespace MG {

template <typename Scalar>
std::vector<Eigen::Triplet<Scalar>> toEigenTriplets(
    Eigen::SparseMatrix<Scalar> &M) {
  std::vector<Eigen::Triplet<Scalar>> v;
  for (int i = 0; i < M.outerSize(); i++)
    for (typename Eigen::SparseMatrix<Scalar>::InnerIterator it(M, i); it; ++it)
      v.emplace_back(it.row(), it.col(), it.value());
  return v;
}

/**
 *  \ingroup MultiGrid
 *  \brief Return the nonzero entries info of the row id of the prolongation
 * matrix using Oslo-Algorithm
 * https://www.uio.no/studier/emner/matnat/ifi/nedlagte-emner/INF-MAT5340/v07/undervisningsmateriale/kap4.pdf
 */
std::tuple<Eigen::MatrixXd, unsigned int> nonZerosInfoRow(
    const std::vector<double> &knots,
    const std::vector<double> &uniform_refined_knots,
    unsigned int polynomial_degree, unsigned int row_id,
    unsigned int guess = 0) {
  unsigned int d = guess;
  while (uniform_refined_knots[row_id] >= knots[d]) ++d;
  --d;
  Eigen::MatrixXd non_zero_alpha =
      Eigen::MatrixXd::Ones(1, polynomial_degree + 1);
  for (auto i = 0; i < polynomial_degree; ++i) {
    Eigen::MatrixXd R = Eigen::MatrixXd::Zero(i + 1, i + 2);
    for (auto j = 0; j < i + 1; ++j) {
      R(j, j) = (knots[d + j + 1] - uniform_refined_knots[row_id + i + 1]) /
                (knots[d + j + 1] - knots[d + j - i]);
      R(j, j + 1) = 1.0 - R(j, j);
    }
    non_zero_alpha.block(0, 0, 1, i + 2) =
        non_zero_alpha.block(0, 0, 1, i + 1) * R;
  }

  return std::make_tuple(non_zero_alpha, d);
}

/**
 *  \ingroup MultiGrid
 *  \brief Return prolongation matrix from  level-1 to level in 1D.
 */
Eigen::SparseMatrix<double> prolongateBsplines1D(
    unsigned int polynomial_degree, unsigned int level,
    unsigned int knotrepetition = 1) {
  assert(level > 0 &&
         "Level must be 1 or larger for build prolongation matrix from level-1 "
         "to level");
  std::vector<double> knots = Bembel::Spl::MakeUniformKnotVector(
      polynomial_degree + 1, ((1 << (level - 1)) - 1), knotrepetition);
  std::vector<double> uniform_refined_knots = Spl::MakeUniformKnotVector(
      polynomial_degree + 1, ((1 << level) - 1), knotrepetition);
  unsigned int num_basis = knots.size() - polynomial_degree - 1;
  unsigned int num_fine_basis =
      uniform_refined_knots.size() - polynomial_degree - 1;
  typedef typename Eigen::Triplet<double> Triplet;
  std::vector<Triplet> tripletList;
  tripletList.reserve(num_fine_basis * (polynomial_degree + 1));
  int d = 0;
  for (int i = 0; i < num_fine_basis; ++i) {
    for (int j = 0; j < polynomial_degree + 1; ++j) {
      Eigen::MatrixXd non_zero_alpha;
      std::tie(non_zero_alpha, d) = nonZerosInfoRow(
          knots, uniform_refined_knots, polynomial_degree, i, d);
      tripletList.push_back(
          Triplet(i, d - polynomial_degree + j, non_zero_alpha(j)));
    }
  }
  Eigen::SparseMatrix<double> P;
  P.resize(num_fine_basis, num_basis);
  P.setFromTriplets(tripletList.begin(), tripletList.end());
  return P;
}

/**
 *  \ingroup MultiGrid
 *  \brief Return prolongation matrix from  level-1 to level in 2D.
 */
Eigen::SparseMatrix<double> prolongateBsplines2D(
    unsigned int polynomial_degree, unsigned int level,
    unsigned int knotrepetition = 1) {
  assert(level > 0 &&
         "Level must be 1 or larger for build prolongation matrix from level-1 "
         "to level");
  Eigen::SparseMatrix<double> P_ =
      prolongateBsplines1D(polynomial_degree, level, knotrepetition);
  Eigen::SparseMatrix<double> P = Eigen::kroneckerProduct(P_, P_);
  return P;
}

/**
 *  \ingroup MultiGrid
 *  \brief Return prolongation matrix from  level-1 to level in 2D in  case of
 * multiple patches.
 */
Eigen::SparseMatrix<double> prolongateBsplinesMultiPatches2D(
    unsigned int polynomial_degree, unsigned int level,
    unsigned int num_patches, unsigned int knotrepetition = 1) {
  assert(level > 0 &&
         "Level must be 1 or larger for build prolongation matrix from level-1 "
         "to level");
  Eigen::SparseMatrix<double> P_ =
      prolongateBsplines2D(polynomial_degree, level, knotrepetition);
  std::vector<Eigen::Triplet<double>> triplets_ = toEigenTriplets(P_);
  std::vector<Eigen::Triplet<double>> triplets;
  Eigen::SparseMatrix<double> P(num_patches * P_.rows(),
                                num_patches * P_.cols());
  for (auto i = 0; i < num_patches; ++i) {
    for (auto entry : triplets_) {
      triplets.push_back(Eigen::Triplet<double>(entry.row() + i * P_.rows(),
                                                entry.col() + i * P_.cols(),
                                                entry.value()));
    }
  }
  P.setFromTriplets(triplets.begin(), triplets.end());
  return P;
}

}  // namespace MG
}  // namespace Bembel
#endif  // BEMBEL_SRC_MULTIGRID_PROLONGATIONBSPLINES_HPP_
