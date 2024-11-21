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

#ifndef BEMBEL_SRC_MULTIGRID_PROLONGATION_HPP_
#define BEMBEL_SRC_MULTIGRID_PROLONGATION_HPP_

namespace Bembel {
namespace MG {

/**
 *    \ingroup MultiGrid
 *    \brief Return the prolongation matrix from level l-1 to level l.
 **/
template <typename Derived>
Eigen::SparseMatrix<double> prolongationMatrix(
    const AnsatzSpace<Derived> &ansatz_space, unsigned int level) {
  assert(level > 0 &&
         "Level must be 1 or larger for build prolongation matrix from level-1 "
         "to level");
  unsigned int knot_repetition = ansatz_space.get_knot_repetition();
  unsigned int polynomial_degree = ansatz_space.get_polynomial_degree();
  unsigned int num_patches = ansatz_space.get_number_of_patches();

  SuperSpace<Derived> coarse_super_space;
  coarse_super_space.init_SuperSpace(
      Geometry(ansatz_space.get_superspace().get_geometry()), level - 1,
      polynomial_degree);
  Projector<Derived> coarse_proj(coarse_super_space, knot_repetition);
  Glue<Derived> coarse_glue(coarse_super_space, coarse_proj);

  SuperSpace<Derived> fine_super_space;
  fine_super_space.init_SuperSpace(
      Geometry(ansatz_space.get_superspace().get_geometry()), level,
      polynomial_degree);
  Projector<Derived> fine_proj(fine_super_space, knot_repetition);
  Glue<Derived> fine_glue(fine_super_space, fine_proj);

  Eigen::VectorXd degree_vector =
      fine_glue.get_glue_matrix().transpose() *
      Eigen::VectorXd::Ones(fine_glue.get_glue_matrix().rows());
  // the reason we multiply 0.5 in the begining is because we use the scaled
  // B-spline basis function in the Bembel
  return 0.5 * degree_vector.cwiseInverse().asDiagonal() *
         fine_glue.get_glue_matrix().transpose() *
         MG::prolongateBsplinesMultiPatches2D(polynomial_degree, level,
                                              num_patches, knot_repetition) *
         coarse_glue.get_glue_matrix();
}

}  // namespace MG
}  // namespace Bembel
#endif  // BEMBEL_SRC_MULTIGRID_PROLONGATION_HPP_
