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

#ifndef BEMBEL_SRC_MULTIGRID_MULTIGRIDSOLVER_HPP_
#define BEMBEL_SRC_MULTIGRID_MULTIGRIDSOLVER_HPP_

namespace Bembel {
namespace MG {

#define MULTIGRID_minLvl 1
#define MULTIGRID_cycleType 2
#define MULTIGRID_preSmooth 3
#define MULTIGRID_postSmooth 3

/**
 *  \ingroup MultiGrid
 *  \brief Add Gauss Seidel smoother for use with multigrid.
 **/
template <typename Derived, typename Derived2, typename Derived3>
void GaussSeidelSmoother(const Derived &S, const Eigen::MatrixBase<Derived2> &b,
                         Eigen::MatrixBase<Derived3> *x,
                         const unsigned int steps = 3) {
  for (auto i = 0; i < steps; ++i)
    x->derived() = S.template triangularView<Eigen::Lower>().solve(
        b - S.template triangularView<Eigen::StrictlyUpper>() * x->derived());
}

/**
 *  \ingroup MultiGrid
 *  \brief Implements a single step of multiplicative multigrid.
 **/
template <typename Derived, typename Derived2, typename Derived3,
          typename Derived4>
void multiplicativeMultiGrid(const Derived &Ss,
                             const Eigen::MatrixBase<Derived2> &f,
                             Eigen::MatrixBase<Derived3> *x, const Derived4 &Ps,
                             const unsigned int lvl = MULTIGRID_minLvl) {
  typedef typename Derived3::Scalar Scalar;
  typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Matrix;
  if (lvl <= MULTIGRID_minLvl) {
    Matrix A = Matrix(Ss[lvl]);
    x->derived() = A.fullPivHouseholderQr().solve(f);
  } else {
    // Pre-smooth
    GaussSeidelSmoother(Ss[lvl], f, x, MULTIGRID_preSmooth);
    Matrix r = f - Ss[lvl] * x->derived();
    // Restrict the residual
    r = Ps[lvl - 1].transpose() * r;
    Matrix y = r;
    y.setZero();
    // Recursively calling multiplicativeMultiGrid for computing the error
    for (auto i = 0; i < MULTIGRID_cycleType; ++i)
      multiplicativeMultiGrid(Ss, r, &y, Ps, lvl - 1);
    // Prolongate the error
    y = Ps[lvl - 1] * y;
    // Correction
    x->derived() += y;
    // Post-smooth
    GaussSeidelSmoother(Ss[lvl], f, x, MULTIGRID_postSmooth);
  }
}

}  // namespace MG
}  // namespace Bembel
#endif  // BEMBEL_SRC_MULTIGRID_MULTIGRIDSOLVER_HPP_
