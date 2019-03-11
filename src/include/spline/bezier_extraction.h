// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef _SPL_PROJECTOR_
#define _SPL_PROJECTOR_

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include "spline/basis.h"
#include "spline/deBoorTP.h"
#include "spline/knots.h"
#include "spline/localise.h"

namespace Spl {
// #include <Eigen/OrderingMethods>
// #include <Eigen/SparseLU>

/**
 *  @brief Implements Bezier extraction via the solution of an interpolation
 * problem.
 *
 */
template <typename T>
Eigen::SparseMatrix<T> make_projection(const std::vector<T> &xknt,
                                       const std::vector<T> &yknt,
                                       const std::vector<T> &xuniq,
                                       const std::vector<T> &yuniq,
                                       const int xp, const int yp) noexcept {
  // In here Phi refers to the BIG bezier basis, and psi refers to the small
  // B-Spline basis given by xknt,yknt and xp,yp.

  // Sparse does weird things. We need not be that accurate, since the right
  // coeefs will be >=.1
  constexpr double sparseTol = 0.001;
  // Number of functions in X,Y and total
  const int nPhiX = (xuniq.size() - 1) * xp;
  const int nPhiY = (yuniq.size() - 1) * yp;
  const int nPhi = nPhiX * nPhiY;
  const int npsix = xknt.size() - xp;
  const int npsiy = yknt.size() - yp;
  const int npsi = npsix * npsiy;

  // Interpolation points
  const auto xmask = make_interpolation_mask(xp);
  const auto ymaks = make_interpolation_mask(yp);

  const auto xpoint = make_interpolation_points(xuniq, xmask);
  const auto ypoint = make_interpolation_points(yuniq, ymaks);

  assert(((int)(xmask.size() * (xuniq.size() - 1))) == nPhiX);
  assert(((int)(xpoint.size())) == nPhiX);
  assert(((int)(ypoint.size())) == nPhiY);

  Eigen::Matrix<T, -1, -1> tempsolver =
      Eigen::Matrix<T, -1, -1>::Zero(nPhi, nPhi);

  double *xcoef = new double[xp];
  double *ycoef = new double[yp];

  for (int i = 0; i < xp; i++) xcoef[i] = 0;

  for (int i = 0; i < yp; i++) ycoef[i] = 0;

  auto bezBasis = [](std::vector<double> uniq, int deg, int pos, double pt,
                     double *coef) {
    std::div_t loc = std::div(pos, deg);
    if (pt > uniq[loc.quot + 1] || pt < uniq[loc.quot]) {
      return 0.;
    } else {
      coef[loc.rem] = 1;
      double out = evalBrnstn(deg - 1, coef,
                              rescale(pt, uniq[loc.quot], uniq[loc.quot + 1]));
      coef[loc.rem] = 0;
      return out;
    }
  };

  // A matrix, each row corresponds to one tp-basisfunction, each col to one
  // tp-interpolation point

  for (int iy = 0; iy < nPhiY; iy++) {
    for (int ix = 0; ix < nPhiX; ix++) {
      for (int y = 0; y < nPhiY; y++) {
        for (int x = 0; x < nPhiX; x++) {
          tempsolver(y * nPhiX + x, ix * nPhiY + iy) =
              bezBasis(xuniq, xp, ix, xpoint[x], xcoef) *
              bezBasis(yuniq, yp, iy, ypoint[y], ycoef);
        }
      }
    }
  }

  delete[] xcoef;
  delete[] ycoef;

  Eigen::FullPivLU<Eigen::Matrix<double, -1, -1>> lu_decomp(tempsolver);
  assert(lu_decomp.rank() == nPhi);
  // The inverse matrix, for the solition of the linear system to come.
  Eigen::Matrix<T, -1, -1> solve = tempsolver.inverse();

  Eigen::Matrix<T, -1, -1> Proj(nPhi, npsi);
  Eigen::Matrix<T, -1, -1> unit = Eigen::Matrix<T, -1, -1>::Zero(npsiy, npsix);

  // For each basisfunction phi_j solves for the right coeef such that \Sum
  // c_k Phi_k = phi_i.
  for (int i = 0; i < npsix; i++)
    for (int j = 0; j < npsiy; j++) {
      unit(j, i) = 1.;
      Eigen::Matrix<T, -1, 1> smallBase =
          (unroll(deBoorTP(unit, xknt, yknt, xpoint, ypoint)));
      Proj.col(i * npsiy + j) = (solve * smallBase).transpose();
      unit(j, i) = 0.;
    }

  // std::cout << nPhi << "  " << npsi << "\n";

  return (Proj.sparseView()).pruned(sparseTol);
}
}  // namespace Spl
#endif
