// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef BEMBEL_SPLINE_BEZIEREXTRACTION_H_
#define BEMBEL_SPLINE_BEZIEREXTRACTION_H_

namespace Bembel {
namespace Spl {

/**
 *  \ingroup Spline
 *  \brief implements Bezier extraction via the solution of an interpolation
 *         problem.
 */
template <typename T>
Eigen::SparseMatrix<T> MakeProjection(const std::vector<T> &x_knots,
                                      const std::vector<T> &y_knots,
                                      const std::vector<T> &x_unique_knots,
                                      const std::vector<T> &y_unique_knots,
                                      const int polynomial_degree_x,
                                      const int polynomial_degree_y) noexcept {
  // In here phi refers to the BIG bezier basis, and psi refers to the small
  // B-Spline basis given by x_knots,y_knots and xp,yp.
  // Sparse does weird things. We need not be that accurate, since the right
  // coeefs will be >=.1
  constexpr double sparseTol = 0.001;
  // Number of functions in X,Y and total
  const int size_phi_x = (x_unique_knots.size() - 1) * polynomial_degree_x;
  const int size_phi_y = (y_unique_knots.size() - 1) * polynomial_degree_y;
  const int size_phi = size_phi_x * size_phi_y;
  const int size_psi_x = x_knots.size() - polynomial_degree_x;
  const int size_psi_y = y_knots.size() - polynomial_degree_y;
  const int size_psi = size_psi_x * size_psi_y;

  // Interpolation points
  const auto xmask = MakeInterpolationMask(polynomial_degree_x);
  const auto ymaks = MakeInterpolationMask(polynomial_degree_y);

  const auto xpoint = MakeInterpolationPoints(x_unique_knots, xmask);
  const auto ypoint = MakeInterpolationPoints(y_unique_knots, ymaks);

  assert(((int)(xmask.size() * (x_unique_knots.size() - 1))) == size_phi_x);
  assert(((int)(xpoint.size())) == size_phi_x);
  assert(((int)(ypoint.size())) == size_phi_y);

  Eigen::Matrix<T, -1, -1> tempsolver =
      Eigen::Matrix<T, -1, -1>::Zero(size_phi, size_phi);

  double *coefficients_x = new double[polynomial_degree_x];
  double *coefficients_y = new double[polynomial_degree_y];

  for (int i = 0; i < polynomial_degree_x; i++) coefficients_x[i] = 0;

  for (int i = 0; i < polynomial_degree_y; i++) coefficients_y[i] = 0;

  auto BezBasis = [](std::vector<double> uniq, int deg, int pos, double pt,
                     double *coef) {
    std::div_t loc = std::div(pos, deg);
    if (pt > uniq[loc.quot + 1] || pt < uniq[loc.quot]) {
      return 0.;
    } else {
      coef[loc.rem] = 1;
      double out = Bembel::Basis::ShapeFunctionHandler::evalCoef(
          deg - 1, coef, Rescale(pt, uniq[loc.quot], uniq[loc.quot + 1]));
      coef[loc.rem] = 0;
      return out;
    }
  };

  // A matrix, each row corresponds to one tp-basisfunction, each col to one
  // tp-interpolation point
  for (int iy = 0; iy < size_phi_y; iy++) {
    for (int ix = 0; ix < size_phi_x; ix++) {
      for (int y = 0; y < size_phi_y; y++) {
        for (int x = 0; x < size_phi_x; x++) {
          tempsolver(y * size_phi_x + x, ix * size_phi_y + iy) =
              BezBasis(x_unique_knots, polynomial_degree_x, ix, xpoint[x],
                       coefficients_x) *
              BezBasis(y_unique_knots, polynomial_degree_y, iy, ypoint[y],
                       coefficients_y);
        }
      }
    }
  }

  delete[] coefficients_x;
  delete[] coefficients_y;

  Eigen::FullPivLU<Eigen::Matrix<double, -1, -1>> lu_decomp(tempsolver);
  assert(lu_decomp.rank() == size_phi);
  /// The inverse matrix, for the solition of the linear system to come.
  Eigen::Matrix<T, -1, -1> solve = tempsolver.inverse();

  Eigen::Matrix<T, -1, -1> Proj(size_phi, size_psi);
  Eigen::Matrix<T, -1, -1> unit =
      Eigen::Matrix<T, -1, -1>::Zero(size_psi_y, size_psi_x);

  /// For each basisfunction phi_j solves for the right coeef such that \Sum
  for (int i = 0; i < size_psi_x; i++)
    for (int j = 0; j < size_psi_y; j++) {
      unit(j, i) = 1.;
      Eigen::Matrix<T, -1, 1> smallBase =
          (Unroll(DeBoorTP(unit, x_knots, y_knots, xpoint, ypoint)));
      Proj.col(i * size_psi_y + j) = (solve * smallBase).transpose();
      unit(j, i) = 0.;
    }

  return (Proj.sparseView()).pruned(sparseTol);
}
}  // namespace Spl
}  // namespace Bembel
#endif
