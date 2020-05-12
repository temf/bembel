// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef BEMBEL_SPLINE_LOCALISE_H_
#define BEMBEL_SPLINE_LOCALISE_H_

namespace Bembel {

/**
 *  \ingroup Spline
 *  \brief The Spl namespace is the backend for basis functions and B-Spline
 * computations, and should not be used by the user in a boundary element code.
 */
namespace Spl {
// At first the a==0 and b==1 seems insane. However, this is the case quite
// often and actually shaves off a couple of % of runtime, since the comparison
// is cheaper than a division.
constexpr inline double Rescale(double x, double a, double b) noexcept {
  return (a == 0 && b == 1) ? x : (x - a) / (b - a);
}

// We use equidistant points for our interpolation problems. This is not
// optimal, but works sufficiently well.
inline std::vector<double> MakeInterpolationMask(
    int polynomial_degree) noexcept {
  std::vector<double> out(polynomial_degree);
  const double h = 1. / (polynomial_degree + 1);
  for (int i = 0; i < polynomial_degree; i++) {
    out[i] = ((i + 1) * h);
  }
  return out;
}

/**
 *  \ingroup Spline
 *  \brief Creates a T-vector of equidistant itnerpolation points for a given
 *         knot vector without any repetitions.
 **/
template <typename T>
inline std::vector<T> MakeInterpolationPoints(
    const std::vector<T> &uniq, const std::vector<T> &mask) noexcept {
  const int size = uniq.size();
  std::vector<T> out;

  out.reserve((size - 1) * mask.size());

  for (int i = 0; i < size - 1; i++) {
    for (auto m : mask) {
      out.push_back(uniq[i] + m * (uniq[i + 1] - uniq[i]));
    }
  }
  return out;
}

/**
 *  \ingroup Spline
 *  \brief returns the coefficients to represent a function in the Bernstein
 *         basis on [0,1].
 **/
inline Eigen::Matrix<double, -1, -1> GetInterpolationMatrix(
    int polynomial_degree, const std::vector<double> &mask) {
  Eigen::Matrix<double, -1, -1> interpolationMatrix(polynomial_degree + 1,
                                                    polynomial_degree + 1);

  double val[Constants::MaxP + 1];
  for (int j = 0; j < polynomial_degree + 1; j++) {
    Bembel::Basis::ShapeFunctionHandler::evalBasis(polynomial_degree, val,
                                                   mask[j]);
    for (int i = 0; i < polynomial_degree + 1; i++)
      interpolationMatrix(j, i) = val[i];
  }

  return interpolationMatrix.inverse();
}

/**
 *  \ingroup Spline
 *  \brief this solves a generic interpolation problem.
 **/
template <typename T>
inline void GetCoefficients(const int incr, const std::vector<double> &mask,
                            const std::vector<T> &values, T *coefs) {
  int polynomial_degree = mask.size() - 1;

  Eigen::Matrix<double, -1, -1> im =
      GetInterpolationMatrix(polynomial_degree, mask);

  Eigen::Matrix<T, -1, 1> rhs(polynomial_degree + 1);

  for (int i = 0; i < incr; i++) {
    for (int j = 0; j < polynomial_degree + 1; j++) {
      rhs[j] = values[i * (polynomial_degree + 1) + j];
    }
    Eigen::Matrix<T, -1, 1> tmp(polynomial_degree + 1);
    tmp = im * rhs;

    for (int j = 0; j < polynomial_degree + 1; j++) {
      coefs[i * (polynomial_degree + 1) + j] = tmp(j);
    }
  }

  return;
}

/**
 *  \ingroup Spline
 *  \brief this solves a generic interpolation problem. 
 **/
template <typename T>
inline std::vector<T> GetCoefficients(const int increments,
                                      const std::vector<double> &mask,
                                      const std::vector<T> &values) {
  std::vector<double> out(increments * (mask.size()));
  GetCoefficients<T>(increments, mask, values, out.data());
  return out;
}
}  // namespace Spl
}  // namespace Bembel
#endif
