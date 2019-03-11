// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef _localiseincluded_
#define _localiseincluded_

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <array>
#include "spline/basis.h"
#include "spline/helper.h"
/**
 *  @brief This header contains many routines useful for
 * localizing/interpolating spline spaces.
 */

namespace Spl {

constexpr inline double rescale(double x, double a, double b) {
  return (x - a) / (b - a);
}

inline std::vector<double> make_interpolation_mask(int deg) {
  std::vector<double> out(deg);
  const double h = 1. / (deg + 1);
  for (int i = 0; i < deg; i++) {
    out[i] = ((i + 1) * h);
  }
  return out;
}

// Creates a T-vector of equidistant itnerpolation points for a given knot
// vector without any repetitions.
template <typename T>
inline std::vector<T> make_interpolation_points(
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

// Returns the coefficients to represent a function in the Bernstein basis on
// [0,1].

inline Eigen::Matrix<double, -1, -1> getInterpolationMatrix(
    int deg, const std::vector<double> &mask) {
  Eigen::Matrix<double, -1, -1> interpolationMatrix(deg + 1, deg + 1);

  for (int j = 0; j < deg + 1; j++) {
    auto val = evalBrnstnBasis(deg, mask[j]);
    for (int i = 0; i < deg + 1; i++) interpolationMatrix(j, i) = val[i];
  }

  return interpolationMatrix.inverse();
}

template <typename T>
inline void get_coeffs(const int incr, const std::vector<double> &mask,
                       const std::vector<T> &vals, T *coefs) {
  int deg = mask.size() - 1;

  Eigen::Matrix<double, -1, -1> im = getInterpolationMatrix(deg, mask);

  Eigen::Matrix<T, -1, 1> rhs(deg + 1);

  for (int i = 0; i < incr; i++) {
    for (int j = 0; j < deg + 1; j++) {
      rhs[j] = vals[i * (deg + 1) + j];
    }
    Eigen::Matrix<T, -1, 1> tmp(deg + 1);
    tmp = im * rhs;

    for (int j = 0; j < deg + 1; j++) {
      coefs[i * (deg + 1) + j] = tmp(j);
    }
  }

  return;
}

template <typename T>
inline std::vector<T> get_coeffs(const int increments,
                                 const std::vector<double> &mask,
                                 const std::vector<T> &vals) {
  std::vector<double> out(increments * (mask.size()));
  get_coeffs<T>(increments, mask, vals, out.data());
  return out;
}
}  // namespace Spl
#endif
