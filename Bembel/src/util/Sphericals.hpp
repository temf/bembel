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

#ifndef BEMBEL_SRC_UTIL_SPHERICALS_HPP_
#define BEMBEL_SRC_UTIL_SPHERICALS_HPP_
#include <Eigen/Dense>

#if !defined pi
#define pi BEMBEL_PI
#endif

namespace Bembel {

double evaluate_sphericals(Eigen::Vector3d x, Eigen::VectorXd cs,
    unsigned int deg);

double evaluate_solid_sphericals(Eigen::Vector3d x, Eigen::VectorXd cs,
    unsigned int deg, bool grad);

Eigen::Vector3d evaluate_dsphericals(Eigen::Vector3d x, Eigen::VectorXd cs,
    unsigned int deg);

Eigen::Vector3d evaluate_dsolid_sphericals(Eigen::Vector3d x,
    Eigen::VectorXd cs, unsigned int deg);

Eigen::Matrix<double, Eigen::Dynamic, 2> spherical_harmonics_full(
    Eigen::Vector3d x, unsigned int N);

inline Eigen::Vector2d spherical_prev(Eigen::Vector3d x, int m, int n,
    Eigen::Vector2d y1, Eigen::Vector2d y2);

Eigen::Matrix<double, 3, Eigen::Dynamic> Dsolid_harmonics_full(
    Eigen::Vector3d x, unsigned int N, Eigen::MatrixXd spherical_val);

Eigen::Matrix<double, 3, 2> dsolid_spherical_prev(Eigen::Vector3d y, int m,
    unsigned int n, Eigen::VectorXd L, double y_re, double y_im);

Eigen::VectorXd legendreFull(unsigned int N, double t);

double constant(int m, unsigned int n);

inline Eigen::Matrix3d functionalMatrix(Eigen::Vector3d z);

inline double pow_int(double x, int n);

/**
 * \brief Evaluates the series \f$ \sum_{n = 0}^{\rm deg} \sum_{m = -n}^n c_m^n Y_m^n(x) \f$ for real coefficients,
 *        with the convenction that \f$ Y_{-m}^n := \overline{Y_m^n} \f$.
 *
 * \see K. Giebermann. _Schnelle Summationsverfahren zur numerischen Lösung von Integralgleichungen
 *  für Streuprobleme im_ \f$ \mathbb{R}^3 \f$. PhD thesis, Universität Karlsruhe (TH), Germany, 1997.
 * \see R. von Rickenbach. _Boundary Element Methods for Shape Optimisation in Homogenisation_.
 *  MSc thesis, Universität Basel, Switzerland, 2021.
 *
 * @param x:		The point of evaluation, a vector with length 1
 * @param cs:		The coefficients stored in the order \f$ [(0,  0), (1, -1), (1, 0), (1, 1)
 * 																			  (2, -2), (2, -1), ... (n, n)] \f$
 * @param deg:	The degree
 */
double evaluate_sphericals(Eigen::Vector3d x, Eigen::VectorXd cs,
    unsigned int deg) {
  unsigned int m, n;
  double z1[2], z2[2], z3[2], z1_start[2];
  double r, fac, rootTimesZ, root_2, root_3;

  assert(abs(x.norm() - 1) < Constants::generic_tolerance);

  r = z1[1] = 0;
  z1[0] = 0.5 / sqrt(pi);
  if (deg <= 1) {
    return cs(0) * z1[0];
  }

  for (m = 0; m < deg - 1; m++) {
    if (m == 0) {
      fac = 1.0;
    } else {
      fac = 2.0;
    }

    z1_start[0] = z1[0];
    z1_start[1] = z1[1];
    rootTimesZ = sqrt(2 * m + 3) * x(2);
    z2[0] = rootTimesZ * z1[0];
    z2[1] = rootTimesZ * z1[1];
    r += fac * cs(m * (m + 1) + m) * z1[0];  // + k[ m   *(m+1)-m]*z1[1];
    r += fac * cs((m + 1) * (m + 2) + m) * z2[0];  // + k[(m+1)*(m+2)-m]*z2[1];
    for (n = m + 2; n < deg; n++) {
      root_2 = sqrt((2 * n + 1.0) / ((n - m) * (n + m)));
      rootTimesZ = sqrt(2 * n - 1) * x(2);
      root_3 = sqrt(((n + m - 1.0) * (n - m - 1.0)) / (2 * n - 3));
      z3[0] = root_2 * (rootTimesZ * z2[0] - root_3 * z1[0]);
      z3[1] = root_2 * (rootTimesZ * z2[1] - root_3 * z1[1]);
      r += fac * cs(n * (n + 1) + m) * z3[0];  // + k[n*(n+1)-m]*z3[1];
      z1[0] = z2[0];
      z1[1] = z2[1];
      z2[0] = z3[0];
      z2[1] = z3[1];
    }
    root_2 = sqrt((2 * m + 3.0) / (2 * m + 2));
    z1[0] = root_2 * (x(0) * z1_start[0] - x(1) * z1_start[1]);
    z1[1] = root_2 * (x(0) * z1_start[1] + x(1) * z1_start[0]);
  }
  r += 2 * cs((deg - 1) * (deg + 1)) * z1[0];  // + k[(nk-1)*(nk-1)]*z1[1];
  return r;
}

/**
 * \brief Evaluates the series \f$ \sum_{n = 0}^{\rm deg} \sum_{m = -n}^n
 *   |x|^n c_m^n Y_m^n( \frac{x}{|x|}) \f$ for real coefficients cs if grad is false, and
 * 		\f$ \sum_{n = 0}^{\rm deg} \sum_{m = -n}^n n |x|^{n-1} c_m^n Y_m^n(\frac{x}{|x|}) \f$
 * 		for real coefficients cs if grad is true
 *
 * @param	x:		The point of evaluation
 * @param cs:		The coefficients stored in the order \f$ [(0,  0), (1, -1), (1, 0), (1, 1)
 * 																			  (2, -2), (2, -1), ... (n, n)] \f$
 * @param deg:	The degree
 * @param grad: distinguish the two cases.
 */
double evaluate_solid_sphericals(Eigen::Vector3d x, Eigen::VectorXd cs,
    unsigned int deg, bool grad) {
  unsigned int m, n;
  double z1[2], z2[2], z3[2], z1_start[2];
  double r, fac, rootTimesZ, root_2, root_3, norm, fac_tot;
  Eigen::Vector3d y;

  r = z1[1] = 0;
  z1[0] = 0.5 / sqrt(pi);
  if (grad && deg <= 1) {
    return 0.0;
  } else if (!grad && deg <= 1) {
    return cs(0) * z1[0];
  }

  norm = x.norm();
  y = x / norm;

  for (m = 0; m < deg - 1; m++) {
    if (m == 0) {
      fac = 1.0;
    } else {
      fac = 2.0;
    }

    if (grad) {
      fac_tot = fac * m * pow_int(norm, m - 1);
    } else {
      fac_tot = fac * pow_int(norm, m);
    }

    z1_start[0] = z1[0];
    z1_start[1] = z1[1];
    rootTimesZ = sqrt(2 * m + 3) * y(2);
    z2[0] = rootTimesZ * z1[0];
    z2[1] = rootTimesZ * z1[1];
    r += fac_tot * cs(m * (m + 1) + m) * z1[0];

    if (grad) {
      fac_tot = fac * (m + 1) * pow_int(norm, m);
    } else {
      fac_tot *= norm;
    }
    r += fac_tot * cs((m + 1) * (m + 2) + m) * z2[0];

    for (n = m + 2; n < deg; n++) {
      if (grad) {
        fac_tot = fac * n * pow_int(norm, n - 1);
      } else {
        fac_tot *= norm;
      }
      root_2 = sqrt((2 * n + 1.0) / ((n - m) * (n + m)));
      rootTimesZ = sqrt(2 * n - 1) * y(2);
      root_3 = sqrt(((n + m - 1.0) * (n - m - 1.0)) / (2 * n - 3));
      z3[0] = root_2 * (rootTimesZ * z2[0] - root_3 * z1[0]);
      z3[1] = root_2 * (rootTimesZ * z2[1] - root_3 * z1[1]);
      r += fac_tot * cs(n * (n + 1) + m) * z3[0];
      z1[0] = z2[0];
      z1[1] = z2[1];
      z2[0] = z3[0];
      z2[1] = z3[1];
    }
    root_2 = sqrt((2 * m + 3.0) / (2 * m + 2));
    z1[0] = root_2 * (y(0) * z1_start[0] - y(1) * z1_start[1]);
    z1[1] = root_2 * (y(0) * z1_start[1] + y(1) * z1_start[0]);
  }

  if (grad) {
    fac_tot = fac * deg * pow_int(norm, deg - 1);
  } else {
    fac_tot *= norm;
  }
  r += fac_tot * cs((deg - 1) * (deg + 1)) * z1[0];
  return r;
}

/**
 * \brief Evaluates the series \f$ \sum_{n = 0}^{\rm deg} \sum_{m = -n}^n c_m^n
 *  \nabla Y_m^n(x) \f$ for real coefficients.
 *
 * @param x:		The point of evaluation, a vector with length 1
 * @param cs:		The coefficients stored in the order \f$ [(0,  0), (1, -1), (1, 0), (1, 1)
 * 																			  (2, -2), (2, -1), ... (n, n)] \f$
 * @param deg:	The degree
 */
Eigen::Vector3d evaluate_dsphericals(Eigen::Vector3d x, Eigen::VectorXd cs,
    unsigned int deg) {
  unsigned int m, n;
  double z1[2], z2[2], z3[2], z1_start[2];
  Eigen::Vector3d dr;
  double fac, rootTimesZ, root_2, root_3;

  assert(abs(x.norm() - 1) < Constants::generic_tolerance);

  z1[0] = sqrt(0.375 / pi);
  z1[1] = 0;

  dr = Eigen::Vector3d(0.0, 0.0, 0.0);

  if (deg == 1) {
    return dr;
  }

  for (m = 1; m < deg - 1; m++) {
    z1_start[0] = z1[0];
    z1_start[1] = z1[1];
    rootTimesZ = sqrt(2 * m + 3) * x(2);
    z2[0] = rootTimesZ * z1[0];
    z2[1] = rootTimesZ * z1[1];
    dr(0) += 2 * m * (cs(m * (m + 1) + m) * z1[0]);
    // + a[ m   *(m+1)-m]*z1[1]) for imaginary coefficients;
    dr(0) += 2 * m * (cs((m + 1) * (m + 2) + m) * z2[0]);
    // + a[(m+1)*(m+2)-m]*z2[1]) for imaginary coefficients;
    dr(1) -= 2 * m * (cs(m * (m + 1) + m) * z1[1]);
    // - a[ m   *(m+1)-m]*z1[0]) for imaginary coefficients;
    dr(1) -= 2 * m * (cs((m + 1) * (m + 2) + m) * z2[1]);
    // - a[(m+1)*(m+2)-m]*z2[0]) for imaginary coefficients;

    if (m == 1) {
      fac = 1.0;
      dr(2) += fac * sqrt(2 * m)
          * (cs(m * (m + 1) + (m - 1)) * z1[0]
              + cs(m * (m + 1) - (m - 1)) * z1[1]);
      dr(2) += fac * sqrt(4 * m + 2)
          * (cs((m + 1) * (m + 2) + (m - 1)) * z2[0]
              + cs((m + 1) * (m + 2) - (m - 1)) * z2[1]);
    } else {
      fac = 2.0;
      dr(2) += fac * sqrt(2 * m) * (cs(m * (m + 1) + (m - 1)) * z1[0]);
      // + a[ m   *(m+1)-(m-1)]*z1[1]) for imaginary coefficients;
      dr(2) += fac * sqrt(4 * m + 2)
          * (cs((m + 1) * (m + 2) + (m - 1)) * z2[0]);
      // + a[(m+1)*(m+2)-(m-1)]*z2[1]) for imaginary coefficients;
    }

    for (n = m + 2; n < deg; n++) {
      root_2 = sqrt((2 * n + 1.0) / ((n - m) * (n + m)));
      rootTimesZ = sqrt(2 * n - 1) * x(2);
      root_3 = sqrt(((n + m - 1.0) * (n - m - 1.0)) / (2 * n - 3));
      z3[0] = root_2 * (rootTimesZ * z2[0] - root_3 * z1[0]);
      z3[1] = root_2 * (rootTimesZ * z2[1] - root_3 * z1[1]);
      dr(0) += 2 * m * (cs(n * (n + 1) + m) * z3[0]);  // + a[n*(n+1)-m]*z3[1]);
      dr(1) -= 2 * m * (cs(n * (n + 1) + m) * z3[1]);  // - a[n*(n+1)-m]*z3[0]);
      dr(2) += fac * sqrt((n + m) * (n - m + 1))
          * (cs(n * (n + 1) + (m - 1)) * z3[0]);  // + a[n*(n+1)-(m-1)]*z3[1]);
      z1[0] = z2[0];
      z1[1] = z2[1];
      z2[0] = z3[0];
      z2[1] = z3[1];
    }
    root_2 = sqrt((2 * m + 3.0) / (2 * m + 2));
    z1[0] = root_2 * (x(0) * z1_start[0] - x(1) * z1_start[1]);
    z1[1] = root_2 * (x(0) * z1_start[1] + x(1) * z1_start[0]);
  }
  dr(0) += (deg - 1) * (2 * cs((deg - 1) * (deg + 1)) * z1[0]);
  // + a[(na-1)*(na-1)]*z1[1]) for imaginary coefficients;
  dr(1) -= (deg - 1) * (2 * cs((deg - 1) * (deg + 1)) * z1[1]);
  // - a[(na-1)*(na-1)]*z1[0]) for imaginary coefficients;
  dr(2) += 2 * sqrt(2 * (deg - 1))
      * (cs(deg * (deg - 1) + (deg - 2)) * z1[0]
          + cs(deg * (deg - 1) - (deg - 2)) * z1[1]);

  return dr;
}

/**
 * \brief Evaluates the series \f$ \sum_{n = 0}^{\rm deg} \sum_{m = -n}^n
 * |x|^{n-3} c_m^n  (\nabla Y_m^n)(\frac{x}{|x|}) \f$ for real coefficients.
 *
 * @param x:		The point of evaluation
 * @param cs:		The coefficients stored in the order \f$ [(0,  0), (1, -1), (1, 0), (1, 1)
 * 																			  (2, -2), (2, -1), ... (n, n)] \f$
 * @param deg:	The degree
 */
Eigen::Vector3d evaluate_dsolid_sphericals(Eigen::Vector3d x,
    Eigen::VectorXd cs, unsigned int deg) {
  unsigned int m, n;
  double z1[2], z2[2], z3[2], z1_start[2];
  Eigen::Vector3d dr, y;
  double fac, rootTimesZ, root_2, root_3, norm, r_n3;

  z1[0] = sqrt(0.375 / pi);
  z1[1] = 0;

  dr = Eigen::Vector3d(0.0, 0.0, 0.0);
  norm = x.norm();
  y = x / norm;

  if (deg == 1) {
    return dr;
  }

  for (m = 1; m < deg - 1; m++) {
    r_n3 = pow_int(norm, m - 3);

    z1_start[0] = z1[0];
    z1_start[1] = z1[1];
    rootTimesZ = sqrt(2 * m + 3) * y(2);
    z2[0] = rootTimesZ * z1[0];
    z2[1] = rootTimesZ * z1[1];
    dr(0) += 2 * m * r_n3 * (cs(m * (m + 1) + m) * z1[0]);
    dr(0) += 2 * m * r_n3 * norm * (cs((m + 1) * (m + 2) + m) * z2[0]);
    dr(1) -= 2 * m * r_n3 * (cs(m * (m + 1) + m) * z1[1]);
    dr(1) -= 2 * m * r_n3 * norm * (cs((m + 1) * (m + 2) + m) * z2[1]);

    if (m == 1) {
      fac = 1.0;
      dr(2) += fac * r_n3 * sqrt(2 * m)
          * (cs(m * (m + 1) + (m - 1)) * z1[0]
              + cs(m * (m + 1) - (m - 1)) * z1[1]);
      dr(2) += fac * r_n3 * norm * sqrt(4 * m + 2)
          * (cs((m + 1) * (m + 2) + (m - 1)) * z2[0]
              + cs((m + 1) * (m + 2) - (m - 1)) * z2[1]);
    } else {
      fac = 2.0;
      dr(2) += fac * r_n3 * sqrt(2 * m) * (cs(m * (m + 1) + (m - 1)) * z1[0]);
      dr(2) += fac * r_n3 * norm * sqrt(4 * m + 2)
          * (cs((m + 1) * (m + 2) + (m - 1)) * z2[0]);
    }

    r_n3 *= norm;

    for (n = m + 2; n < deg; n++) {
      r_n3 *= norm;
      root_2 = sqrt((2 * n + 1.0) / ((n - m) * (n + m)));
      rootTimesZ = sqrt(2 * n - 1) * y(2);
      root_3 = sqrt(((n + m - 1.0) * (n - m - 1.0)) / (2 * n - 3));
      z3[0] = root_2 * (rootTimesZ * z2[0] - root_3 * z1[0]);
      z3[1] = root_2 * (rootTimesZ * z2[1] - root_3 * z1[1]);
      dr(0) += 2 * m * r_n3 * (cs(n * (n + 1) + m) * z3[0]);
      dr(1) -= 2 * m * r_n3 * (cs(n * (n + 1) + m) * z3[1]);
      dr(2) += fac * r_n3 * sqrt((n + m) * (n - m + 1))
          * (cs(n * (n + 1) + (m - 1)) * z3[0]);
      z1[0] = z2[0];
      z1[1] = z2[1];
      z2[0] = z3[0];
      z2[1] = z3[1];
    }
    root_2 = sqrt((2 * m + 3.0) / (2 * m + 2));
    z1[0] = root_2 * (y(0) * z1_start[0] - y(1) * z1_start[1]);
    z1[1] = root_2 * (y(0) * z1_start[1] + y(1) * z1_start[0]);
  }
  r_n3 *= norm;
  dr(0) += (deg - 1) * r_n3 * (2 * cs((deg - 1) * (deg + 1)) * z1[0]);
  dr(1) -= (deg - 1) * r_n3 * (2 * cs((deg - 1) * (deg + 1)) * z1[1]);
  dr(2) += 2 * r_n3 * sqrt(2 * (deg - 1))
      * (cs(deg * (deg - 1) + (deg - 2)) * z1[0]
          + cs(deg * (deg - 1) - (deg - 2)) * z1[1]);

  return dr;
}

/**
 * \brief Calculates the the spherical harmonics \f$ Y_n^m(\frac{x}{|x|}) \f$,
 *  ordered by \f$ [Y_0^0, \, Y_1^{-1}, \, Y_1^0, \, Y_1^1, \,
 *  Y_2^{-2}, ..., Y_N^N] \f$
 *
 *  @param  x:    the point of evaluation
 *  @param  N:    the maximal degree
 */
Eigen::Matrix<double, Eigen::Dynamic, 2> spherical_harmonics_full(
    Eigen::Vector3d x, unsigned int N) {

  Eigen::VectorXd real((N + 1) * (N + 1));
  Eigen::VectorXd imag((N + 1) * (N + 1));

  Eigen::Vector3d y = x / x.norm();

  int m, n;

  /* handle n = 0 separately */
  real(0) = 0.5 / sqrt(pi);
  imag(0) = 0.0;

  /* assemble two temporary vectors */
  Eigen::Vector2d z1, z2, tmp;

  for (n = 1; n <= N; n++) {
    for (m = 0; m < n - 1; m++) {
      z1(0) = real((n - 1) * (n - 1) + n - 1 + m);
      z1(1) = real((n - 2) * (n - 2) + n - 2 + m);
      z2(0) = imag((n - 1) * (n - 1) + n - 1 + m);
      z2(1) = imag((n - 2) * (n - 2) + n - 2 + m);

      tmp = spherical_prev(y, m, n, z1, z2);

      real(n * n + n + m) = tmp(0);
      imag(n * n + n + m) = tmp(1);

      // if m > 0, copy the value with the prefactor
      if (m > 0) {
        real(n * n + n - m) = tmp(0);
        imag(n * n + n - m) = -tmp(1);
      }
    }

    /* for m = n-1 */
    m = n - 1;
    z1(0) = real(m * m + m + m);
    z2(0) = imag(m * m + m + m);

    tmp = spherical_prev(y, m, n, z1, z2);
    real(n * n + n + m) = tmp(0);
    imag(n * n + n + m) = tmp(1);

    if (m > 0) {
      real(n * n + n - m) = tmp(0);
      imag(n * n + n - m) = -tmp(1);
    }

    /* for m = n */
    m = n;
    z1(0) = real((n - 1) * (n - 1) + n - 1 + n - 1);
    z2(0) = imag((n - 1) * (n - 1) + n - 1 + n - 1);

    tmp = spherical_prev(y, m, n, z1, z2);
    real(n * n + n + m) = tmp(0);
    imag(n * n + n + m) = tmp(1);

    if (m > 0) {
      real(n * n + n - m) = tmp(0);
      imag(n * n + n - m) = -tmp(1);
    }
  }

  Eigen::MatrixXd res((N + 1) * (N + 1), 2);
  res.col(0) = real;
  res.col(1) = imag;

  return res;
}

/**
 * \brief Calculates the spherical harmonic \f$ Y_n^m(x) \f$ based on previous values
 */
inline Eigen::Vector2d spherical_prev(Eigen::Vector3d x, int m, int n,
    Eigen::Vector2d y1, Eigen::Vector2d y2) {
  Eigen::Vector2d z;

  assert(abs(x.norm() - 1) < 1e-14);

  if ((m == 0) && (n == 0)) {
    z(0) = 0.5 / sqrt(pi);
    z(1) = 0;
  } else if (m == n) {
    z(0) = sqrt((2 * m + 1.0) / (2 * m)) * (x(0) * y1(0) - x(1) * y2(0));
    z(1) = sqrt((2 * m + 1.0) / (2 * m)) * (x(0) * y2(0) + x(1) * y1(0));
  } else if (m + 1 == n) {
    z(0) = y1(0) * sqrt(2 * m + 3) * x(2);
    z(1) = y2(0) * sqrt(2 * m + 3) * x(2);
  } else {
    z(0) = sqrt((2 * n + 1.0) / ((n - m) * (n + m)))
        * (sqrt(2 * n - 1) * x(2) * y1(0)
            - sqrt(((n + m - 1.0) * (n - m - 1.0)) / (2 * n - 3)) * y1(1));
    z(1) = sqrt((2 * n + 1.0) / ((n - m) * (n + m)))
        * (sqrt(2 * n - 1) * x(2) * y2(0)
            - sqrt(((n + m - 1.0) * (n - m - 1.0)) / (2 * n - 3)) * y2(1));
  }
  return z;
}

/**
 * \brief Calculates all gradients of the solid harmonics, given the Legendre Coefficients L.
 */
Eigen::Matrix<double, 3, Eigen::Dynamic> Dsolid_harmonics_full(
    Eigen::Vector3d x, unsigned int N, Eigen::MatrixXd spherical_val) {

  Eigen::VectorXd L = legendreFull(N, x(2) / x.norm());
  Eigen::Matrix<double, 3, 2> z;

  Eigen::MatrixXd reals(3, (N + 1) * (N + 1)), imags(3, (N + 1) * (N + 1));

  int m, n;

  for (n = 0; n <= N; n++) {
    for (m = 0; m <= n; m++) {
      z = dsolid_spherical_prev(x, m, n, L, spherical_val(n * n + n + m, 0),
          spherical_val(n * n + n + m, 1));

      reals(0, n * n + n + m) = z(0, 0);
      reals(1, n * n + n + m) = z(1, 0);
      reals(2, n * n + n + m) = z(2, 0);

      imags(0, n * n + n + m) = z(0, 1);
      imags(1, n * n + n + m) = z(1, 1);
      imags(2, n * n + n + m) = z(2, 1);

      if (m > 0) {
        reals(0, n * n + n - m) = z(0, 0);
        reals(1, n * n + n - m) = z(1, 0);
        reals(2, n * n + n - m) = z(2, 0);

        imags(0, n * n + n - m) = -z(0, 1);
        imags(1, n * n + n - m) = -z(1, 1);
        imags(2, n * n + n - m) = -z(2, 1);
      }
    }
  }
  return reals;
}

/**
 * \brief Calculates \f$ \nabla Y_m^n(x) \f$ based on the Legendre Coefficients L.
 */
inline Eigen::Matrix<double, 3, 2> dspherical_prev(Eigen::Vector3d x, int m,
    unsigned int n, Eigen::VectorXd L) {

  double c;
  unsigned int i;

  assert(abs(x.norm() - 1) < 1e-14);

  Eigen::Matrix<double, 3, 2> z;

  if (m == 0) {
    z(0, 0) = z(0, 1) = z(2, 1) = 0;
    if (1 > n) {
      z(2, 0) = 0;
    } else {
      z(2, 0) = L((n * (n + 1)) / 2 + 1);
    }
  } else {
    if (m > n) {
      z(0, 0) = 0;
    } else {
      z(0, 0) = L((n * (n + 1)) / 2 + m);
    }

    if (m + 1 > n) {
      z(2, 0) = 0;
    } else {
      z(2, 0) = L((n * (n + 1)) / 2 + m + 1);
    }
    z(0, 1) = z(2, 1) = 0;

    c = x(0) * z(2, 0) - x(1) * z(2, 1);
    z(2, 1) = x(0) * z(2, 1) + x(1) * z(2, 0);
    z(2, 0) = c;

    for (i = 1; i < m; i++) {
      c = x(0) * z(0, 0) - x(1) * z(0, 1);
      z(0, 1) = x(0) * z(0, 1) + x(1) * z(0, 0);
      z(0, 0) = c;

      c = x(0) * z(2, 0) - x(1) * z(2, 1);
      z(2, 1) = x(0) * z(2, 1) + x(1) * z(2, 0);
      z(2, 0) = c;
    }
  }

  c = constant(m, n);
  z(0, 0) *= m * c;
  z(0, 1) *= m * c;
  z(1, 0) = -z(0, 1);
  z(1, 1) = +z(0, 0);
  z(2, 0) *= c;
  z(2, 1) *= c;

  return z;
}

/**
 * \brief Calculates the derivative of the solid harmonics function \f$ \phi_n^m \f$ at
 *  the point \f$ y \f$ with the spherical values.
 */
Eigen::Matrix<double, 3, 2> dsolid_spherical_prev(Eigen::Vector3d y, int m,
    unsigned int n, Eigen::VectorXd L, double y_re, double y_im) {

  //  normalise y
  // x denotes the normalised y and is passed to the original functions
  double r = y.norm();
  Eigen::Vector3d x = y / r;
  Eigen::Matrix<double, 3, 2> z_harm;
  Eigen::Matrix3d A;
  unsigned int m_tilde;

  // part with the gradient of Y
  if (m >= 0) {
    m_tilde = m;
    z_harm = dspherical_prev(x, m_tilde, n, L);

    // multiply the gradient of Y by r^(n-3) and the matrix A
    double r_n3;
    if (n > 3) {
      r_n3 = pow(r, n - 3);
    } else if (n == 3) {
      r_n3 = 1.0;
    } else {
      r_n3 = pow(1.0 / r, 3 - n);
    }

    A = functionalMatrix(y);
    z_harm = r_n3 * (A * z_harm);

  } else {
    m_tilde = -m;
    z_harm = dspherical_prev(x, m_tilde, n, L);

    // define r^(n-3) as above
    double r_n3;
    if (n > 3) {
      r_n3 = pow(r, n - 3);
    } else if (n == 3) {
      r_n3 = 1.0;
    } else {
      r_n3 = pow(1 / r, 3 - n);
    }

    // get the part of the gradient of Y,
    // conjugate and multiply by the functional matrix
    A = functionalMatrix(y);
    z_harm = r_n3 * (A * z_harm);

    // get the complex conjugation
    z_harm.col(1) *= -1.0;
  }

  // calculate the radial part and multiply with n*r^{n-1}
  Eigen::Matrix<double, 3, 2> z_rad;

  // case n = 0: r^n = 1 is constant
  if (n == 0) {
    z_rad.setZero();
  } else {
    // the gradient of r^n is n*r^(n-2)*y = nr^(n-1)*x
    // again adapt the power
    double nr_n1;
    if (n == 1) {
      nr_n1 = 1.0;
    } else {
      nr_n1 = n * pow(r, n - 1);
    }

    double y_imt;
    if (m >= 0) {
      y_imt = y_im;
    } else {
      y_imt = -y_im;
    }

    z_rad.col(0) = y_re * nr_n1 * x;
    z_rad.col(1) = y_imt * nr_n1 * x;
  }
  return z_harm + z_rad;
}

/**
 * \brief Returns the values of the spherical polynomials \f$ P_n^m(t) \f$
 */
Eigen::VectorXd legendreFull(unsigned int N, double t) {
  Eigen::VectorXd L(((N + 1) * (N + 2)) / 2);
  int n, s, m;

  /* n = 0 */
  L(0) = 1;

  /* n = 1 */
  L(1) = t;
  L(2) = 1;

  for (n = 2; n <= N; n++) {
    s = (n * (n + 1)) / 2;
    for (m = 0; m < n - 1; m++) {
      L(s + m) = ((2 * n - 1) * t * L((n * (n - 1)) / 2 + m)
          - (n + m - 1) * L(((n - 1) * (n - 2)) / 2 + m)) / (n - m);
    }

    /* m = n-1 */
    m = n - 1;
    L(s + m) = (2 * m + 1) * t * L((m * (m + 1)) / 2 + m);

    /* m = n */
    m = n;
    L(s + m) = (2 * m - 1) * L((n * (n - 1)) / 2 + n - 1);
  }
  return L;
}

/**
 *  \brief gives the norming factor of the spherical polynomials
 */
double constant(int m, unsigned int n) {
  unsigned int i;
  double c = sqrt((2 * n + 1) / (4 * pi));

  for (i = 0; i < m; i++) {
    c *= 1 / sqrt((n + i + 1) * (n - i));
  }

  return c;
}

/**
 * \brief gives \f$ |z|^3 \f$ times the Jacobi Matrix of the transformation \f$ z \mapsto \frac{z}{|z|} \f$
 */
inline Eigen::Matrix3d functionalMatrix(Eigen::Vector3d z) {
  Eigen::Matrix3d M;

  double r = z.norm();

  M(0, 0) = r * r - z(0) * z(0);
  M(1, 1) = r * r - z(1) * z(1);
  M(2, 2) = r * r - z(2) * z(2);

  M(1, 0) = M(0, 1) = -z(0) * z(1);
  M(2, 0) = M(0, 2) = -z(0) * z(2);
  M(2, 1) = M(1, 2) = -z(1) * z(2);

  return M;
}

/**
 * \brief returns the \f$ n \f$-th power of \f$ x \f$ without using pow.
 */
inline double pow_int(double x, int n) {
  if (n == 0) {
    return 1.0;
  }

  unsigned int ul, k;
  if (n > 0) {
    ul = n;
  } else {
    ul = -n;
  }

  double res = x;
  for (k = 1; k < ul; k++) {
    res *= x;
  }

  if (n > 0) {
    return res;
  } else {
    return 1.0 / res;
  }
}

}  // namespace Bembel

#endif  // BEMBEL_SRC_UTIL_SPHERICALS_HPP_
