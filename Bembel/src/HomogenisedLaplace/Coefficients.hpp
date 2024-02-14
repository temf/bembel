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

#ifndef BEMBEL_SRC_HOMOGENISEDLAPLACE_COEFFICIENTS_HPP_
#define BEMBEL_SRC_HOMOGENISEDLAPLACE_COEFFICIENTS_HPP_

#ifndef POINT_DEGREE
#define POINT_DEGREE 20 /** the number of points on the surface of the cube */
#endif

#include <Eigen/Dense>
#include <Bembel/HomogenisedLaplace>
#include <Bembel/Quadrature>

#include <functional>

namespace Bembel {

Eigen::VectorXd getCoefficients(double precision);

inline unsigned int getDegree(double precision);

Eigen::VectorXd getDisplacement(Eigen::MatrixXd ps_l, Eigen::MatrixXd ps_f,
    Eigen::MatrixXd ps_b);

inline double k_mod(Eigen::Vector3d in);

inline Eigen::Vector3d Dk_mod(Eigen::Vector3d in);

double calculateFirstCoefficient(Eigen::VectorXd cs, unsigned int deg,
    Eigen::MatrixXd ps_l, Eigen::MatrixXd ps_f, Eigen::MatrixXd ps_b);

/**
 * \brief Calculates the coefficients for the solid harmonics expansion of
 *  the periodic kernel \f$ k(x) = \sum_{m \in \{-1, 0, 1\}^3} \frac{1}{4 \pi |x-m|}
 *  + \frac{|x|^2}{6} + \sum_{n = 0}^{N} \sum_{m = -n}^n c_m^n \phi_n^m(x) \f$.
 *
 * @param precision:  The desired mean error from periodicity.
 *
 * \see A. Barnett and L. Greengard. A new integral representation for quasi-periodic
 *  fields and its application to two-dimensional band structure calculations.
 *  _Journal of Computational Physics_, 229(19):6898--6914, 2010.
 *
 * \see P. Cazeaux and O. Zahm. A fast boundary element method
 *  for the solution of periodic many-inclusion problems via
 *  hierarchical matrix techniques. _ESAIM: Proceedings and Surveys_,
 *  48:156--168, 2015.
 *
 * \see R. von Rickenbach. _Boundary Element Methods for Shape
 *  Optimisation in Homogenisation_. MSc thesis, Universität Basel, 2021.
 */
Eigen::VectorXd getCoefficients(double precision) {
  Eigen::VectorXd diff;
  Eigen::RowVectorXd difft;
  Eigen::Vector3d v;
  Eigen::MatrixXd spherical_values_pre, spherical_values_pst, Dsolid_values_pre,
      Dsolid_values_pst;

  unsigned int deg = getDegree(precision);
  unsigned int Msquare = (POINT_DEGREE + 1) * (POINT_DEGREE + 1);

  unsigned int m, n, k;
  double scale, fac, norm;

  GaussSquare<POINT_DEGREE> GS;
  Eigen::MatrixXd xs = GS[POINT_DEGREE].xi_;
  xs -= 0.5 * Eigen::MatrixXd::Ones(xs.rows(), xs.cols());

  Eigen::Vector3d ex(1.0, 0.0, 0.0);
  Eigen::Vector3d ey(0.0, 1.0, 0.0);
  Eigen::Vector3d ez(0.0, 0.0, 1.0);

  Eigen::MatrixXd ps_left(3, xs.cols());
  ps_left.row(0) = -0.5 * Eigen::VectorXd::Ones(xs.cols());
  ps_left.block(1, 0, 2, xs.cols()) = xs.block(0, 0, 2, xs.cols());

  Eigen::MatrixXd ps_front(3, xs.cols());
  ps_front.row(0) = xs.row(0);
  ps_front.row(1) = -0.5 * Eigen::VectorXd::Ones(xs.cols());
  ps_front.row(2) = xs.row(1);

  Eigen::MatrixXd ps_bottom(3, xs.cols());
  ps_bottom.block(0, 0, 2, xs.cols()) = xs.block(0, 0, 2, xs.cols());
  ps_bottom.row(2) = -0.5 * Eigen::VectorXd::Ones(xs.cols());

  Eigen::VectorXd displacement = getDisplacement(ps_left, ps_front, ps_bottom);

  Eigen::MatrixXd systemMatrix(6 * Msquare, ((deg + 1) * (deg + 2)) / 2 - 1);

  for (k = 0; k < Msquare; k++) {
    /* left - right difference */
    v = ps_left.col(k);
    norm = v.norm();

    spherical_values_pre = spherical_harmonics_full(v, deg);
    spherical_values_pst = spherical_harmonics_full(v + ex, deg);

    Dsolid_values_pre = Dsolid_harmonics_full(v, deg, spherical_values_pre);
    Dsolid_values_pst = Dsolid_harmonics_full(v + ex, deg,
        spherical_values_pst);

    for (n = 1; n <= deg; n++) {
      scale = pow(norm, n);
      diff = scale
          * (spherical_values_pre.block(n * n + n, 0, n + 1, 1)
              - spherical_values_pst.block(n * n + n, 0, n + 1, 1));
      difft = Dsolid_values_pre.block(0, n * n + n, 1, n + 1)
          - Dsolid_values_pst.block(0, n * n + n, 1, n + 1);

      /* adapt the scaling */
      diff.segment(1, n) *= 2.0;
      difft.segment(1, n) *= 2.0;

      systemMatrix.block(k, (n * (n + 1)) / 2 - 1, 1, n + 1) = diff.transpose();
      systemMatrix.block(k + Msquare, (n * (n + 1)) / 2 - 1, 1, n + 1) = difft;
    }

    /* front - back difference */
    v = ps_front.col(k);
    norm = v.norm();

    spherical_values_pre = spherical_harmonics_full(v, deg);
    spherical_values_pst = spherical_harmonics_full(v + ey, deg);

    Dsolid_values_pre = Dsolid_harmonics_full(v, deg, spherical_values_pre);
    Dsolid_values_pst = Dsolid_harmonics_full(v + ey, deg,
        spherical_values_pst);

    for (n = 1; n <= deg; n++) {
      scale = pow(norm, n);
      diff = scale
          * (spherical_values_pre.block(n * n + n, 0, n + 1, 1)
              - spherical_values_pst.block(n * n + n, 0, n + 1, 1));
      difft = Dsolid_values_pre.block(1, n * n + n, 1, n + 1)
          - Dsolid_values_pst.block(1, n * n + n, 1, n + 1);

      /* adapt the scaling */
      diff.segment(1, n) *= 2.0;
      difft.segment(1, n) *= 2.0;

      systemMatrix.block(k + 2 * Msquare, (n * (n + 1)) / 2 - 1, 1, n + 1) =
          diff.transpose();
      systemMatrix.block(k + 3 * Msquare, (n * (n + 1)) / 2 - 1, 1, n + 1) =
          difft;
    }

    /* bottom - top difference */
    v = ps_bottom.col(k);
    norm = v.norm();

    spherical_values_pre = spherical_harmonics_full(v, deg);
    spherical_values_pst = spherical_harmonics_full(v + ez, deg);

    Dsolid_values_pre = Dsolid_harmonics_full(v, deg, spherical_values_pre);
    Dsolid_values_pst = Dsolid_harmonics_full(v + ez, deg,
        spherical_values_pst);

    for (n = 1; n <= deg; n++) {
      scale = pow(norm, n);
      diff = scale
          * (spherical_values_pre.block(n * n + n, 0, n + 1, 1)
              - spherical_values_pst.block(n * n + n, 0, n + 1, 1));
      difft = Dsolid_values_pre.block(2, n * n + n, 1, n + 1)
          - Dsolid_values_pst.block(2, n * n + n, 1, n + 1);

      /* adapt the scaling */
      diff.segment(1, n) *= 2.0;
      difft.segment(1, n) *= 2.0;

      systemMatrix.block(k + 4 * Msquare, (n * (n + 1)) / 2 - 1, 1, n + 1) =
          diff.transpose();
      systemMatrix.block(k + 5 * Msquare, (n * (n + 1)) / 2 - 1, 1, n + 1) =
          difft;
    }
  }

  /* solve the system */
  Eigen::VectorXd coeffs(((deg + 1) * (deg + 2)) / 2);
  coeffs.setZero();
  coeffs.segment(1, ((deg + 1) * (deg + 2)) / 2 - 1) =
      systemMatrix.colPivHouseholderQr().solve(-displacement);

  /* Copy the stuff into the full Coefficient list */
  Eigen::VectorXd coeffs_full((deg + 1) * (deg + 1));
  coeffs_full(0) = 0;
  for (n = 1; n <= deg; n++) {
    coeffs_full(n * n + n) = coeffs((n * (n + 1)) / 2);
    for (m = 1; m <= n; m++) {
      coeffs_full(n * n + n + m) = coeffs((n * (n + 1)) / 2 + m);
      coeffs_full(n * n + n - m) = coeffs((n * (n + 1)) / 2 + m);
    }
  }

  /* calculate the first coefficient */
  coeffs_full(0) = calculateFirstCoefficient(coeffs_full, deg, ps_left,
      ps_front, ps_bottom);

  return coeffs_full;
}

/**
 * \brief Returns the degree of the sphericals expansion
 * given a precision. Can be extended, use even numbers only!
 *
 * \see P. Cazeaux and O. Zahm. A fast boundary element method
 *  for the solution of periodic many-inclusion problems via
 *  hierarchical matrix techniques. _ESAIM: Proceedings and Surveys_,
 *  48:156--168, 2015.
 *
 * \see R. von Rickenbach. _Boundary Element Methods for Shape
 *  Optimisation in Homogenisation_. MSc thesis, Universität Basel, 2021.
 */
inline unsigned int getDegree(double precision) {
  if (precision > 1e-4) {
    return 4;
  } else if (precision > 1e-6) {
    return 8;
  } else {
    return 12;
  }
}

/**
 * \brief Returns the right-hand side for the homogenised coefficient calculation.
 *
 * @param ps_l: Points on the left side, i.e. \f$ \subseteq \{-0.5\} \times [-0.5, 0.5] \times [-0.5, 0.5] \f$.
 * @param ps_f: Points on the front side, i.e. \f$ \subseteq [-0.5, 0.5] \times \{-0.5\} \times [-0.5, 0.5] \f$.
 * @param ps_b: Points on the bottom side, i.e. \f$ \subseteq [-0.5, 0.5] \times [-0.5, 0.5] \times \{-0.5\} \f$.
 */
Eigen::VectorXd getDisplacement(Eigen::MatrixXd ps_l, Eigen::MatrixXd ps_f,
    Eigen::MatrixXd ps_b) {

  std::function<double(Eigen::Vector3d)> u = [](Eigen::Vector3d in) {
    return k_mod(in);
  };
  std::function<Eigen::Vector3d(Eigen::Vector3d)> Du = [](Eigen::Vector3d in) {
    return Dk_mod(in);
  };

  unsigned int Msquare = (POINT_DEGREE + 1) * (POINT_DEGREE + 1);
  unsigned int k;

  Eigen::Vector3d ex(1.0, 0.0, 0.0);
  Eigen::Vector3d ey(0.0, 1.0, 0.0);
  Eigen::Vector3d ez(0.0, 0.0, 1.0);

  Eigen::Vector3d tmp;
  Eigen::VectorXd d(6 * Msquare); /* the return vector */

  /* left - right displacement */
  for (k = 0; k < Msquare; k++) {
    tmp = Du(ps_l.col(k)) - Du(ps_l.col(k) + ex);
    d(k) = u(ps_l.col(k)) - u(ps_l.col(k) + ex);
    d(k + Msquare) = tmp(0);
  }

  /* front - back displacement */
  for (k = 0; k < Msquare; k++) {
    tmp = Du(ps_f.col(k)) - Du(ps_f.col(k) + ey);
    d(k + 2 * Msquare) = u(ps_f.col(k)) - u(ps_f.col(k) + ey);
    d(k + 3 * Msquare) = tmp(1);
  }

  /* bottom - top displacement */
  for (k = 0; k < Msquare; k++) {
    tmp = Du(ps_b.col(k)) - Du(ps_b.col(k) + ez);
    d(k + 4 * Msquare) = u(ps_b.col(k)) - u(ps_b.col(k) + ez);
    d(k + 5 * Msquare) = tmp(2);
  }

  return d;
}

/**
 * \brief Returns the modified kernel \f$ \sum_{m \in \{-1, 0, 1 \}^3} \frac{1}{4 \pi |x-m|}
 *  + \frac{|x|^2}{6} \f$.
 */
inline double k_mod(Eigen::Vector3d in) {
  double r = 0.0;
  int i, j, k;
  Eigen::Vector3d m;

  for (i = -1; i <= 1; i++) {
    for (j = -1; j <= 1; j++) {
      for (k = -1; k <= 1; k++) {
        m = Eigen::Vector3d(i, j, k);
        r += 1.0 / ((in - m).norm());
      }
    }
  }

  r /= (4 * M_PI);

  /* the part to ensure the vanishing mean on the Laplacian */
  r += (in.dot(in)) / 6.0;

  return r;
}

/**
 * \brief Returns the gradient of the modified kernel \f$ - \sum_{m \in \{-1, 0, 1\}^3}
 *  \frac{x-m}{|x-m|^3} + \frac{x}{3} \f$.
 */
inline Eigen::Vector3d Dk_mod(Eigen::Vector3d in) {
  Eigen::Vector3d r, s;
  double snorm;
  r.setZero();

  int i, j, k;
  for (i = -1; i <= 1; i++) {
    for (j = -1; j <= 1; j++) {
      for (k = -1; k <= 1; k++) {
        s = in - Eigen::Vector3d(i, j, k);
        snorm = s.norm();
        r -= s / (snorm * snorm * snorm);
      }
    }
  }

  r /= (4.0 * M_PI);

  /* the part to ensure the vanishing mean on the Laplacian */
  r += in / 3.0;

  return r;
}

/**
 * \brief Calculates the first coefficient such that the mean of the
 *  kernel vanishes.
 *
 * @param cs:   The already calculated other coefficients
 * @param deg:  The degree
 * @param ps_l: Gauss quadrature points on the left side of the cube
 * @param ps_f: Gauss quadrature points on the front side of the cube
 * @param ps_b: Gauss quadrature points on the bottom side of the cube
 */
double calculateFirstCoefficient(Eigen::VectorXd cs, unsigned int deg,
    Eigen::MatrixXd ps_l, Eigen::MatrixXd ps_f, Eigen::MatrixXd ps_b) {
  double res = 0.0;

  GaussSquare<POINT_DEGREE> GS;
  Eigen::VectorXd ws = GS[POINT_DEGREE].w_;

  Eigen::VectorXd cs_tmp(cs.rows());
  cs_tmp.setZero();

  unsigned int k, n;
  double norm;

  for (k = 0; k < ps_l.cols(); k++) {
    norm = ps_l.col(k).norm(); /* equal for ps_f and ps_b */
    for (n = 1; n <= deg; n++) {
      cs_tmp.segment(n * n, 2 * n + 1) = pow(norm, n)
          * cs.segment(n * n, 2 * n + 1);
    }

    res += ws(k)
        * (k_mod(ps_l.col(k))
            + evaluate_sphericals(ps_l.col(k) / norm, cs_tmp, deg));
    res += ws(k)
        * (k_mod(ps_f.col(k))
            + evaluate_sphericals(ps_f.col(k) / norm, cs_tmp, deg));
    res += ws(k)
        * (k_mod(ps_b.col(k))
            + evaluate_sphericals(ps_b.col(k) / norm, cs_tmp, deg));
  }

  res /= 3.0;

  res += (1.0 / 24);

  return -2.0 * sqrt(M_PI) * res;
}

} /* namespace Bembel */

#endif  // BEMBEL_SRC_HOMOGENISEDLAPLACE_COEFFICIENTS_HPP_
