// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef __BEMBEL_CONVERSIONS_
#define __BEMBEL_CONVERSIONS_

#include <Eigen/Dense>
#include "myvector.h"

namespace Bembel {

inline vector3 vector3_make(const Eigen::Vector3d &in) {
  vector3 c;
  c.x = in(0);
  c.y = in(1);
  c.z = in(2);
  return c;
}

inline Eigen::VectorXd ptr2eigen(double *ptr, int l) {
  Eigen::VectorXd out(l);
  for (int i = 0; i < l; i++) {
    out(i) = ptr[i];
  }
  return out;
}

inline Eigen::VectorXcd cmplxptr2eigen(double *ptr, int l) {
  const int sz = l / 2;
  Eigen::VectorXcd out(sz);
  for (int i = 0; i < sz; i++) {
    out(i) = std::complex<double>(ptr[i], ptr[i + sz]);
  }
  return out;
}

inline double *eigen2ptr(const Eigen::VectorXd &in) {
  double *out;
  const int sz = in.rows();
  out = (double *)malloc(sizeof(double) * sz);
  for (int i = 0; i < sz; i++) {
    out[i] = in(i);
  }
  return out;
}

inline double *eigen2cmplxptr(const Eigen::VectorXcd &in) {
  double *out;
  const int sz = in.rows();
  out = (double *)malloc(sizeof(double) * sz * 2);
  for (int i = 0; i < sz; i++) {
    out[i] = in(i).real();
  }
  for (int i = 0; i < sz; i++) {
    out[i + sz] = in(i).imag();
  }
  return out;
}

}  // namespace Bembel
#endif