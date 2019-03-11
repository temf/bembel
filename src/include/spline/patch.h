// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef _PATCH_INCLUDED_
#define _PATCH_INCLUDED_

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
#include "spline/basis.h"
#include "spline/bezier_extraction.h"
#include "spline/knots.h"

#ifdef _BEMBEL_OLD_
#include "myvector.h"
#endif

namespace Spl {

// The patch class.
// data stores the controll points in Superspace format, i.e., Bezier
// extraction.
// Attention! turnNormal and flip : are for reorientation of the
// patch without reparametrisation. They only affect the basel-wrapper
//
// About the organisation of this file:
// First are the class member vars. Then, there are two larger methods for
// evalutaion of points and the jacobian. Later one, those are are called by
// either Eigen- or basel-wrapper

class Patch {
 public:
  std::vector<double> data;  // Controllpoints in Bezier-Extracted Format.
  int xp;                    // Degree in x
  int yp;                    // Degree in y
  std::vector<double>
      xuniqueknt;  // The knot vectors, where each knot is unique
  std::vector<double>
      yuniqueknt;  // The knot vectors, where each knot is unique
  inline void init(const std::vector<Eigen::Matrix<double, -1, -1>> &mats,
                   const std::vector<double> &knotx,
                   const std::vector<double> &knoty) {
    initHom(mats, knotx, knoty);
  }
  /* initHom() is just a Bezier Extraction via superspace embedding. Data gets
   * stored in the std::vector data. */
  void initHom(const std::vector<Eigen::Matrix<double, -1, -1>> &mats,
               const std::vector<double> &knotx,
               const std::vector<double> &knoty);

  /* eval() evaluates the geometry. I look up the position in the knot vector,
   * scale the input arguments, evaluate the 1D basis functions and sum over
   * them with the controll points from data. */

  Eigen::Vector3d evalCore(double x, double y) const;

  Eigen::Matrix<double, 3, 2> jacobianCore(const double x,
                                           const double y) const;

  inline Eigen::Matrix<double, 3, 1> eval(const double x,
                                          const double y) const {
    return evalCore(x, y);
  }

  // Eigen::Matrix<double, 3, 1> eval(const double x, const double y) const {
  //   std::array<double, 3> tmp = evalCore(x, y);
  //   Eigen::Matrix<double, 3, 1> out;
  //   for (int k = 0; k < 3; k++)
  //     out(k) = tmp[k];
  //   return out;
  // }

  inline Eigen::Matrix<double, 3, 2> jacobian(const double x,
                                              const double y) const {
    return jacobianCore(x, y);
  }

  // Eigen::Matrix<double, 3, 2> jacobian(const double x, const double y)
  // const
  // {
  //   Eigen::Matrix<double, 3, 2> out;
  //   std::array<double, 6> jc = jacobianCore(x, y);

  //   for (int j = 0; j < 3; j++) {
  //     out(j, 0) = jc[j];
  //     out(j, 1) = jc[j + 3];
  //   }
  //   return out;
  // }

  inline Eigen::Matrix<double, 3, 1> evaln(const double x,
                                           const double y) const {
    Eigen::Matrix<double, 3, 2> jac = jacobianCore(x, y);
    return jac.col(0).cross(jac.col(1));
  }

#ifdef _BEMBEL_OLD_
  inline Bembel::vector3 f(const Bembel::vector2 &in) const {
    Eigen::Vector3d tmp = evalCore(in.x, in.y);
    return Bembel::vector3_make(tmp(0), tmp(1), tmp(2));
  }
  inline Bembel::vector3 n_f(const Bembel::vector2 &in) const {
    Eigen::Vector3d tmp = evaln(in.x, in.y);
    return Bembel::vector3_make(tmp(0), tmp(1), tmp(2));
  }
  // std::array<double, 6> jac_f(const vector2 &in) const {
  //   return flip ? jacobianCore(in.y, in.x) : jacobianCore(in.x, in.y);
  // }

  inline Bembel::vector3 df_dx(const Bembel::vector2 &in) const {
    const auto tmp = jacobianCore(in.x, in.y);
    return Bembel::vector3_make(tmp(0, 0), tmp(1, 0), tmp(2, 0));
  }
  inline Bembel::vector3 df_dy(const Bembel::vector2 &in) const {
    const auto tmp = jacobianCore(in.x, in.y);
    return Bembel::vector3_make(tmp(0, 1), tmp(1, 1), tmp(2, 1));
  }
  inline std::array<Bembel::vector3, 2> jacobian(
      const Bembel::vector2 &in) const {
    const auto tmp = jacobianCore(in.x, in.y);
    return {Bembel::vector3_make(tmp(0, 0), tmp(1, 0), tmp(2, 0)),
            Bembel::vector3_make(tmp(0, 1), tmp(1, 1), tmp(2, 1))};
  }
#endif
};

}  // namespace Spl
#endif
