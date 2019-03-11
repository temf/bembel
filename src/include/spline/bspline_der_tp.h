// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef _BSPLN_DERTP_INCLUDED_
#define _BSPLN_DERTP_INCLUDED_
#include "BSpline/bspline_eval.h"

// using namespace std;
// using namespace Eigen;
namespace Spl {

// A "by the book" implementation of the derivatives and TP-algos based on the
// DeBoor Recursion.

// Simple Derivative

template <typename T>
Matrix<T, -1, -1> deBoorDer(Eigen::Matrix<T, -1, -1> const &ctrl,
                            std::vector<double> const &knot,
                            vector<double> const &eval) noexcept {
  // constexpr int -1 = -1 == -1 ? -1 : -1 - 1;
  // assert( -1 == -1 || -1 > 0 );
  int lengthCtrl = ctrl.cols();
  int poldeg = knot.size() - lengthCtrl - 1;
  Eigen::Matrix<T, -1, -1> temp(ctrl.rows(), lengthCtrl - 1);
  // assert( lengthCtrl - 1 + ( poldeg + 1 ) == (poldeg - 1 + lengthCtrl) );
  for (int i = lengthCtrl - 1; 0 <= --i;) {
    temp.col(i) = (poldeg) / (knot[i + poldeg + 1] - knot[i + 1]) *
                  (ctrl.col(i + 1) - ctrl.col(i));
  }

  std::vector<double> tempknt(knot.begin() + 1, knot.end() - 1);
  return deBoor(temp, tempknt, eval);
}

template <typename T>
vector<Matrix<T, -1, -1>> deBoorDer(
    std::vector<Eigen::Matrix<T, -1, -1>> const &ctrl,
    std::vector<double> const &knot, std::vector<double> const &eval) noexcept {
  // constexpr int -1 = -1 == -1 ? -1 : -1 - 1;
  // assert( -1 == -1 || -1 > 0 );
  int lengthCtrl = ctrl[0].cols();
  int DIM = ctrl.size();
  for (int ll = 0; ll < DIM; ll++) {
    assert(ctrl[0].cols() == ctrl[ll].cols() &&
           ctrl[0].rows() == ctrl[ll].rows());
  };
  int poldeg = knot.size() - lengthCtrl - 1;
  std::vector<Eigen::Matrix<T, -1, -1>> temp(DIM);
  for (int ll = 0; ll < DIM; ll++) {
    temp[ll].resize(ctrl[ll].rows(), lengthCtrl - 1);
  };
  // assert( lengthCtrl - 1 + ( poldeg + 1 ) == (poldeg - 1 + lengthCtrl) );

  for (int i = lengthCtrl - 1; 0 <= --i;) {
    double factor = (poldeg) / (knot[i + poldeg + 1] - knot[i + 1]);
    for (int ll = 0; ll < DIM; ll++)
      temp[ll].col(i) = factor * (ctrl[ll].col(i + 1) - ctrl[ll].col(i));
  }

  std::vector<double> tempknt(knot.begin() + 1, knot.end() - 1);

  return deBoor(temp, tempknt, eval);
}

template <typename T>
vector<Matrix<T, -1, -1>> deBoorDerGiveData(
    std::vector<Matrix<T, -1, -1>> const &ctrl,
    std::vector<double> const &knot) noexcept {
  int lengthCtrl = ctrl[0].cols();
  int DIM = ctrl.size();
  for (int ll = 0; ll < DIM; ll++) {
    assert(ctrl[0].cols() == ctrl[ll].cols() &&
           ctrl[0].rows() == ctrl[ll].rows());
  };
  int poldeg = knot.size() - lengthCtrl - 1;
  std::vector<Eigen::Matrix<T, -1, -1>> temp(DIM);
  for (int ll = 0; ll < DIM; ll++) {
    temp[ll].resize(ctrl[ll].rows(), lengthCtrl - 1);
  };
  // assert( lengthCtrl - 1 + ( poldeg + 1 ) == (poldeg - 1 + lengthCtrl) );

  for (int i = lengthCtrl - 1; 0 <= --i;) {
    double factor = (poldeg) / (knot[i + poldeg + 1] - knot[i + 1]);
    for (int ll = 0; ll < DIM; ll++)
      temp[ll].col(i) = factor * (ctrl[ll].col(i + 1) - ctrl[ll].col(i));
  }

  return temp;
}

// Simple TP Bsplines

template <typename T>
Eigen::Matrix<T, -1, -1> deBoorTP(Eigen::Matrix<T, -1, -1> const &ctrl,
                                  std::vector<double> const &knotx,
                                  std::vector<double> const &knoty,
                                  std::vector<double> const &evalx,
                                  std::vector<double> const &evaly) noexcept {
  Eigen::Matrix<T, -1, -1> tmp = deBoor(ctrl, knotx, evalx).transpose();
  return deBoor(tmp, knoty, evaly);
}

template <typename T>
std::vector<Eigen::Matrix<T, -1, -1>> deBoorTP(
    std::vector<Matrix<T, -1, -1>> const &ctrl,
    std::vector<double> const &knotx, std::vector<double> const &knoty,
    std::vector<double> const &evalx,
    std::vector<double> const &evaly) noexcept {
  std::vector<Matrix<T, -1, -1>> tmp = deBoor(ctrl, knotx, evalx);
  for (int ll = tmp.size() - 1; ll >= 0; ll--)
    tmp[ll] = tmp[ll].transpose().eval();
  return deBoor(tmp, knoty, evaly);
}

// TP Der with direction declaration

template <typename T>
Eigen::Matrix<T, -1, -1> deBoorTPDer(Eigen::Matrix<T, -1, -1> const &ctrl,
                                     std::vector<double> const &knotx,
                                     std::vector<double> const &knoty,
                                     std::vector<double> const &evalx,
                                     std::vector<double> const &evaly,
                                     bool const &xder,
                                     bool const &yder) noexcept {
  Eigen::Matrix<T, -1, -1> tmp = xder
                                     ? deBoorDer(ctrl, knotx, evalx).transpose()
                                     : deBoor(ctrl, knotx, evalx).transpose();
  return yder ? deBoorDer(tmp, knoty, evaly) : deBoor(tmp, knoty, evaly);
}

template <typename T>
std::vector<Eigen::Matrix<T, -1, -1>> deBoorTPDer(
    std::vector<Matrix<T, -1, -1>> const &ctrl,
    std::vector<double> const &knotx, std::vector<double> const &knoty,
    std::vector<double> const &evalx, std::vector<double> const &evaly,
    bool const &xder, bool const &yder) noexcept {
  std::vector<Matrix<T, -1, -1>> tmp =
      xder ? deBoorDer(ctrl, knotx, evalx) : deBoor(ctrl, knotx, evalx);
  for (int ll = tmp.size() - 1; ll >= 0; ll--)
    tmp[ll] = tmp[ll].transpose().eval();
  return yder ? deBoorDer(tmp, knoty, evaly) : deBoor(tmp, knoty, evaly);
}

}  // namespace Spl
#endif