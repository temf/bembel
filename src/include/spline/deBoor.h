// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef _include_deBoor_
#define _include_deBoor_

#include <assert.h>
#include <Eigen/Dense>
#include <vector>
namespace Spl {
//"By the book" implementations of the Cox-DeBoor formula. Inefficient, do not
// use at bottlenecks.
template <typename T>
Eigen::Matrix<T, -1, -1> deBoor(Eigen::Matrix<T, -1, -1> const &ctrl,
                                const std::vector<double> &knots,
                                const std::vector<double> &eval) noexcept {
  const int lengthCtrl = ctrl.cols();
  const int depthCtrl = ctrl.rows();
  const int poldeg = knots.size() - lengthCtrl - 1;
  const int numEval = eval.size();
  assert(poldeg >= 0);
  Eigen::Matrix<T, -1, -1> out(depthCtrl, numEval);

  for (int j = numEval - 1; j >= 0; j--) {
    const double &x = eval[j];

    assert(x >= knots[poldeg] && x <= knots[lengthCtrl + poldeg]);
    int l = 0;
    {
      while (knots[l] <= x && l != lengthCtrl) l++;
      l = l - 1;
    }

    assert(l >= 0);
    Eigen::Matrix<T, -1, -1> temp =
        ctrl.block(0, l - poldeg, depthCtrl, poldeg + 1);

    // ISO C++ forbids variable length array
    assert(poldeg <= 18);
    double ws[18];
    // Iterators remain the same size, order is correct
    for (int k = poldeg; 0 != k; k--) {
      for (int i = 0; i < k; i++) {
        ws[i] = (x - knots[l - k + 1 + i]) /
                (knots[l + 1 + i] - knots[l - k + 1 + i]);
      }
      for (int i = 0; i < k; i++) {
        temp.col(i) =
            ((ws[i] * temp.col(1 + i)) + ((1 - ws[i]) * temp.col(i))).eval();
      }
    }
    out.col(j) = temp.col(0);
  }
  return out;
}

template <typename T>
std::vector<T> deBoor(std::vector<T> const &ctrl,
                      const std::vector<double> &knots,
                      const std::vector<double> &eval) noexcept {
  const int lengthCtrl = ctrl.size();
  const int poldeg = knots.size() - lengthCtrl - 1;
  const int numEval = eval.size();
  assert(poldeg >= 0);
  std::vector<T> out(numEval);

  for (int j = numEval - 1; j >= 0; j--) {
    const double &x = eval[j];

    assert(x >= knots[poldeg] && x <= knots[lengthCtrl + poldeg]);
    int l = 0;
    {
      while (knots[l] <= x && l != lengthCtrl) l++;
      l = l - 1;
    }

    assert(l >= 0);

    std::vector<T> temp(ctrl.begin() + l - poldeg, ctrl.begin() + l + 1);

    // ISO C++ forbids variable length array
    assert(poldeg <= 19 &&
           "Polynomial degree not supprted. See spline.h and polynom.h");
    double ws[20];
    // Iterators remain the same size, order is correct
    for (int k = poldeg; 0 != k; k--) {
      for (int i = 0; i < k; i++) {
        ws[i] = (x - knots[l - k + 1 + i]) /
                (knots[l + 1 + i] - knots[l - k + 1 + i]);
      }
      for (int i = 0; i < k; i++) {
        temp[i] = ((ws[i] * temp[1 + i]) + ((1 - ws[i]) * temp[i]));
      }
    }
    out[j] = temp[0];
  }
  return out;
}

template <typename T>
std::vector<Eigen::Matrix<T, -1, -1>> deBoor(
    std::vector<Eigen::Matrix<T, -1, -1>> const &ctrl,
    std::vector<double> const &knots,
    std::vector<double> const &eval) noexcept {
  const int dim = ctrl.size();
  const int lengthCtrl = ctrl[0].cols();
  const int depthCtrl = ctrl[0].rows();
  const int poldeg = knots.size() - lengthCtrl - 1;
  const int numEval = eval.size();
  assert(poldeg >= 0);
  std::vector<Eigen::Matrix<T, -1, -1>> out(dim);

  for (int l = 0; l < dim; l++) out[l].resize(depthCtrl, numEval);

  for (int j = numEval - 1; j >= 0; j--) {
    const double &x = eval[j];

    assert(x >= knots[poldeg] && x <= knots[lengthCtrl + poldeg]);
    int l = 0;
    {
      while (knots[l] <= x && l != lengthCtrl) l++;
      l = l - 1;
    }

    assert(l >= 0);
    std::vector<Eigen::Matrix<T, -1, -1>> temp(dim);
    for (int h = 0; h < dim; h++)
      temp[h] = ctrl[h].block(0, l - poldeg, depthCtrl, poldeg + 1);

    // ISO C++ forbids variable length array
    assert(poldeg <= 18);
    double ws[18];
    // Iterators remain the same size, order is correct
    for (int k = poldeg; 0 != k; k--) {
      for (int i = 0; i < k; i++) {
        ws[i] = (x - knots[l - k + 1 + i]) /
                (knots[l + 1 + i] - knots[l - k + 1 + i]);
      }

      for (int h = 0; h < dim; h++) {
        for (int i = 0; i < k; i++) {
          temp[h].col(i) =
              ((ws[i] * temp[h].col(1 + i)) + ((1 - ws[i]) * temp[h].col(i)))
                  .eval();
        }
      }
    }
    for (int h = 0; h < dim; h++) out[h].col(j) = temp[h].col(0);
  }
  return out;
}
}  // namespace Spl
#endif
