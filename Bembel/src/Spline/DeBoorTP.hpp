// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef BEMBEL_SPLINE_DEBOORTP_H_
#define BEMBEL_SPLINE_DEBOORTP_H_

namespace Bembel {
namespace Spl {
/**
 *  \ingroup Spline
 *  \brief A "by the book" implementation of the derivatives and TP-algos based 
 *         on the DeBoor Recursion.
 **/
template <typename T>
Eigen::Matrix<T, -1, -1> DeBoorDer(
    Eigen::Matrix<T, -1, -1> const &control_points,
    std::vector<double> const &knot,
    std::vector<double> const &evaluation_points) noexcept {
  const int control_points_cols = control_points.cols();
  const int polynomial_degree = knot.size() - control_points_cols - 1;
  Eigen::Matrix<T, -1, -1> temp(control_points.rows(), control_points_cols - 1);
  for (int i = control_points_cols - 1; 0 <= --i;) {
    temp.col(i) = (polynomial_degree) /
                  (knot[i + polynomial_degree + 1] - knot[i + 1]) *
                  (control_points.col(i + 1) - control_points.col(i));
  }

  std::vector<double> tempknt(knot.begin() + 1, knot.end() - 1);
  return DeBoor(temp, tempknt, evaluation_points);
}

template <typename T>
std::vector<Eigen::Matrix<T, -1, -1>> DeBoorDer(
    std::vector<Eigen::Matrix<T, -1, -1>> const &control_points,
    std::vector<double> const &knot,
    std::vector<double> const &evaluation_points) noexcept {
  const int control_points_cols = control_points[0].cols();
  const int dimension = control_points.size();
  for (int ll = 0; ll < dimension; ll++) {
    assert(control_points[0].cols() == control_points[ll].cols() &&
           control_points[0].rows() == control_points[ll].rows());
  };
  const int polynomial_degree = knot.size() - control_points_cols - 1;
  std::vector<Eigen::Matrix<T, -1, -1>> temp(dimension);
  for (int ll = 0; ll < dimension; ll++) {
    temp[ll].resize(control_points[ll].rows(), control_points_cols - 1);
  };

  for (int i = control_points_cols - 1; 0 <= --i;) {
    double factor =
        (polynomial_degree) / (knot[i + polynomial_degree + 1] - knot[i + 1]);
    for (int ll = 0; ll < dimension; ll++)
      temp[ll].col(i) =
          factor * (control_points[ll].col(i + 1) - control_points[ll].col(i));
  }

  std::vector<double> tempknt(knot.begin() + 1, knot.end() - 1);

  return DeBoor(temp, tempknt, evaluation_points);
}

template <typename T>
std::vector<Eigen::Matrix<T, -1, -1>> deBoorDerGiveData(
    std::vector<Eigen::Matrix<T, -1, -1>> const &control_points,
    std::vector<double> const &knot) noexcept {
  const int control_points_cols = control_points[0].cols();
  const int dimension = control_points.size();
  for (int ll = 0; ll < dimension; ll++) {
    assert(control_points[0].cols() == control_points[ll].cols() &&
           control_points[0].rows() == control_points[ll].rows());
  };
  const int polynomial_degree = knot.size() - control_points_cols - 1;
  std::vector<Eigen::Matrix<T, -1, -1>> temp(dimension);
  for (int ll = 0; ll < dimension; ll++) {
    temp[ll].resize(control_points[ll].rows(), control_points_cols - 1);
  };

  for (int i = control_points_cols - 1; 0 <= --i;) {
    double factor =
        (polynomial_degree) / (knot[i + polynomial_degree + 1] - knot[i + 1]);
    for (int ll = 0; ll < dimension; ll++)
      temp[ll].col(i) =
          factor * (control_points[ll].col(i + 1) - control_points[ll].col(i));
  }

  return temp;
}
////////////////////////////////////////////////////////////////////////////////
/// Simple TP Bsplines
////////////////////////////////////////////////////////////////////////////////
template <typename T>
Eigen::Matrix<T, -1, -1> DeBoorTP(
    Eigen::Matrix<T, -1, -1> const &control_points,
    std::vector<double> const &knots_x, std::vector<double> const &knots_y,
    std::vector<double> const &evaluation_points_x,
    std::vector<double> const &evaluation_points_y) noexcept {
  Eigen::Matrix<T, -1, -1> tmp =
      DeBoor(control_points, knots_x, evaluation_points_x).transpose();
  return DeBoor(tmp, knots_y, evaluation_points_y);
}

template <typename T>
std::vector<Eigen::Matrix<T, -1, -1>> DeBoorTP(
    std::vector<Eigen::Matrix<T, -1, -1>> const &control_points,
    std::vector<double> const &knots_x, std::vector<double> const &knots_y,
    std::vector<double> const &evaluation_points_x,
    std::vector<double> const &evaluation_points_y) noexcept {
  std::vector<Eigen::Matrix<T, -1, -1>> tmp =
      DeBoor(control_points, knots_x, evaluation_points_x);
  for (int ll = tmp.size() - 1; ll >= 0; ll--)
    tmp[ll] = tmp[ll].transpose().eval();
  return DeBoor(tmp, knots_y, evaluation_points_y);
}
////////////////////////////////////////////////////////////////////////////////
// TP Der with direction declaration
////////////////////////////////////////////////////////////////////////////////
template <typename T>
Eigen::Matrix<T, -1, -1> DeBoorTPDer(
    Eigen::Matrix<T, -1, -1> const &control_points,
    std::vector<double> const &knots_x, std::vector<double> const &knots_y,
    std::vector<double> const &evaluation_points_x,
    std::vector<double> const &evaluation_points_y, bool const &x_to_be_derived,
    bool const &y_to_be_derived) noexcept {
  Eigen::Matrix<T, -1, -1> tmp =
      x_to_be_derived
          ? DeBoorDer(control_points, knots_x, evaluation_points_x).transpose()
          : DeBoor(control_points, knots_x, evaluation_points_x).transpose();
  return y_to_be_derived ? DeBoorDer(tmp, knots_y, evaluation_points_y)
                         : DeBoor(tmp, knots_y, evaluation_points_y);
}

template <typename T>
std::vector<Eigen::Matrix<T, -1, -1>> DeBoorTPDer(
    std::vector<Eigen::Matrix<T, -1, -1>> const &control_points,
    std::vector<double> const &knots_x, std::vector<double> const &knots_y,
    std::vector<double> const &evaluation_points_x,
    std::vector<double> const &evaluation_points_y, bool const &x_to_be_derived,
    bool const &y_to_be_derived) noexcept {
  std::vector<Eigen::Matrix<T, -1, -1>> tmp =
      x_to_be_derived ? DeBoorDer(control_points, knots_x, evaluation_points_x)
                      : DeBoor(control_points, knots_x, evaluation_points_x);
  for (int ll = tmp.size() - 1; ll >= 0; ll--)
    tmp[ll] = tmp[ll].transpose().eval();
  return y_to_be_derived ? DeBoorDer(tmp, knots_y, evaluation_points_y)
                         : DeBoor(tmp, knots_y, evaluation_points_y);
}
}  // namespace Spl
}  // namespace Bembel
#endif
