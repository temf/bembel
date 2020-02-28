// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef BEMBEL_SPLINE_DEBOOR_H_
#define BEMBEL_SPLINE_DEBOOR_H_

namespace Bembel {
namespace Spl {
/**
 *  \ingroup Spline
 *  \brief "By the book" implementations of the Cox-DeBoor formula. Inefficient,
 *         do not use at bottlenecks.
 **/
template <typename T>
Eigen::Matrix<T, -1, -1> DeBoor(
    Eigen::Matrix<T, -1, -1> const &control_points,
    const std::vector<double> &knot_vector,
    const std::vector<double> &evaluation_points) noexcept {
  const int cols_control_points = control_points.cols();
  const int rows_control_points = control_points.rows();
  const int polynomial_degree = knot_vector.size() - cols_control_points - 1;
  const int size_evaluation_points = evaluation_points.size();
  assert(polynomial_degree >= 0);
  Eigen::Matrix<T, -1, -1> out(rows_control_points, size_evaluation_points);

  for (int j = size_evaluation_points - 1; j >= 0; j--) {
    const double &x = evaluation_points[j];

    assert(x >= knot_vector[polynomial_degree] &&
           x <= knot_vector[cols_control_points + polynomial_degree]);
    int l = 0;
    {
      while (knot_vector[l] <= x && l != cols_control_points) l++;
      l = l - 1;
    }

    assert(l >= 0);
    Eigen::Matrix<T, -1, -1> temp = control_points.block(
        0, l - polynomial_degree, rows_control_points, polynomial_degree + 1);

    /// ISO C++ forbids variable length array
    assert(polynomial_degree <= 18);
    double ws[18];
    // Iterators remain the same size, order is correct
    for (int k = polynomial_degree; 0 != k; k--) {
      for (int i = 0; i < k; i++) {
        ws[i] = (x - knot_vector[l - k + 1 + i]) /
                (knot_vector[l + 1 + i] - knot_vector[l - k + 1 + i]);
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
std::vector<T> DeBoor(std::vector<T> const &control_points,
                      const std::vector<double> &knot_vector,
                      const std::vector<double> &evaluation_points) noexcept {
  const int cols_control_points = control_points.size();
  const int polynomial_degree = knot_vector.size() - cols_control_points - 1;
  const int size_evaluation_points = evaluation_points.size();
  assert(polynomial_degree >= 0);
  std::vector<T> out(size_evaluation_points);

  for (int j = size_evaluation_points - 1; j >= 0; j--) {
    const double &x = evaluation_points[j];

    assert(x >= knot_vector[polynomial_degree] &&
           x <= knot_vector[cols_control_points + polynomial_degree]);
    int l = 0;
    {
      while (knot_vector[l] <= x && l != cols_control_points) l++;
      l = l - 1;
    }

    assert(l >= 0);

    std::vector<T> temp(control_points.begin() + l - polynomial_degree,
                        control_points.begin() + l + 1);

    // ISO C++ forbids variable length array
    assert(polynomial_degree <= 19 &&
           "Polynomial degree not supported in DeBoor.");
    double ws[20];
    // Iterators remain the same size, order is correct
    for (int k = polynomial_degree; 0 != k; k--) {
      for (int i = 0; i < k; i++) {
        ws[i] = (x - knot_vector[l - k + 1 + i]) /
                (knot_vector[l + 1 + i] - knot_vector[l - k + 1 + i]);
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
std::vector<Eigen::Matrix<T, -1, -1>> DeBoor(
    std::vector<Eigen::Matrix<T, -1, -1>> const &control_points,
    std::vector<double> const &knot_vector,
    std::vector<double> const &evaluation_points) noexcept {
  const int dimension = control_points.size();
  const int cols_control_points = control_points[0].cols();
  const int rows_control_points = control_points[0].rows();
  const int polynomial_degree = knot_vector.size() - cols_control_points - 1;
  const int size_evaluation_points = evaluation_points.size();
  assert(polynomial_degree >= 0);
  std::vector<Eigen::Matrix<T, -1, -1>> out(dimension);

  for (int l = 0; l < dimension; l++)
    out[l].resize(rows_control_points, size_evaluation_points);

  for (int j = size_evaluation_points - 1; j >= 0; j--) {
    const double &x = evaluation_points[j];

    assert(x >= knot_vector[polynomial_degree] &&
           x <= knot_vector[cols_control_points + polynomial_degree]);
    int l = 0;
    {
      while (knot_vector[l] <= x && l != cols_control_points) l++;
      l = l - 1;
    }

    assert(l >= 0);
    std::vector<Eigen::Matrix<T, -1, -1>> temp(dimension);
    for (int h = 0; h < dimension; h++)
      temp[h] = control_points[h].block(
          0, l - polynomial_degree, rows_control_points, polynomial_degree + 1);

    // ISO C++ forbids variable length array
    assert(polynomial_degree <= 18);
    double ws[18];
    // Iterators remain the same size, order is correct
    for (int k = polynomial_degree; 0 != k; k--) {
      for (int i = 0; i < k; i++) {
        ws[i] = (x - knot_vector[l - k + 1 + i]) /
                (knot_vector[l + 1 + i] - knot_vector[l - k + 1 + i]);
      }

      for (int h = 0; h < dimension; h++) {
        for (int i = 0; i < k; i++) {
          temp[h].col(i) =
              ((ws[i] * temp[h].col(1 + i)) + ((1 - ws[i]) * temp[h].col(i)))
                  .eval();
        }
      }
    }
    for (int h = 0; h < dimension; h++) out[h].col(j) = temp[h].col(0);
  }
  return out;
}
}  // namespace Spl
}  // namespace Bembel
#endif
