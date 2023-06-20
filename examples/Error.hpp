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
#ifndef EXAMPLES_ERROR_HPP_
#define EXAMPLES_ERROR_HPP_

/**
 * @brief Routines for the evalutation of pointwise errors.
 */

namespace Bembel {

template <typename Scalar>
inline Eigen::Matrix<double, Eigen::Dynamic, 1> errors(
    const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> &pot,
    const Eigen::MatrixXd &grid,
    const std::function<Scalar(Eigen::Vector3d)> &fun) {
  const int gridsz = grid.rows();
  assert((std::max(pot.cols(), pot.rows()) == grid.rows()) &&
         ("The size does not match!"));
  assert((grid.cols() == 3) &&
         "The grid must be a Matrix with a 3d point in each row!");
  assert((std::min(pot.rows(), pot.cols())) && ("Potential must be Vector!"));

  Eigen::Matrix<double, 1, Eigen::Dynamic> errors(gridsz);

  for (int i = 0; i < gridsz; i++) {
    errors(i) = std::abs(pot(i) - fun((grid.row(i).transpose()).eval()));
  }
  return errors;
}

inline Eigen::Matrix<double, Eigen::Dynamic, 1> errors(
    const Eigen::MatrixXcd &pot, const Eigen::MatrixXd &grid,
    const std::function<Eigen::Vector3cd(Eigen::Vector3d, std::complex<double>)>
        &fun,
    std::complex<double> kappa) {
  const int gridsz = grid.rows();
  assert((std::max(pot.cols(), pot.rows()) == grid.rows()) &&
         ("The size does not match!"));
  assert((grid.cols() == 3) &&
         "The grid must be a Matrix with a 3d point in each row!");
  assert((pot.cols() == 3) &&
         ("Potential must be a Matrix with a Vector3cd in each col!"));
  Eigen::Matrix<double, 1, Eigen::Dynamic> errors(gridsz);

  for (int i = 0; i < gridsz; i++) {
    errors(i) =
        (pot.col(i) - fun((grid.row(i).transpose()).eval(), kappa)).norm();
  }
  return errors;
}

template <typename Scalar>
inline double maxPointwiseError(
    const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> &pot,
    const Eigen::MatrixXd &grid,
    const std::function<Scalar(Eigen::Vector3d)> &fun) {
  const int gridsz = grid.rows();
  assert((std::max(pot.cols(), pot.rows()) == grid.rows()) &&
         ("The size does not match!"));
  assert((grid.cols() == 3) &&
         "The grid must be a Matrix with a 3d point in each row!");
  assert((std::min(pot.cols(), pot.rows()) == 1) &&
         ("Potential must be a vector!"));
  double error = 0;

  for (int i = 0; i < gridsz; i++) {
    double tmp = std::abs(pot(i) - fun((grid.row(i).transpose()).eval()));
    error = tmp > error ? tmp : error;
  }
  return error;
}

inline double maxPointwiseError(
    const Eigen::MatrixXcd &pot, const Eigen::MatrixXd &grid,
    const std::function<Eigen::Vector3cd(Eigen::Vector3d)> &fun) {
  const int gridsz = grid.rows();
  assert(pot.cols() == grid.cols());
  assert((std::max(pot.cols(), pot.rows()) == grid.rows()) &&
         ("The size does not match!"));
  assert((grid.cols() == 3) &&
         "The grid must be a Matrix with a 3d point in each row!");
  assert((pot.cols() == 3) &&
         ("Must be a Matrix with a point solution in each row!"));
  double error = 0;

  for (int i = 0; i < gridsz; i++) {
    double tmp = (pot.row(i) - fun(grid.row(i)).transpose()).norm();
    error = tmp > error ? tmp : error;
  }
  return error;
}

double estimateRateOfConvergence(const Eigen::VectorXd &errors) {
  Eigen::MatrixXd A(errors.rows(), 2);
  A << Eigen::VectorXd::Ones(errors.rows()),
      Eigen::VectorXd::LinSpaced(errors.rows(), 0, errors.rows() - 1);
  Eigen::VectorXd b = errors.array().abs().log() / std::log(2);
  Eigen::VectorXd x = A.colPivHouseholderQr().solve(b);
  return -x(1);
}

bool checkRateOfConvergence(const Eigen::VectorXd &errors,
                            const int expected_rate, const double tol_factor,
                            double *rate_of_convergence_out = NULL) {
  double rate_of_convergence = estimateRateOfConvergence(errors);
  if (rate_of_convergence_out) rate_of_convergence_out[0] = rate_of_convergence;
  std::cout << "Estimated rate of convergence:" << rate_of_convergence
            << std::endl;
  return (rate_of_convergence > tol_factor * expected_rate);
}

}  // namespace Bembel

#endif  // EXAMPLES_ERROR_HPP_
