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

#ifndef BEMBEL_SRC_QUADRATURE_QUADRATUREVECTOR_HPP_
#define BEMBEL_SRC_QUADRATURE_QUADRATUREVECTOR_HPP_

namespace Bembel {

/**
 *  \ingroup Quadrature
 *  \brief this struct wraps all the defined quadrature Rules in a nice
 *         structure overloading the [] operator such that they can
 *         be accessed within a loop during runtime
 **/
template <template <unsigned int qrOrder> class QuadratureRule,
          unsigned int Order>
struct QuadratureVector {
  QuadratureVector() {
    QuadratureRule<Order + 1> QR;
    Q_.xi_ = Eigen::Map<Eigen::VectorXd>(QR.xi_.data(), QR.xi_.size());
    Q_.w_ = Eigen::Map<Eigen::VectorXd>(QR.w_.data(), QR.w_.size());
  }
  QuadratureVector<QuadratureRule, Order - 1> remainingQuadratures_;
  const Quadrature<1> &operator[](unsigned int i) const {
    return (i == Order) ? Q_ : remainingQuadratures_[i];
  }

  Quadrature<1> Q_;
};

template <template <unsigned int qrOrder> class QuadratureRule>
struct QuadratureVector<QuadratureRule, 0> {
  QuadratureVector() {
    QuadratureRule<1> QR;
    Q_.xi_ = Eigen::Map<Eigen::VectorXd>(QR.xi_.data(), QR.xi_.size());
    Q_.w_ = Eigen::Map<Eigen::VectorXd>(QR.w_.data(), QR.w_.size());
  }
  Quadrature<1> Q_;
  const Quadrature<1> &operator[](unsigned int i) const { return Q_; }
};
}  // namespace Bembel
#endif  // BEMBEL_SRC_QUADRATURE_QUADRATUREVECTOR_HPP_
