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
//
#ifndef TESTS_DUFFYTRICK_TEST_INTEGRATE3_HPP_
#define TESTS_DUFFYTRICK_TEST_INTEGRATE3_HPP_

namespace Test {
namespace DuffyTrick {

template <typename Derived, unsigned int maxqdeg>
bool test_integrate3(const Bembel::AnsatzSpace<Derived> &ansatz_space,
                     const Bembel::LinearOperatorBase<Derived> &linOp) {
  Bembel::GaussSquare<maxqdeg + 1> GS;
  auto Q = GS[maxqdeg];

  std::vector<SurfacePoint> ffield_qnodes;
  Eigen::MatrixXd intval;
  Eigen::MatrixXd axis;
  intval.resize(1, 1);
  axis.resize(2, 4);

  double maxError = 0;
  double error = 0;
  auto begin =
      ansatz_space.get_superspace().get_mesh().get_element_tree().cpbegin();
  auto end =
      ansatz_space.get_superspace().get_mesh().get_element_tree().cpend();
  double h = begin->get_h();
  // we test all identical elements...

  for (auto it = begin; it != end; ++it)
    for (auto it2 = begin; it2 != end; ++it2) {
      double dist = 0;
      auto cp = Bembel::DuffyTrick::compareElements(*it, *it2, &dist);
      if (cp(2) == 3) {
        intval.setZero();
        error = 0;
        ////////////////////////////////////////////////////////////////////////
        Bembel::DuffyTrick::integrate3(linOp, ansatz_space.get_superspace(),
                                       *it, cp(0), *it2, cp(1), ffield_qnodes,
                                       ffield_qnodes, Q, &intval);
        ////////////////////////////////////////////////////////////////////////
        axis.col(0) << it->llc_(0), it->llc_(0) + h;
        axis.col(1) << it->llc_(1), it->llc_(1) + h;
        axis.col(2) << it2->llc_(0), it2->llc_(0) + h;
        axis.col(3) << it2->llc_(1), it2->llc_(1) + h;
        auto exactInt =
            integrate3_test_function_integral(axis) / BEMBEL_SQUARED_(h);
        ////////////////////////////////////////////////////////////////////////
        error = std::abs(exactInt - intval(0, 0)) / std::abs(exactInt);
        maxError = error > maxError ? error : maxError;
      }
    }
  std::cout << maxError << std::endl;
  return maxError < Constants::integrate3_tolerance;
}
}  // namespace DuffyTrick
}  // namespace Test
#endif  // TESTS_DUFFYTRICK_TEST_INTEGRATE3_HPP_
