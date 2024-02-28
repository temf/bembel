// This file is part of Bembel, the higher order C++ boundary element library.
//
// Copyright (C) 2024 see <http://www.bembel.eu>
//
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.

/**
 * This unit test checks if the integration routines implemented in the
 * HelmholtzSingleLayerOperator and Potential are working correctly.
 */

#include <Bembel/AnsatzSpace>
#include <Bembel/Helmholtz>
#include <Eigen/Dense>

#include "tests/Test.hpp"

int main() {
  using namespace Bembel;

  int polynomial_degree = 0;
  int refinement_level = 0;

  const std::complex<double> wavenumber(1., 0.);

  Test::TestGeometryWriter::writeScreen();
  Bembel::Geometry geometry("test_Screen.dat");

  AnsatzSpace<HelmholtzSingleLayerOperator> ansatz_space(
      geometry, refinement_level, polynomial_degree);

  const SuperSpace<HelmholtzSingleLayerOperator>& super_space =
      ansatz_space.get_superspace();
  const ElementTree& element_tree =
      ansatz_space.get_superspace().get_mesh().get_element_tree();

  // SingleLayerOperator
  {
    HelmholtzSingleLayerOperator linOp;
    linOp.set_wavenumber(wavenumber);

    Eigen::MatrixXcd intval = Eigen::MatrixXd::Zero(1, 1);

    Eigen::Vector2d ref1(0., 0.);
    Eigen::Vector2d ref2(0., 1.);

    SurfacePoint p1, p2;
    super_space.map2surface(*(element_tree.pbegin()), ref1, 1, &p1);
    super_space.map2surface(*(element_tree.pbegin()), ref2, 1, &p2);

    linOp.evaluateIntegrand_impl(super_space, p1, p2, &intval);

    BEMBEL_TEST_IF(
        std::abs(intval(0, 0) -
                 std::exp(-std::complex<double>(0., 1.) * wavenumber) /
                     (4 * M_PI)) < Constants::generic_tolerance);
  }

  // SingleLayerPotential
  {
    HelmholtzSingleLayerPotential<HelmholtzSingleLayerOperator> pot;
    pot.set_wavenumber(wavenumber);

    FunctionEvaluator<HelmholtzSingleLayerOperator> fun_ev(ansatz_space);
    fun_ev.set_function(Eigen::VectorXd::Ones(1));

    Eigen::Vector2d pt(0.5, 0.5);
    Eigen::Vector3d pt_eval(0.5, 0.5, 1.0);

    SurfacePoint sp;
    super_space.map2surface(*(element_tree.pbegin()), pt, 1, &sp);

    auto intval = pot.evaluateIntegrand_impl(fun_ev, *(element_tree.pbegin()),
                                             pt_eval, sp);

    BEMBEL_TEST_IF(
        std::abs(intval(0, 0) -
                 std::exp(-std::complex<double>(0., 1.) * wavenumber) /
                     (4 * M_PI)) < Constants::generic_tolerance);
  }

  return 0;
}
