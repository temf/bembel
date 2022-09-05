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

#include "tests/DuffyTrick/test_DuffyTrick"

#include "tests/TestGeometries.hpp"

int main() {
  Test::TestGeometryWriter::writeScreen();
  Bembel::Geometry geometry("test_Screen.dat");
  Bembel::AnsatzSpace<Bembel::DummyOperator> ansatz_space(
      geometry, Test::Constants::DuffyTrick_refinement_level, 0);
  // test integrate0
  {
    Bembel::DummyOperator Op0(Test::DuffyTrick::integrate0_test_function);
    bool test0 = Test::DuffyTrick::test_integrate0<
        Bembel::DummyOperator, Test::Constants::DuffyTrick_quadrature_degree>(
        ansatz_space, Op0);
    BEMBEL_TEST_IF(test0);
  }
  // test integrate1
  {
    Bembel::DummyOperator Op1(Test::DuffyTrick::integrate1_test_function);
    bool test1 = Test::DuffyTrick::test_integrate1<
        Bembel::DummyOperator, Test::Constants::DuffyTrick_quadrature_degree>(
        ansatz_space, Op1);
    BEMBEL_TEST_IF(test1);
  }
  // test integrate2
  {
    Bembel::DummyOperator Op2(Test::DuffyTrick::integrate2_test_function);
    bool test2 = Test::DuffyTrick::test_integrate2<
        Bembel::DummyOperator, Test::Constants::DuffyTrick_quadrature_degree>(
        ansatz_space, Op2);
    BEMBEL_TEST_IF(test2);
  }
  // test integrate3
  {
    Bembel::DummyOperator Op3(Test::DuffyTrick::integrate3_test_function);
    bool test3 = Test::DuffyTrick::test_integrate3<
        Bembel::DummyOperator, Test::Constants::DuffyTrick_quadrature_degree>(
        ansatz_space, Op3);
    BEMBEL_TEST_IF(test3);
  }
  // test integrate4
  {
    Bembel::DummyOperator Op4(Test::DuffyTrick::integrate4_test_function);
    bool test4 = Test::DuffyTrick::test_integrate4<
        Bembel::DummyOperator, Test::Constants::DuffyTrick_quadrature_degree>(
        ansatz_space, Op4);
    BEMBEL_TEST_IF(test4);
  }

  return 0;
}
