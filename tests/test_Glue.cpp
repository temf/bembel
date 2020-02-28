// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.

#include "Test.hpp"
#include <Bembel/AnsatzSpace>

class TestOperatorDivC;
template <>
struct Bembel::LinearOperatorTraits<TestOperatorDivC> {
  typedef Eigen::VectorXd EigenType;
  typedef Eigen::VectorXd::Scalar Scalar;
  enum { OperatorOrder = 0, Form = DifferentialForm::DivConforming };
};

class TestOperatorC;
template <>
struct Bembel::LinearOperatorTraits<TestOperatorC> {
  typedef Eigen::VectorXd EigenType;
  typedef Eigen::VectorXd::Scalar Scalar;
  enum { OperatorOrder = 0, Form = DifferentialForm::Continuous };
};

int main() {
  /* Idea of this test: We want to test the divergence-conformity of the global
     space. For this, the normal component must be continuously glued. Since
     only geometries with uniform orientation are allowed, there are 16 allowed
     edge-to-edge combinations. We will induce all of them via four geometries,
     which look like this:


  //
  //      ----
  //      |4 |
  //   ----------
  //   |2 |0 |1 |
  //   ----------
  //      |3 |
  //      ----
  //


     The 0-patch is rotated in a different way in each of the upcoming
  parametrisations.

  */
  using namespace Bembel;

  Test::TestGeometryWriter::writeEdgeCase1();
  Test::TestGeometryWriter::writeEdgeCase2();
  Test::TestGeometryWriter::writeEdgeCase3();
  Test::TestGeometryWriter::writeEdgeCase4();

  Bembel::Geometry g1("test_EdgeCase1.dat");
  Bembel::Geometry g2("test_EdgeCase2.dat");
  Bembel::Geometry g3("test_EdgeCase3.dat");
  Bembel::Geometry g4("test_EdgeCase4.dat");

  // We check that the orientation is correct, i.e., that all allowed cases are
  // present and no mistake has been made in the rotation of the center-patch.
  for (int i = 0; i < 5; ++i) {
    assert(g1.get_geometry()[i].evalNormal(Eigen::Vector2d(.5, .5))(2) > 0);
    assert(g2.get_geometry()[i].evalNormal(Eigen::Vector2d(.5, .5))(2) > 0);
    assert(g3.get_geometry()[i].evalNormal(Eigen::Vector2d(.5, .5))(2) > 0);
    assert(g4.get_geometry()[i].evalNormal(Eigen::Vector2d(.5, .5))(2) > 0);
  }

  // Now we test for several discretisations.
  for (int P = 1; P < 5; ++P) {
    for (int M = 0; M < 4; ++M) {
      // We initialize the points on which the vektor field will be evaluated
      // along the edges. These are in the middle of each element. Note that
      // testing on element boundaries is dangerous: The normal component ist
      // NOT contiuous along the boundary, so the evaluation can jump and
      // trigger a false positive.
      std::vector<double> pts;
      const int elements = (1 << M);
      const double h = 1. / double(elements);
      for (int i = 0; i < elements; ++i) {
        pts.push_back((i + 0.5) * h);
      }
      pts.shrink_to_fit();

      AnsatzSpace<TestOperatorDivC> a1(g1, M, P);
      AnsatzSpace<TestOperatorDivC> a2(g2, M, P);
      AnsatzSpace<TestOperatorDivC> a3(g3, M, P);
      AnsatzSpace<TestOperatorDivC> a4(g4, M, P);

      AnsatzSpace<TestOperatorC> a1c(g1, M, P);
      AnsatzSpace<TestOperatorC> a2c(g2, M, P);
      AnsatzSpace<TestOperatorC> a3c(g3, M, P);
      AnsatzSpace<TestOperatorC> a4c(g4, M, P);

      const int dofsDiv = a1.get_number_of_dofs();
      const int dofsCon = a1c.get_number_of_dofs();

      // Each dof gets a unique coefficient, such that we can exclude symmetry
      // effects and the like.
      Eigen::VectorXd coefsDiv(dofsDiv);
      for (int i = 0; i < dofsDiv; ++i) {
        coefsDiv(i) = i - dofsDiv / 2;
      }
      Eigen::VectorXd coefsCon(dofsCon);
      for (int i = 0; i < dofsCon; ++i) {
        coefsCon(i) = i - dofsCon / 2;
      }

      FunctionEvaluator<TestOperatorDivC> fe1(a1, coefsDiv);
      FunctionEvaluator<TestOperatorDivC> fe2(a2, coefsDiv);
      FunctionEvaluator<TestOperatorDivC> fe3(a3, coefsDiv);
      FunctionEvaluator<TestOperatorDivC> fe4(a4, coefsDiv);

      FunctionEvaluator<TestOperatorC> fe1c(a1c, coefsCon);
      FunctionEvaluator<TestOperatorC> fe2c(a2c, coefsCon);
      FunctionEvaluator<TestOperatorC> fe3c(a3c, coefsCon);
      FunctionEvaluator<TestOperatorC> fe4c(a4c, coefsCon);

      for (int j = 0; j < pts.size(); ++j) {
        // We check if the normal component of the vector field is continuous.
        // First, we define the points on which geometries will be evaluated.

        const double pt_p = pts[j];
        const double pt_n = 1.0 - pt_p;

        const Eigen::Vector2d e1p(pt_p, 0);
        const Eigen::Vector2d e2p(1, pt_p);
        const Eigen::Vector2d e3p(pt_p, 1);
        const Eigen::Vector2d e4p(0, pt_p);

        const Eigen::Vector2d e1n(pt_n, 0);
        const Eigen::Vector2d e2n(1, pt_n);
        const Eigen::Vector2d e3n(pt_n, 1);
        const Eigen::Vector2d e4n(0, pt_n);

        // Here we check that the points we identify in physical space are
        // indeed correct.
        BEMBEL_TEST_IF(
            (g1.get_geometry()[0].eval(e2p) - g1.get_geometry()[1].eval(e4p))
                .norm() < Test::Constants::test_tolerance_geometry);
        BEMBEL_TEST_IF(
            (g1.get_geometry()[0].eval(e4p) - g1.get_geometry()[2].eval(e2p))
                .norm() < Test::Constants::test_tolerance_geometry);
        BEMBEL_TEST_IF(
            (g1.get_geometry()[0].eval(e1p) - g1.get_geometry()[3].eval(e3p))
                .norm() < Test::Constants::test_tolerance_geometry);
        BEMBEL_TEST_IF(
            (g1.get_geometry()[0].eval(e3p) - g1.get_geometry()[4].eval(e1p))
                .norm() < Test::Constants::test_tolerance_geometry);

        BEMBEL_TEST_IF(
            (g2.get_geometry()[0].eval(e4n) - g2.get_geometry()[1].eval(e4p))
                .norm() < Test::Constants::test_tolerance_geometry);
        BEMBEL_TEST_IF(
            (g2.get_geometry()[0].eval(e2n) - g2.get_geometry()[2].eval(e2p))
                .norm() < Test::Constants::test_tolerance_geometry);
        BEMBEL_TEST_IF(
            (g2.get_geometry()[0].eval(e3n) - g2.get_geometry()[3].eval(e3p))
                .norm() < Test::Constants::test_tolerance_geometry);
        BEMBEL_TEST_IF(
            (g2.get_geometry()[0].eval(e1n) - g2.get_geometry()[4].eval(e1p))
                .norm() < Test::Constants::test_tolerance_geometry);

        BEMBEL_TEST_IF(
            (g3.get_geometry()[0].eval(e1p) - g3.get_geometry()[1].eval(e4p))
                .norm() < Test::Constants::test_tolerance_geometry);
        BEMBEL_TEST_IF(
            (g3.get_geometry()[0].eval(e3p) - g3.get_geometry()[2].eval(e2p))
                .norm() < Test::Constants::test_tolerance_geometry);
        BEMBEL_TEST_IF(
            (g3.get_geometry()[0].eval(e4n) - g3.get_geometry()[3].eval(e3p))
                .norm() < Test::Constants::test_tolerance_geometry);
        BEMBEL_TEST_IF(
            (g3.get_geometry()[0].eval(e2n) - g3.get_geometry()[4].eval(e1p))
                .norm() < Test::Constants::test_tolerance_geometry);

        BEMBEL_TEST_IF(
            (g4.get_geometry()[0].eval(e3n) - g4.get_geometry()[1].eval(e4p))
                .norm() < Test::Constants::test_tolerance_geometry);
        BEMBEL_TEST_IF(
            (g4.get_geometry()[0].eval(e1n) - g4.get_geometry()[2].eval(e2p))
                .norm() < Test::Constants::test_tolerance_geometry);
        BEMBEL_TEST_IF(
            (g4.get_geometry()[0].eval(e2p) - g4.get_geometry()[3].eval(e3p))
                .norm() < Test::Constants::test_tolerance_geometry);
        BEMBEL_TEST_IF(
            (g4.get_geometry()[0].eval(e4p) - g4.get_geometry()[4].eval(e1p))
                .norm() < Test::Constants::test_tolerance_geometry);

        BEMBEL_TEST_IF((std::abs((fe1c.evaluateOnPatch(0, e2p)(0) -
                                  fe1c.evaluateOnPatch(1, e4p)(0))) <
                        Test::Constants::test_tolerance_geometry))
        BEMBEL_TEST_IF((std::abs((fe1c.evaluateOnPatch(0, e4p)(0) -
                                  fe1c.evaluateOnPatch(2, e2p)(0))) <
                        Test::Constants::test_tolerance_geometry))
        BEMBEL_TEST_IF((std::abs((fe1c.evaluateOnPatch(0, e1p)(0) -
                                  fe1c.evaluateOnPatch(3, e3p)(0))) <
                        Test::Constants::test_tolerance_geometry))
        BEMBEL_TEST_IF((std::abs((fe1c.evaluateOnPatch(0, e3p)(0) -
                                  fe1c.evaluateOnPatch(4, e1p)(0))) <
                        Test::Constants::test_tolerance_geometry))
        //------------------
        BEMBEL_TEST_IF((std::abs((fe2c.evaluateOnPatch(0, e4n)(0) -
                                  fe2c.evaluateOnPatch(1, e4p)(0))) <
                        Test::Constants::test_tolerance_geometry))
        // std::cout << (fe2c.evaluateOnPatch(0, e4n) - fe2c.evaluateOnPatch(1,
        // e4p)).transpose() << std::endl;

        BEMBEL_TEST_IF((std::abs((fe2c.evaluateOnPatch(0, e2n)(0) -
                                  fe2c.evaluateOnPatch(2, e2p)(0))) <
                        Test::Constants::test_tolerance_geometry))
        BEMBEL_TEST_IF((std::abs((fe2c.evaluateOnPatch(0, e3n)(0) -
                                  fe2c.evaluateOnPatch(3, e3p)(0))) <
                        Test::Constants::test_tolerance_geometry))
        BEMBEL_TEST_IF((std::abs((fe2c.evaluateOnPatch(0, e1n)(0) -
                                  fe2c.evaluateOnPatch(4, e1p)(0))) <
                        Test::Constants::test_tolerance_geometry))
        //------------------
        BEMBEL_TEST_IF((std::abs((fe3c.evaluateOnPatch(0, e1p)(0) -
                                  fe3c.evaluateOnPatch(1, e4p)(0))) <
                        Test::Constants::test_tolerance_geometry))
        BEMBEL_TEST_IF((std::abs((fe3c.evaluateOnPatch(0, e3p)(0) -
                                  fe3c.evaluateOnPatch(2, e2p)(0))) <
                        Test::Constants::test_tolerance_geometry))
        BEMBEL_TEST_IF((std::abs((fe3c.evaluateOnPatch(0, e4n)(0) -
                                  fe3c.evaluateOnPatch(3, e3p)(0))) <
                        Test::Constants::test_tolerance_geometry))
        BEMBEL_TEST_IF((std::abs((fe3c.evaluateOnPatch(0, e2n)(0) -
                                  fe3c.evaluateOnPatch(4, e1p)(0))) <
                        Test::Constants::test_tolerance_geometry))
        //------------------
        BEMBEL_TEST_IF((std::abs((fe4c.evaluateOnPatch(0, e3n)(0) -
                                  fe4c.evaluateOnPatch(1, e4p)(0))) <
                        Test::Constants::test_tolerance_geometry))
        BEMBEL_TEST_IF((std::abs((fe4c.evaluateOnPatch(0, e1n)(0) -
                                  fe4c.evaluateOnPatch(2, e2p)(0))) <
                        Test::Constants::test_tolerance_geometry))
        BEMBEL_TEST_IF((std::abs((fe4c.evaluateOnPatch(0, e2p)(0) -
                                  fe4c.evaluateOnPatch(3, e3p)(0))) <
                        Test::Constants::test_tolerance_geometry))
        BEMBEL_TEST_IF((std::abs((fe4c.evaluateOnPatch(0, e4p)(0) -
                                  fe4c.evaluateOnPatch(4, e1p)(0))) <
                        Test::Constants::test_tolerance_geometry))

        // The div-conforming glue must yield continuity of the normal
        // component, i.e., for out above-defined geometries one of the
        // following two unit vectors.
        const Eigen::Vector3d uv1(1, 0, 0);
        const Eigen::Vector3d uv2(0, 1, 0);

        // Now, we evaluate the vector field on the matching points above and
        // compare the normal components.
        BEMBEL_TEST_IF((
            std::abs((fe1.evaluateOnPatch(0, e2p) - fe1.evaluateOnPatch(1, e4p))
                         .dot(uv1)) < Test::Constants::test_tolerance_geometry))
        BEMBEL_TEST_IF((
            std::abs((fe1.evaluateOnPatch(0, e4p) - fe1.evaluateOnPatch(2, e2p))
                         .dot(uv1)) < Test::Constants::test_tolerance_geometry))
        BEMBEL_TEST_IF((
            std::abs((fe1.evaluateOnPatch(0, e1p) - fe1.evaluateOnPatch(3, e3p))
                         .dot(uv2)) < Test::Constants::test_tolerance_geometry))
        BEMBEL_TEST_IF((
            std::abs((fe1.evaluateOnPatch(0, e3p) - fe1.evaluateOnPatch(4, e1p))
                         .dot(uv2)) < Test::Constants::test_tolerance_geometry))
        //------------------
        BEMBEL_TEST_IF((
            std::abs((fe2.evaluateOnPatch(0, e4n) - fe2.evaluateOnPatch(1, e4p))
                         .dot(uv1)) < Test::Constants::test_tolerance_geometry))
        BEMBEL_TEST_IF((
            std::abs((fe2.evaluateOnPatch(0, e2n) - fe2.evaluateOnPatch(2, e2p))
                         .dot(uv1)) < Test::Constants::test_tolerance_geometry))
        BEMBEL_TEST_IF((
            std::abs((fe2.evaluateOnPatch(0, e3n) - fe2.evaluateOnPatch(3, e3p))
                         .dot(uv2)) < Test::Constants::test_tolerance_geometry))
        BEMBEL_TEST_IF((
            std::abs((fe2.evaluateOnPatch(0, e1n) - fe2.evaluateOnPatch(4, e1p))
                         .dot(uv2)) < Test::Constants::test_tolerance_geometry))
        //------------------
        BEMBEL_TEST_IF((
            std::abs((fe3.evaluateOnPatch(0, e1p) - fe3.evaluateOnPatch(1, e4p))
                         .dot(uv1)) < Test::Constants::test_tolerance_geometry))
        BEMBEL_TEST_IF((
            std::abs((fe3.evaluateOnPatch(0, e3p) - fe3.evaluateOnPatch(2, e2p))
                         .dot(uv1)) < Test::Constants::test_tolerance_geometry))
        BEMBEL_TEST_IF((
            std::abs((fe3.evaluateOnPatch(0, e4n) - fe3.evaluateOnPatch(3, e3p))
                         .dot(uv2)) < Test::Constants::test_tolerance_geometry))
        BEMBEL_TEST_IF((
            std::abs((fe3.evaluateOnPatch(0, e2n) - fe3.evaluateOnPatch(4, e1p))
                         .dot(uv2)) < Test::Constants::test_tolerance_geometry))
        //------------------
        BEMBEL_TEST_IF((
            std::abs((fe4.evaluateOnPatch(0, e3n) - fe4.evaluateOnPatch(1, e4p))
                         .dot(uv1)) < Test::Constants::test_tolerance_geometry))
        BEMBEL_TEST_IF((
            std::abs((fe4.evaluateOnPatch(0, e1n) - fe4.evaluateOnPatch(2, e2p))
                         .dot(uv1)) < Test::Constants::test_tolerance_geometry))
        BEMBEL_TEST_IF((
            std::abs((fe4.evaluateOnPatch(0, e2p) - fe4.evaluateOnPatch(3, e3p))
                         .dot(uv2)) < Test::Constants::test_tolerance_geometry))
        BEMBEL_TEST_IF((
            std::abs((fe4.evaluateOnPatch(0, e4p) - fe4.evaluateOnPatch(4, e1p))
                         .dot(uv2)) < Test::Constants::test_tolerance_geometry))
      }
    }
  }
  return 0;
}
