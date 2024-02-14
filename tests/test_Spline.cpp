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

#include <Bembel/Spline>

#include "tests/Test.hpp"

int main() {
  using namespace Bembel;

  // We test the Bernstein Basis of the BasisHandler against the deBoor code
  for (int p = 0; p <= Bembel::Constants::MaxP; ++p) {
    for (auto x : Test::Constants::eq_points) {
      Eigen::VectorXd result1 = Eigen::VectorXd::Zero(p + 1);
      Eigen::VectorXcd result1_complex = Eigen::VectorXcd::Zero(p + 1);
      Basis::BasisHandler<double>::phi(p, &result1, 1, x);
      Basis::BasisHandler<std::complex<double>>::phi(p, &result1_complex, 1, x);

      Eigen::VectorXd result2 = Eigen::VectorXd::Zero(p + 1);
      Eigen::MatrixXd coef = Eigen::VectorXd::Zero(p + 1).transpose();
      for (int i = 0; i < p + 1; ++i) {
        coef(i) = 1;
        std::vector<double> v = {x};
        result2(i) = Spl::DeBoor(coef, Spl::MakeBezierKnotVector(p + 1), v)(0);
        coef(i) = 0;
      }
      BEMBEL_TEST_IF((result1 - result2).norm() <
                     Test::Constants::coefficient_accuracy);
      BEMBEL_TEST_IF((result1_complex - result2).norm() <
                     Test::Constants::coefficient_accuracy);
    }
  }

  // Now, we do the same for the derivatives
  for (int p = 1; p <= Bembel::Constants::MaxP; ++p) {
    for (auto x : Test::Constants::eq_points) {
      Eigen::VectorXd result1 = Eigen::VectorXd::Zero(p + 1);
      Eigen::VectorXcd result1_complex = Eigen::VectorXcd::Zero(p + 1);
      Basis::BasisHandler<double>::phiDx(p, &result1, 1, x);
      Basis::BasisHandler<std::complex<double>>::phiDx(p, &result1_complex, 1,
                                                       x);

      Eigen::VectorXd result2 = Eigen::VectorXd::Zero(p + 1);
      Eigen::MatrixXd coef = Eigen::VectorXd::Zero(p + 1).transpose();
      for (int i = 0; i < p + 1; ++i) {
        coef(i) = 1;
        std::vector<double> v = {x};
        result2(i) =
            Spl::DeBoorDer(coef, Spl::MakeBezierKnotVector(p + 1), v)(0);
        coef(i) = 0;
      }
      BEMBEL_TEST_IF((result1 - result2).norm() <
                     Test::Constants::coefficient_accuracy);
      BEMBEL_TEST_IF((result1_complex - result2).norm() <
                     Test::Constants::coefficient_accuracy);
    }
  }

  // We test the tensor product Bernstein Basis of the BasisHandler against the
  // deBoor code
  for (int p = 0; p <= Bembel::Constants::MaxP; ++p) {
    for (auto x : Test::Constants::eq_points)
      for (auto y : Test::Constants::eq_points) {
        Eigen::VectorXd result1 = Eigen::VectorXd::Zero((p + 1) * (p + 1));
        Eigen::VectorXcd result1_complex =
            Eigen::VectorXcd::Zero((p + 1) * (p + 1));
        Basis::BasisHandler<double>::phiPhi(p, &result1, 1,
                                            Eigen::Vector2d(x, y));
        Basis::BasisHandler<std::complex<double>>::phiPhi(
            p, &result1_complex, 1, Eigen::Vector2d(x, y));

        Eigen::VectorXd result2_x = Eigen::VectorXd::Zero(p + 1);
        Eigen::MatrixXd coef_x = Eigen::VectorXd::Zero(p + 1).transpose();
        Eigen::VectorXd result2_y = Eigen::VectorXd::Zero(p + 1);
        Eigen::MatrixXd coef_y = Eigen::VectorXd::Zero(p + 1).transpose();
        for (int i = 0; i < p + 1; ++i) {
          coef_x(i) = 1;
          std::vector<double> v = {x};
          result2_x(i) =
              Spl::DeBoor(coef_x, Spl::MakeBezierKnotVector(p + 1), v)(0);
          coef_x(i) = 0;
        }
        for (int i = 0; i < p + 1; ++i) {
          coef_y(i) = 1;
          std::vector<double> v = {y};
          result2_y(i) =
              Spl::DeBoor(coef_y, Spl::MakeBezierKnotVector(p + 1), v)(0);
          coef_y(i) = 0;
        }

        Eigen::VectorXd result2 = Eigen::VectorXd::Zero((p + 1) * (p + 1));
        for (int iy = 0; iy < p + 1; ++iy)
          for (int ix = 0; ix < p + 1; ++ix)
            result2(iy * (p + 1) + ix) = result2_x(ix) * result2_y(iy);

        BEMBEL_TEST_IF((result1 - result2).norm() <
                       Test::Constants::coefficient_accuracy);
        BEMBEL_TEST_IF((result1_complex - result2).norm() <
                       Test::Constants::coefficient_accuracy);
      }
  }

  const double x1 = Test::Constants::eq_points[3];
  const double y1 = Test::Constants::eq_points[7];
  const double x2 = Test::Constants::eq_points[2];
  const double y2 = Test::Constants::eq_points[8];
  // We test the interactions of two phiphis at exemplary points
  for (int p = 0; p <= Bembel::Constants::MaxP; ++p) {
    Eigen::MatrixXd result1 =
        Eigen::MatrixXd::Zero((p + 1) * (p + 1), (p + 1) * (p + 1));

    Eigen::MatrixXcd result1_complex =
        Eigen::MatrixXcd::Zero((p + 1) * (p + 1), (p + 1) * (p + 1));
    Basis::BasisHandler<double>::phiTimesPhi(
        p, &result1, 1, Eigen::Vector2d(x1, y1), Eigen::Vector2d(x2, y2));
    Basis::BasisHandler<std::complex<double>>::phiTimesPhi(
        p, &result1_complex, 1, Eigen::Vector2d(x1, y1),
        Eigen::Vector2d(x2, y2));

    Eigen::VectorXd result2_x1 = Eigen::VectorXd::Zero(p + 1);
    Eigen::MatrixXd coef_x1 = Eigen::VectorXd::Zero(p + 1).transpose();
    Eigen::VectorXd result2_y1 = Eigen::VectorXd::Zero(p + 1);
    Eigen::MatrixXd coef_y1 = Eigen::VectorXd::Zero(p + 1).transpose();
    for (int i = 0; i < p + 1; ++i) {
      coef_x1(i) = 1;
      std::vector<double> v = {x1};
      result2_x1(i) =
          Spl::DeBoor(coef_x1, Spl::MakeBezierKnotVector(p + 1), v)(0);
      coef_x1(i) = 0;
    }
    for (int i = 0; i < p + 1; ++i) {
      coef_y1(i) = 1;
      std::vector<double> v = {y1};
      result2_y1(i) =
          Spl::DeBoor(coef_y1, Spl::MakeBezierKnotVector(p + 1), v)(0);
      coef_y1(i) = 0;
    }
    Eigen::VectorXd result2_x2 = Eigen::VectorXd::Zero(p + 1);
    Eigen::MatrixXd coef_x2 = Eigen::VectorXd::Zero(p + 1).transpose();
    Eigen::VectorXd result2_y2 = Eigen::VectorXd::Zero(p + 1);
    Eigen::MatrixXd coef_y2 = Eigen::VectorXd::Zero(p + 1).transpose();
    for (int i = 0; i < p + 1; ++i) {
      coef_x2(i) = 1;
      std::vector<double> v = {x2};
      result2_x2(i) =
          Spl::DeBoor(coef_x2, Spl::MakeBezierKnotVector(p + 1), v)(0);
      coef_x2(i) = 0;
    }
    for (int i = 0; i < p + 1; ++i) {
      coef_y2(i) = 1;
      std::vector<double> v = {y2};
      result2_y2(i) =
          Spl::DeBoor(coef_y2, Spl::MakeBezierKnotVector(p + 1), v)(0);
      coef_y2(i) = 0;
    }

    Eigen::VectorXd result2_1 = Eigen::VectorXd::Zero((p + 1) * (p + 1));
    for (int iy = 0; iy < p + 1; ++iy)
      for (int ix = 0; ix < p + 1; ++ix)
        result2_1(iy * (p + 1) + ix) = result2_x1(ix) * result2_y1(iy);

    Eigen::VectorXd result2_2 = Eigen::VectorXd::Zero((p + 1) * (p + 1));
    for (int iy = 0; iy < p + 1; ++iy)
      for (int ix = 0; ix < p + 1; ++ix)
        result2_2(iy * (p + 1) + ix) = result2_x2(ix) * result2_y2(iy);

    BEMBEL_TEST_IF((result1 - result2_1 * result2_2.transpose()).norm() <
                   Test::Constants::coefficient_accuracy);
    BEMBEL_TEST_IF(
        (result1_complex - result2_1 * result2_2.transpose()).norm() <
        Test::Constants::coefficient_accuracy);
  }

  // We test the interactions of the div of two phiphis at exemplary points
  for (int p = 1; p <= Bembel::Constants::MaxP; ++p) {
    Eigen::MatrixXd result1 =
        Eigen::MatrixXd::Zero(2 * (p + 1) * (p + 1), 2 * (p + 1) * (p + 1));

    Eigen::MatrixXcd result1_complex =
        Eigen::MatrixXcd::Zero(2 * (p + 1) * (p + 1), 2 * (p + 1) * (p + 1));
    Basis::BasisHandler<double>::divPhiTimesDivPhi(
        p, &result1, 1, Eigen::Vector2d(x1, y1), Eigen::Vector2d(x2, y2));
    Basis::BasisHandler<std::complex<double>>::divPhiTimesDivPhi(
        p, &result1_complex, 1, Eigen::Vector2d(x1, y1),
        Eigen::Vector2d(x2, y2));

    Eigen::VectorXd result2_x1_der = Eigen::VectorXd::Zero(p + 1);
    Eigen::MatrixXd coef_x1_der = Eigen::VectorXd::Zero(p + 1).transpose();
    Eigen::VectorXd result2_y1_der = Eigen::VectorXd::Zero(p + 1);
    Eigen::MatrixXd coef_y1_der = Eigen::VectorXd::Zero(p + 1).transpose();
    for (int i = 0; i < p + 1; ++i) {
      coef_x1_der(i) = 1;
      std::vector<double> v = {x1};
      result2_x1_der(i) =
          Spl::DeBoorDer(coef_x1_der, Spl::MakeBezierKnotVector(p + 1), v)(0);
      coef_x1_der(i) = 0;
    }
    for (int i = 0; i < p + 1; ++i) {
      coef_y1_der(i) = 1;
      std::vector<double> v = {y1};
      result2_y1_der(i) =
          Spl::DeBoorDer(coef_y1_der, Spl::MakeBezierKnotVector(p + 1), v)(0);
      coef_y1_der(i) = 0;
    }
    Eigen::VectorXd result2_x1 = Eigen::VectorXd::Zero(p + 1);
    Eigen::MatrixXd coef_x1 = Eigen::VectorXd::Zero(p + 1).transpose();
    Eigen::VectorXd result2_y1 = Eigen::VectorXd::Zero(p + 1);
    Eigen::MatrixXd coef_y1 = Eigen::VectorXd::Zero(p + 1).transpose();
    for (int i = 0; i < p + 1; ++i) {
      coef_x1(i) = 1;
      std::vector<double> v = {x1};
      result2_x1(i) =
          Spl::DeBoor(coef_x1, Spl::MakeBezierKnotVector(p + 1), v)(0);
      coef_x1(i) = 0;
    }
    for (int i = 0; i < p + 1; ++i) {
      coef_y1(i) = 1;
      std::vector<double> v = {y1};
      result2_y1(i) =
          Spl::DeBoor(coef_y1, Spl::MakeBezierKnotVector(p + 1), v)(0);
      coef_y1(i) = 0;
    }

    Eigen::VectorXd result2_x2_der = Eigen::VectorXd::Zero(p + 1);
    Eigen::MatrixXd coef_x2_der = Eigen::VectorXd::Zero(p + 1).transpose();
    Eigen::VectorXd result2_y2_der = Eigen::VectorXd::Zero(p + 1);
    Eigen::MatrixXd coef_y2_der = Eigen::VectorXd::Zero(p + 1).transpose();
    for (int i = 0; i < p + 1; ++i) {
      coef_x2_der(i) = 1;
      std::vector<double> v = {x2};
      result2_x2_der(i) =
          Spl::DeBoorDer(coef_x2_der, Spl::MakeBezierKnotVector(p + 1), v)(0);
      coef_x2_der(i) = 0;
    }
    for (int i = 0; i < p + 1; ++i) {
      coef_y2_der(i) = 1;
      std::vector<double> v = {y2};
      result2_y2_der(i) =
          Spl::DeBoorDer(coef_y2_der, Spl::MakeBezierKnotVector(p + 1), v)(0);
      coef_y2_der(i) = 0;
    }
    Eigen::VectorXd result2_x2 = Eigen::VectorXd::Zero(p + 1);
    Eigen::MatrixXd coef_x2 = Eigen::VectorXd::Zero(p + 1).transpose();
    Eigen::VectorXd result2_y2 = Eigen::VectorXd::Zero(p + 1);
    Eigen::MatrixXd coef_y2 = Eigen::VectorXd::Zero(p + 1).transpose();
    for (int i = 0; i < p + 1; ++i) {
      coef_x2(i) = 1;
      std::vector<double> v = {x2};
      result2_x2(i) =
          Spl::DeBoor(coef_x2, Spl::MakeBezierKnotVector(p + 1), v)(0);
      coef_x2(i) = 0;
    }
    for (int i = 0; i < p + 1; ++i) {
      coef_y2(i) = 1;
      std::vector<double> v = {y2};
      result2_y2(i) =
          Spl::DeBoor(coef_y2, Spl::MakeBezierKnotVector(p + 1), v)(0);
      coef_y2(i) = 0;
    }

    Eigen::VectorXd result2_dx1 = Eigen::VectorXd::Zero((p + 1) * (p + 1));
    for (int iy = 0; iy < p + 1; ++iy)
      for (int ix = 0; ix < p + 1; ++ix)
        result2_dx1(iy * (p + 1) + ix) = result2_x1_der(ix) * result2_y1(iy);

    Eigen::VectorXd result2_dy1 = Eigen::VectorXd::Zero((p + 1) * (p + 1));
    for (int iy = 0; iy < p + 1; ++iy)
      for (int ix = 0; ix < p + 1; ++ix)
        result2_dy1(iy * (p + 1) + ix) = result2_x1(ix) * result2_y1_der(iy);

    Eigen::VectorXd result2_dx2 = Eigen::VectorXd::Zero((p + 1) * (p + 1));
    for (int iy = 0; iy < p + 1; ++iy)
      for (int ix = 0; ix < p + 1; ++ix)
        result2_dx2(iy * (p + 1) + ix) = result2_x2_der(ix) * result2_y2(iy);

    Eigen::VectorXd result2_dy2 = Eigen::VectorXd::Zero((p + 1) * (p + 1));
    for (int iy = 0; iy < p + 1; ++iy)
      for (int ix = 0; ix < p + 1; ++ix)
        result2_dy2(iy * (p + 1) + ix) = result2_x2(ix) * result2_y2_der(iy);

    Eigen::MatrixXd result2 =
        Eigen::MatrixXd::Zero(2 * (p + 1) * (p + 1), 2 * (p + 1) * (p + 1));

    const int pp1s = (p + 1) * (p + 1);
    result2.block(0, 0, pp1s, pp1s) = result2_dx1 * result2_dx2.transpose();
    result2.block(0, pp1s, pp1s, pp1s) = result2_dx1 * result2_dy2.transpose();
    result2.block(pp1s, 0, pp1s, pp1s) = result2_dy1 * result2_dx2.transpose();
    result2.block(pp1s, pp1s, pp1s, pp1s) =
        result2_dy1 * result2_dy2.transpose();

    BEMBEL_TEST_IF((result1 - result2).norm() <
                   Test::Constants::coefficient_accuracy);
    BEMBEL_TEST_IF((result1_complex - result2).norm() <
                   Test::Constants::coefficient_accuracy);
  }

  return 0;
}
