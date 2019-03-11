// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#include "spltest.h"

namespace {
template <int deg>
int testDegree(std::vector<double> xs) {
  using namespace Spl;
  constexpr double h = 0.0000001;
  for (auto x : xs) {
    double vals1[deg + 1];
    double vals2[deg + 1];
    double vals3[deg + 1];
    double vals4[deg + 1];
    double vals5[deg + 1];
    double vals6[deg + 1];
    double coefs[deg + 1];

    for (int i = 0; i < deg + 1; i++) {
      coefs[i] = 0;
      vals1[i] = 0;
      vals2[i] = 0;
      vals3[i] = 0;
      vals4[i] = 0;
      vals5[i] = 0;
      vals6[i] = 0;
    }

    for (int i = 0; i < deg + 1; i++) {
      coefs[i] = 1;
      vals1[i] = evalBrnstnDer<double, deg>(coefs, x);
      vals2[i] = (evalBrnstn<double, deg>(coefs, x + h) -
                  evalBrnstn<double, deg>(coefs, x)) /
                 h;
      vals5[i] = evalBrnstn<double, deg>(coefs, x);
      coefs[i] = 0;
    }

    evalBrnstnDerBasis<double, deg>(vals3, x);

    evalBrnstnBasis<double, deg>(vals4, x + h);
    evalBrnstnBasis<double, deg>(vals6, x);

    for (int i = 0; i < deg + 1; i++) {
      vals4[i] = (vals4[i] - vals6[i]) / h;
    }

    for (int i = 0; i < deg + 1; i++) {
      if (std::abs(vals1[i] - vals2[i]) < 0.0001) {
      } else {
        return 1;
      };
      if (std::abs(vals2[i] - vals3[i]) < 0.0001) {
      } else {
        return 1;
      };
      if (std::abs(vals3[i] - vals4[i]) < 0.0001) {
      } else {
        return 1;
      };
      if (std::abs(vals4[i] - vals1[i]) < 0.0001) {
      } else {
        return 1;
      };
      if (std::abs(vals5[i] - vals6[i]) < 0.0000000001) {
      } else {
        return 1;
      };
    }
  }
  return 0;
}
}  // namespace

int test_bernstein() {
  using namespace Spl;

  constexpr double h = 0.0000001;
  constexpr int numxs = 2500;

  std::vector<double> xs;
  xs.reserve(numxs);
  std::default_random_engine generator(time(0));
  std::uniform_real_distribution<double> distribution(0.0, 1.0);

  for (int i = 0; i < numxs; i++) xs.push_back(distribution(generator));

  // First some analytical tests

  for (auto x : xs) {
    // Values
    double array[5];
    double testarray[5];

    for (int i = 0; i < 5; i++) testarray[i] = 0;

    evalBrnstnBasis<double, 2>(testarray, x);

    array[0] = (1. - x) * (1. - x);
    array[1] = (2. * x) * (1. - x);
    array[2] = x * x;

    for (int i = 0; i < 3; i++) {
      if (std::abs(array[i] - testarray[i]) < 0.00000001 &&
          "test_bernstein.cpp encountert an inaccuracy. There must be "
          "something wrong with the evaluation of the bernstein polynomials.") {
      } else {
        return 1;
      };
    }

    for (int i = 0; i < 5; i++) testarray[i] = 0;

    evalBrnstnBasis<double, 4>(testarray, x);

    array[0] = (1. - x) * (1. - x) * (1. - x) * (1. - x);
    array[1] = 4. * x * (1. - x) * (1. - x) * (1. - x);
    array[2] = 6. * x * x * (1. - x) * (1. - x);
    array[3] = 4. * x * x * x * (1. - x);
    array[4] = x * x * x * x;

    for (int i = 0; i < 5; i++) {
      if (std::abs(array[i] - testarray[i]) < 0.00000001 &&
          "test_bernstein.cpp encountert an inaccuracy. There must be "
          "something wrong with the evaluation of the bernstein polynomials.") {
      } else {
        return 1;
      };
    }
    // Derivatives

    for (int i = 0; i < 5; i++) testarray[i] = 0;

    evalBrnstnDerBasis<double, 2>(testarray, x);

    array[0] = 2. * x - 2;
    array[1] = 2. - 4. * x;
    array[2] = 2. * x;

    for (int i = 0; i < 3; i++) {
      if (std::abs(array[i] - testarray[i]) < 0.00000001 &&
          "test_bernstein.cpp encountert an inaccuracy. There must be "
          "something wrong with the evaluation of the bernstein polynomials.") {
      } else {
        return 1;
      };
    }

    for (int i = 0; i < 5; i++) testarray[i] = 0;

    evalBrnstnDerBasis<double, 4>(testarray, x);

    array[0] = -4 * (1. - x) * (1. - x) * (1. - x);
    array[1] = -16. * (1. - x) * (1. - x) * (-0.25 + x);
    array[2] = x * (12. - 36. * x + 24. * x * x);
    array[3] = (12 - 16 * x) * x * x;
    array[4] = 4. * x * x * x;

    for (int i = 0; i < 5; i++) {
      if (std::abs(array[i] - testarray[i]) < 0.00000001 &&
          "test_bernstein.cpp encountert an inaccuracy. There must be "
          "something wrong with the evaluation of the bernstein polynomials.") {
      } else {
        return 1;
      };
    }
  }

  // Now comparisons of different methods that should yield the same result.
  // First, hardcoded templates
  int testint = 0;
  testint += testDegree<0>(xs);
  testint += testDegree<1>(xs);
  testint += testDegree<2>(xs);
  testint += testDegree<3>(xs);
  testint += testDegree<4>(xs);
  testint += testDegree<5>(xs);
  testint += testDegree<6>(xs);
  testint += testDegree<7>(xs);
  testint += testDegree<8>(xs);
  testint += testDegree<9>(xs);
  testint += testDegree<10>(xs);
  testint += testDegree<11>(xs);
  testint += testDegree<12>(xs);
  testint += testDegree<13>(xs);
  testint += testDegree<14>(xs);
  testint += testDegree<15>(xs);
  testint += testDegree<16>(xs);
  testint += testDegree<17>(xs);
  testint += testDegree<18>(xs);
  if (testint > 0) return 1;

  // Now, integer lookup

  for (int deg = 0; deg < 19; deg++) {
    constexpr double h = 0.0000001;
    for (auto x : xs) {
      double vals1[deg + 1];
      double vals2[deg + 1];
      double vals3[deg + 1];
      double vals4[deg + 1];
      double vals5[deg + 1];
      double vals6[deg + 1];
      double coefs[deg + 1];

      for (int i = 0; i < deg + 1; i++) {
        coefs[i] = 0;
        vals1[i] = 0;
        vals2[i] = 0;
        vals3[i] = 0;
        vals4[i] = 0;
        vals5[i] = 0;
        vals6[i] = 0;
      }

      for (int i = 0; i < deg + 1; i++) {
        coefs[i] = 1;
        vals1[i] = evalBrnstnDer(deg, coefs, x);
        vals2[i] =
            (evalBrnstn(deg, coefs, x + h) - evalBrnstn(deg, coefs, x)) / h;

        vals5[i] = evalBrnstn(deg, coefs, x);
        coefs[i] = 0;
      }

      evalBrnstnDerBasis(deg, vals3, x);

      evalBrnstnBasis(deg, vals4, x + h);
      evalBrnstnBasis(deg, vals6, x);

      for (int i = 0; i < deg + 1; i++) {
        vals4[i] = (vals4[i] - vals6[i]) / h;
      }

      for (int i = 0; i < deg + 1; i++) {
        if (std::abs(vals1[i] - vals2[i]) < 0.0001) {
        } else {
          return 1;
        };
        if (std::abs(vals2[i] - vals3[i]) < 0.0001) {
        } else {
          return 1;
        };
        if (std::abs(vals3[i] - vals4[i]) < 0.0001) {
        } else {
          return 1;
        };
        if (std::abs(vals4[i] - vals1[i]) < 0.0001) {
        } else {
          return 1;
        };
        if (std::abs(vals5[i] - vals6[i]) < 0.0000000001) {
        } else {
          return 1;
        };
      }
    }
  }

  // Now, we check if PEQ and EQ, i.e., the basis lookup with = and += work as
  // expected.
  for (auto x : xs) {
    for (int deg = 0; deg < 19; deg++) {
      double* eq;
      double* peq;
      eq = (double*)malloc((deg + 1) * sizeof(double));
      peq = (double*)malloc((deg + 1) * sizeof(double));

      for (int i = 0; i < deg + 1; i++) {
        eq[i] = 1.3;
        peq[i] = 1.3;
      }

      evalBrnstnBasis(deg, eq, x);
      evalBrnstnBasisPEQ(deg, peq, x);

      for (int i = 0; i < deg + 1; i++) {
        if (std::abs(eq[i] - (peq[i] - 1.3)) < 0.0000000001) {
        } else {
          return 1;
        }
      }
      free(eq);
      free(peq);
    }
  }

  return 0;
}