// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information. This file is part of Bembel, the higher order C++ boundary
// element library. It was written as part of a cooperation of J. Doelz, H.
// Harbrecht, S. Kurz, M. Multerer, S. Schoeps, and F. Wolf at Technische
// Universtaet Darmstadt, Universitaet Basel, and Universita della Svizzera
// italiana, Lugano. This source code is subject to the GNU General Public
// License version 3 and provided WITHOUT ANY WARRANTY, see
// <http://www.bembel.eu> for further information.
#include "spltest.h"

// This test assumes correctness of the bernstein polynomials and thus should be
// the second test. It checks weather the deBoor algorithm generates Bernstein
// Polynomials.

int test_deBoor() {
  using namespace Spl;
  constexpr double TOL = 0.000000001;
  constexpr int testpts = 1000;
  constexpr double h = 1. / testpts;

  std::vector<double> xs;

  for (int i = 0; i < testpts + 1; i++) xs.push_back(i * h);

  {
    constexpr int N = 1;
    std::vector<double> ctrl(N, 0);

    std::vector<double> knt;
    for (int i = 0; i < N; i++) knt.push_back(0);
    for (int i = 0; i < N; i++) knt.push_back(1);
    for (auto x : xs) {
      for (int i = 0; i < N; i++) {
        ctrl[i] = 0;
        // std::cout << std::abs(deBoor<double>(ctrl, knt, {x})[0] -
        //                 evalBrnstn<double, N-1>(ctrl.data(), x)) << "\n";
        if (not(std::abs(Spl::deBoor<double>(ctrl, knt, {x})[0] -
                         evalBrnstn<double, N - 1>(ctrl.data(), x)) < TOL))
          return 1;
        ctrl[i] = 1;
      }
    }
  }

  {
    constexpr int N = 2;
    std::vector<double> ctrl(N, 0);

    std::vector<double> knt;
    for (int i = 0; i < N; i++) knt.push_back(0);
    for (int i = 0; i < N; i++) knt.push_back(1);
    for (auto x : xs) {
      for (int i = 0; i < N; i++) {
        ctrl[i] = 0;
        // std::cout << std::abs(deBoor<double>(ctrl, knt, {x})[0] -
        //                 evalBrnstn<double, N-1>(ctrl.data(), x)) << "\n";
        if (not(std::abs(Spl::deBoor<double>(ctrl, knt, {x})[0] -
                         evalBrnstn<double, N - 1>(ctrl.data(), x)) < TOL))
          return 1;
        ctrl[i] = 1;
      }
    }
  }

  {
    constexpr int N = 6;
    std::vector<double> ctrl(N, 0);

    std::vector<double> knt;
    for (int i = 0; i < N; i++) knt.push_back(0);
    for (int i = 0; i < N; i++) knt.push_back(1);
    for (auto x : xs) {
      for (int i = 0; i < N; i++) {
        ctrl[i] = 0;
        // std::cout << std::abs(deBoor<double>(ctrl, knt, {x})[0] -
        //                 evalBrnstn<double, N-1>(ctrl.data(), x)) << "\n";
        if ((not(std::abs(Spl::deBoor<double>(ctrl, knt, {x})[0] -
                          evalBrnstn<double, N - 1>(ctrl.data(), x)) < TOL)))
          return 1;
        ctrl[i] = 1;
      }
    }
  }

  {
    constexpr int N = 9;
    std::vector<double> ctrl(N, 0);

    std::vector<double> knt;
    for (int i = 0; i < N; i++) knt.push_back(0);
    for (int i = 0; i < N; i++) knt.push_back(1);
    for (auto x : xs) {
      for (int i = 0; i < N; i++) {
        ctrl[i] = 0;
        // std::cout << std::abs(deBoor<double>(ctrl, knt, {x})[0] -
        //                 evalBrnstn<double, N-1>(ctrl.data(), x)) << "\n";
        if (not(std::abs(Spl::deBoor<double>(ctrl, knt, {x})[0] -
                         evalBrnstn<double, N - 1>(ctrl.data(), x)) < TOL))
          return 1;
        ctrl[i] = 1;
      }
    }
  }

  {
    constexpr int N = 19;
    std::vector<double> ctrl(N, 0);

    std::vector<double> knt;
    for (int i = 0; i < N; i++) knt.push_back(0);
    for (int i = 0; i < N; i++) knt.push_back(1);
    for (auto x : xs) {
      for (int i = 0; i < N; i++) {
        ctrl[i] = 0;
        // std::cout << std::abs(deBoor<double>(ctrl, knt, {x})[0] -
        //                 evalBrnstn<double, N-1>(ctrl.data(), x)) << "\n";
        if (not(std::abs(Spl::deBoor<double>(ctrl, knt, {x})[0] -
                         evalBrnstn<double, N - 1>(ctrl.data(), x)) < TOL))
          return 1;
        ctrl[i] = 1;
      }
    }
  }

  return 0;
}
