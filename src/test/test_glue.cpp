// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#include "bemtest.h"

using namespace Bembel;

std::vector<Spl::Patch> make_doublescreen(int pos_2ndpatch, int turns) {
  Spl::Patch center, second;

  std::vector<double> knot = {0.0, 0.0, 1.0, 1.0};
  Eigen::MatrixXd X_c(2, 2), X(2, 2), Y_c(2, 2), Y(2, 2), Z(2, 2), ws(2, 2);
  Z << 0, 0, 0, 0;
  X << 0, 0, 1, 1;
  Y << 0, 1, 0, 1;
  ws << 1, 1, 1, 1;

  double x, y;

  switch (pos_2ndpatch) {
    case (0): {
      x = 0;
      y = -1;
    } break;
    case (1): {
      x = 1;
      y = 0;
    } break;
    case (2): {
      x = 0;
      y = 1;
    } break;
    default: {
      assert(pos_2ndpatch == 3);
      x = -1;
      y = 0;
    }
  }

  X_c << x, x, x + 1, x + 1;
  Y_c << y, y + 1, y, y + 1;

  switch (turns) {
    case (0):
      X << 0, 0, 1, 1;
      Y << 0, 1, 0, 1;
    case (1): {
      X << 0, 1, 0, 1;
      Y << 1, 1, 0, 0;
    } break;
    case (2): {
      X << 1, 1, 0, 0;
      Y << 1, 0, 1, 0;
    } break;
    default: {
      assert(turns == 3);
      X << 1, 0, 1, 0;
      Y << 0, 0, 1, 1;
    }
  }

  center.initHom({X_c, Y_c, Z, ws}, knot, knot);
  second.initHom({X, Y, Z, ws}, knot, knot);

  return {center, second};
}
/**
 *  \author        Felix Wolf
 */
std::vector<int> find_local_edge_dofs_laplace_dc(int order, int edge,
                                                 int shift) {
  assert(order == 1 or order == 2);

  assert(edge >= 0 and edge < 4);
  if (order == 1) {
    switch (edge) {
      case (0):
        return {14 + shift, 12 + shift, 2 + shift, 0 + shift};
      case (1):
        return {11 + shift, 10 + shift, 15 + shift, 14 + shift};
      case (2):
        return {5 + shift, 7 + shift, 9 + shift, 11 + shift};
      default:
        return {0 + shift, 1 + shift, 4 + shift, 5 + shift};
    }
  } else {
    switch (edge) {
      case (0):
        return {33 + shift, 30 + shift, 27 + shift,
                6 + shift,  3 + shift,  0 + shift};
      case (1):
        return {26 + shift, 25 + shift, 24 + shift,
                35 + shift, 34 + shift, 33 + shift};
      case (2):
        return {11 + shift, 14 + shift, 17 + shift,
                20 + shift, 23 + shift, 26 + shift};
      default:
        return {0 + shift, 1 + shift,  2 + shift,
                9 + shift, 10 + shift, 11 + shift};
    }
  }
}

struct edge_dofs {
  std::array<std::vector<int>, 4> FRST_X;
  std::array<std::vector<int>, 4> FRST_Y;
  std::array<std::vector<int>, 4> SCND_X;
  std::array<std::vector<int>, 4> SCND_Y;
};
/**
 *  \author        Felix Wolf
 */
inline edge_dofs get_index_handwritten(int p, int M) {
  assert((M == 1 or M == 2) and (p == 1 or p == 2));

  edge_dofs ed;

  // p = 1; M = 1
  if (p == 1 and M == 1) {
    ed.FRST_X[0] = {3, 0};
    ed.FRST_X[1] = {5, 4, 3};
    ed.FRST_X[2] = {2, 5};
    ed.FRST_X[3] = {0, 1, 2};

    ed.SCND_X[0] = {9, 6};
    ed.SCND_X[1] = {11, 10, 9};
    ed.SCND_X[2] = {8, 11};
    ed.SCND_X[3] = {6, 7, 8};

    ed.FRST_Y[0] = {16, 14, 12};
    ed.FRST_Y[1] = {17, 16};
    ed.FRST_Y[2] = {13, 15, 17};
    ed.FRST_Y[3] = {12, 13};

    ed.SCND_Y[0] = {22, 20, 18};
    ed.SCND_Y[1] = {23, 22};
    ed.SCND_Y[2] = {19, 21, 23};
    ed.SCND_Y[3] = {18, 19};
  }

  // p = 2; M = 1
  if (p == 2 and M == 1) {
    ed.FRST_X[0] = {8, 4, 0};
    ed.FRST_X[1] = {11, 10, 9, 8};
    ed.FRST_X[2] = {3, 7, 11};
    ed.FRST_X[3] = {0, 1, 2, 3};

    ed.SCND_X[0] = {20, 16, 12};
    ed.SCND_X[1] = {23, 22, 21, 20};
    ed.SCND_X[2] = {15, 19, 23};
    ed.SCND_X[3] = {12, 13, 14, 15};

    ed.FRST_Y[0] = {33, 30, 27, 24};
    ed.FRST_Y[1] = {35, 34, 33};
    ed.FRST_Y[2] = {26, 29, 32, 35};
    ed.FRST_Y[3] = {24, 25, 26};

    ed.SCND_Y[0] = {45, 42, 39, 36};
    ed.SCND_Y[1] = {47, 46, 45};
    ed.SCND_Y[2] = {38, 41, 44, 47};
    ed.SCND_Y[3] = {36, 37, 38};
  }

  // p = 1; M = 2
  if (p == 1 and M == 2) {
    ed.FRST_X[0] = {15, 10, 5, 0};
    ed.FRST_X[1] = {19, 18, 17, 16, 15};
    ed.FRST_X[2] = {4, 9, 14, 19};
    ed.FRST_X[3] = {0, 1, 2, 3, 4};

    ed.SCND_X[0] = {35, 30, 25, 20};
    ed.SCND_X[1] = {39, 38, 37, 36, 35};
    ed.SCND_X[2] = {24, 29, 34, 39};
    ed.SCND_X[3] = {20, 21, 22, 23, 24};

    ed.FRST_Y[0] = {56, 52, 48, 44, 40};
    ed.FRST_Y[1] = {59, 58, 57, 56};
    ed.FRST_Y[2] = {43, 47, 51, 55, 59};
    ed.FRST_Y[3] = {40, 41, 42, 43};

    ed.SCND_Y[0] = {76, 72, 68, 64, 60};
    ed.SCND_Y[1] = {79, 78, 77, 76};
    ed.SCND_Y[2] = {63, 67, 71, 75, 79};
    ed.SCND_Y[3] = {60, 61, 62, 63};
  }

  // p = 2; M = 2
  if (p == 2 and M == 2) {
    ed.FRST_X[0] = {24, 18, 12, 6, 0};
    ed.FRST_X[1] = {29, 28, 27, 26, 25, 24};
    ed.FRST_X[2] = {5, 11, 17, 23, 29};
    ed.FRST_X[3] = {0, 1, 2, 3, 4, 5};

    ed.SCND_X[0] = {54, 48, 42, 36, 30};
    ed.SCND_X[1] = {59, 58, 57, 56, 55, 54};
    ed.SCND_X[2] = {35, 41, 47, 53, 59};
    ed.SCND_X[3] = {30, 31, 32, 33, 34, 35};

    ed.FRST_Y[0] = {85, 80, 75, 70, 65, 60};
    ed.FRST_Y[1] = {89, 88, 87, 86, 85};
    ed.FRST_Y[2] = {64, 69, 74, 79, 84, 89};
    ed.FRST_Y[3] = {60, 61, 62, 63, 64};

    ed.SCND_Y[0] = {115, 110, 105, 100, 95, 90};
    ed.SCND_Y[1] = {119, 118, 117, 116, 115};
    ed.SCND_Y[2] = {94, 99, 104, 109, 114, 119};
    ed.SCND_Y[3] = {90, 91, 92, 93, 94};
  }

  return ed;
}
/**
 *  \author        Felix Wolf
 */
int get_complex_shift(int p, int M) {
  if (p == 1 and M == 1) {
    return 24;
  }
  if (p == 2 and M == 1) {
    return 48;
  }
  if (p == 1 and M == 2) {
    return 80;
  }
  if (p == 2 and M == 2) {
    return 120;
  }

  assert(false && "This case is not covered");
  return 0;
}

template <typename T>
inline bool isIn(T x, std::vector<T> xs) {
  for (auto t : xs) {
    if (t == x) {
      return true;
    }
  }
  return false;
}

template <typename T>
inline void sort(bool reverse, std::vector<T> &in) {
  std::sort(in.begin(), in.end());
  if (reverse) {
    const int sz = in.size();
    const std::vector<T> tmp = in;
    for (int i = 0; i < sz; i++) in[i] = tmp[sz - i - 1];
  }
  return;
}
/**
 *  \author        Felix Wolf
 */
bool positionIsAsExpected(int gluecond1, int gluecond2, bool reverse,
                          std::vector<Spl::Patch> geo) {
  /*
matchpn =>
0 -> x = 0 kante
1 -> y = 1 kante
2 -> x = 1 kante
3 -> y = 0 kante
*/

  const int geosize = geo.size();

  assert(geo.size() == 2);

  Eigen::Vector3d p1, p2, p10, p11, p20, p21;

  switch (gluecond1) {
    case (0):
      p1 = geo[0].eval(0, .5);
      p10 = geo[0].eval(0, 0);
      p11 = geo[0].eval(0, 1);
      break;
    case (1):
      p1 = geo[0].eval(.5, 1);
      p10 = geo[0].eval(0, 1);
      p11 = geo[0].eval(1, 1);
      break;
    case (2):
      p1 = geo[0].eval(1, .5);
      p10 = geo[0].eval(1, 0);
      p11 = geo[0].eval(1, 1);
      break;
    default:
      assert(gluecond1 == 3);
      p1 = geo[0].eval(.5, 0);
      p10 = geo[0].eval(0, 0);
      p11 = geo[0].eval(1, 0);
      break;
  }

  switch (gluecond2) {
    case (0):
      p2 = geo[1].eval(0, .5);
      p20 = geo[1].eval(0, 0);
      p21 = geo[1].eval(0, 1);
      break;
    case (1):
      p2 = geo[1].eval(.5, 1);
      p20 = geo[1].eval(0, 1);
      p21 = geo[1].eval(1, 1);
      break;
    case (2):
      p2 = geo[1].eval(1, .5);
      p20 = geo[1].eval(1, 0);
      p21 = geo[1].eval(1, 1);
      break;
    default:
      assert(gluecond2 == 3);
      p2 = geo[1].eval(.5, 0);
      p20 = geo[1].eval(0, 0);
      p21 = geo[1].eval(1, 0);
      break;
  }

  constexpr double tol = 0.00000001;
  bool sameedge = (p1 - p2).norm() < tol;

  bool same0 = reverse ? (p10 - p21).norm() < tol : (p10 - p20).norm() < tol;
  bool same1 = reverse ? (p11 - p20).norm() < tol : (p11 - p21).norm() < tol;

  if (sameedge == false) {
    std::cout << "Wrong edge found" << std::endl;
  }

  bool sameorientation = same0 and same1;
  if (sameorientation == false) {
    std::cout << "Orientation Mismatch" << std::endl;
  }
  return (sameedge and sameorientation);
}
/**
 *	\brief				 This functions tests if the glue-matrix
 *is assembled correctly \author        Felix Wolf
 */
int Test::test_glue() {
  int failcount = 0;

  for (int M = 1; M < 3; M++)
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++)
        for (int p = 1; p < 3; p++) {
          auto geo = make_doublescreen(i, j);

          // Spl::geometry(geo, 5);

          auto gls = boundary_match(geo);
          auto ed = get_index_handwritten(p, M);
          const int complexshift = get_complex_shift(p, M);

          // only one edge should be detected

          assert(gls.size() == 1);

          auto glue = gls[0];
#ifdef _BEMBEL_TEST_WITH_ASSERT_
          assert(glue.patch1 == 0 and glue.patch2 == 1);
#endif
          if (not(glue.patch1 == 0 and glue.patch2 == 1)) {
            failcount = 1;
          }
          // edges, see top for meaning of case
          const int ec1 = glue.matchp1;
          const int ec2 = glue.matchp2;
          const bool reverse = glue.reverse;

          // Check if the direction has been identified correctly
          assert(glue.gluecond == 1 or glue.gluecond == -1);

          switch (ec1) {
            case (0): {
              switch (ec2) {
                case (0):
                  assert(glue.gluecond == -1);
                  break;
                case (1):
                  assert(glue.gluecond == 1);
                  break;
                case (2):
                  assert(glue.gluecond == 1);
                  break;
                default:
                  assert(glue.gluecond == -1);
                  break;
              }
            } break;
            case (1): {
              switch (ec2) {
                case (0):
                  assert(glue.gluecond == 1);
                  break;
                case (1):
                  assert(glue.gluecond == -1);
                  break;
                case (2):
                  assert(glue.gluecond == -1);
                  break;
                default:
                  assert(glue.gluecond == 1);
                  break;
              }
            } break;
            case (2): {
              switch (ec2) {
                case (0):
                  assert(glue.gluecond == 1);
                  break;
                case (1):
                  assert(glue.gluecond == -1);
                  break;
                case (2):
                  assert(glue.gluecond == -1);
                  break;
                default:
                  assert(glue.gluecond == 1);
                  break;
              }
            } break;
            default: {
              switch (ec2) {
                case (0):
                  assert(glue.gluecond == -1);
                  break;
                case (1):
                  assert(glue.gluecond == 1);
                  break;
                case (2):
                  assert(glue.gluecond == 1);
                  break;
                default:
                  assert(glue.gluecond == -1);
                  break;
              }
            } break;
          }

#if 0
          std::cout << "ec1 : " << ec1 << "  ec2 : " << ec2
                    << "   reverse : " << reverse << " p : " << p
                    << "  M : " << M << std::endl;
#endif
#ifdef _BEMBEL_TEST_WITH_ASSERT_
          assert(positionIsAsExpected(ec1, ec2, reverse, geo));
#endif
          if (not positionIsAsExpected(ec1, ec2, reverse, geo)) {
            failcount = 1;
          }

          MaxwellSingle maxwell;
          maxwell._pde.kappa[0] = 1;
          maxwell._pde.kappa[1] = 0;
          Discretization<MaxwellSingle> myDisc;
          Geometry myGeom(geo);

          myDisc.init_Discretization(myGeom, maxwell, p, 1, M);

          discretization disc = myDisc.get_disc();

          // pdeproblem pde = pdeproblem_MaxwellSingle();
          // pde.kappa[0] = 1;
          // pde.kappa[1] = 0;
          // meshdata mesh = get_meshdata(geo, M);
          // discretization disc = get_discretization_NB(p + 1, 1, &pde, &mesh);

          Gluebucket gb(&disc, gls, p, M, 1);

          std::vector<int> from = gb.gluefrom;
          std::vector<int> to = gb.glueto;
#ifdef _BEMBEL_TEST_WITH_ASSERT_
          assert(from.size() == to.size());
#endif
          if (not(from.size() == to.size())) {
            failcount = 1;
          }
          auto vectors = get_index_handwritten(p, M);

          /*
        matchpn =>
        0 -> x = 0 kante
        1 -> y = 1 kante
        2 -> x = 1 kante
        3 -> y = 0 kante
        */

          std::vector<int> myfrom, myto;

          switch (ec1) {
            case (0):
              myfrom = vectors.FRST_X[ec1];
              break;
            case (1):
              myfrom = vectors.FRST_Y[ec1];
              break;
            case (2):
              myfrom = vectors.FRST_X[ec1];
              break;
            default:

              assert(ec1 == 3);
              myfrom = vectors.FRST_Y[ec1];
              break;
          }

          switch (ec2) {
            case (0):
              myto = vectors.SCND_X[ec2];
              break;
            case (1):
              myto = vectors.SCND_Y[ec2];
              break;
            case (2):
              myto = vectors.SCND_X[ec2];
              break;
            default:

              assert(ec2 == 3);
              myto = vectors.SCND_Y[ec2];
              break;
          }

          std::sort(myfrom.begin(), myfrom.end());
          std::sort(myto.begin(), myto.end());
          sort(false, myfrom);
          sort(reverse, myto);

#if 0
          std::cout << "  from : ";
          for (auto x : from) {
            std::cout << x << ", ";
          }
          std::cout << std::endl;

          std::cout << "myfrom : ";
          for (auto x : myfrom) {
            std::cout << x << ", ";
          }
          std::cout << std::endl;
          std::cout << "  to   : ";
          for (auto x : to) {
            std::cout << x << ", ";
          }
          std::cout << std::endl;

          std::cout << "myto   : ";
          for (auto x : myto) {
            std::cout << x << ", ";
          }
          std::cout << std::endl;
#endif

#ifdef _BEMBEL_TEST_WITH_ASSERT_
          assert(myfrom.size() * 2 == from.size());
          assert(myto.size() * 2 == to.size());
          assert(myto.size() == to.size() / 2);

#endif
          if (not(myfrom.size() * 2 == from.size() and
                  myto.size() * 2 == to.size())) {
            failcount = 1;
          }

          for (auto x : from) {
#ifdef _BEMBEL_TEST_WITH_ASSERT_
            assert(isIn(x, myfrom) or isIn(x - complexshift, myfrom));
#endif
            if (not(isIn(x, myfrom) or isIn(x - complexshift, myfrom))) {
              failcount = 1;
            }
          }

          if (not(myto.size() == to.size() / 2)) {
            failcount = 1;
          }
          for (auto x : to) {
#ifdef _BEMBEL_TEST_WITH_ASSERT_
            assert(isIn(x, myto) or isIn(x - complexshift, myto));
#endif
            if (not(isIn(x, myto) or isIn(x - complexshift, myto))) {
              failcount = 1;
            }
          }

          // Check wheather the 1-to-1 relation has been identified correctly.

          const int sz = myto.size();
          for (int i = 0; i < sz; i++) {
// if (not gb.getPartnerId(myfrom[i])[0] == myto[i]) {
// std::cout << myfrom[i]<< " ++ "<< gb.getPartnerId(myfrom[i])[0]
// << " == " << myto[i] << std::endl;
#ifdef _BEMBEL_TEST_WITH_ASSERT_
            assert(gb.getPartnerId(myfrom[i])[0] == myto[i]);

            // }
            // if (not gb.getPartnerId(myfrom[i] + complexshift)[0] ==
            // (myto[i] + complexshift)) {
            // std::cout << gb.getPartnerId(myfrom[i] + complexshift)[0]
            // << " == " << (myto[i] + complexshift) << std::endl;

            assert(gb.getPartnerId(myfrom[i] + complexshift)[0] ==
                   (myto[i] + complexshift));
#endif
            if (not(gb.getPartnerId(myfrom[i])[0] == myto[i] and
                    gb.getPartnerId(myfrom[i] + complexshift)[0] ==
                        (myto[i] + complexshift))) {
              failcount = 1;
            }
            // }
          }

          // free_discretization(&disc);
          // free_meshdata(&mesh);
        }

  return failcount;
}
