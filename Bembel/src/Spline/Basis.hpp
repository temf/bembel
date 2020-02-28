// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.

#ifndef BEMBEL_SPLINE_BASIS_H_
#define BEMBEL_SPLINE_BASIS_H_

namespace Bembel {
/**
 *  \ingroup Spline
 *  \brief The Basis namespace contains classes and functions that are to be
 * used as an interface between the BEM code and the functions in the Spl
 * namespace.
 */
namespace Basis {

// These typedefs are required for the superspace to store the correct functions
template <typename Scalar>
using funptr_voidOut_scalarptrScalarDoubleIn =
    void (*)(Eigen::Matrix<Scalar, Eigen::Dynamic, 1> *, Scalar, double);
template <typename Scalar>
using funptr_voidOut_scalarptrScalarVec2In = void (*)(
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> *, Scalar, Eigen::Vector2d);
template <typename Scalar>
using funptr_voidOut_scalarptrScalarVec2Vec2In =
    void (*)(Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> *, Scalar,
             Eigen::Vector2d, Eigen::Vector2d);

// These typedefs are a convenience to make the above human-readable
template <typename Scalar>
using funptr_phi = funptr_voidOut_scalarptrScalarDoubleIn<Scalar>;
template <typename Scalar>
using funptr_phidx = funptr_voidOut_scalarptrScalarDoubleIn<Scalar>;
template <typename Scalar>
using funptr_phiphi = funptr_voidOut_scalarptrScalarVec2In<Scalar>;
template <typename Scalar>
using funptr_phiphidx = funptr_voidOut_scalarptrScalarVec2In<Scalar>;
template <typename Scalar>
using funptr_phiphidy = funptr_voidOut_scalarptrScalarVec2In<Scalar>;
template <typename Scalar>
using funptr_phitimesphi = funptr_voidOut_scalarptrScalarVec2Vec2In<Scalar>;
template <typename Scalar>
using funptr_divphitimesdivphi =
    funptr_voidOut_scalarptrScalarVec2Vec2In<Scalar>;

/**
 * \ingroup Spline
 * \brief evaluates the 1D basis at x weighted with a quadrature weight w
 **/
template <int P, typename Scalar>
inline void phi_(Eigen::Matrix<Scalar, Eigen::Dynamic, 1> *c, Scalar w,
                 double x) {
  constexpr int I = P + 1;
  double base[I];
  PSpecificShapeFunctionHandler<I - 1>::evalBasis(I - 1, base, x);
  for (int i = 0; i < I; i++) (*c)(i) += w * base[i];
  return;
}

/**
 * \ingroup Spline
 * \brief evaluates the derivative of phi
 **/
template <int P, typename Scalar>
inline void phi_dx_(Eigen::Matrix<Scalar, Eigen::Dynamic, 1> *c, Scalar w,
                    double x) {
  constexpr int I = P + 1;
  double base[I];
  PSpecificShapeFunctionHandler<I - 1>::evalDerBasis(I - 1, base, x);
  for (int i = 0; i < I; i++) (*c)(i) += w * base[i];
  return;
}

/**
 * \ingroup Spline
 * \brief evaluates the 2D tensor product basis at a point a in [0,1]^2
 **/
template <int P, typename Scalar>
inline void phiphi_(Eigen::Matrix<Scalar, Eigen::Dynamic, 1> *c, Scalar w,
                    Eigen::Vector2d a) {
  constexpr int I = P + 1;
  double X[I], Y[I];
  PSpecificShapeFunctionHandler<I - 1>::evalBasis(I - 1, X, a(0));
  PSpecificShapeFunctionHandler<I - 1>::evalBasis(I - 1, Y, a(1));

  for (int iy = 0; iy < I; iy++)
    for (int ix = 0; ix < I; ix++) (*c)(iy * I + ix) += w * X[ix] * Y[iy];

  return;
}

/**
 * \ingroup Spline
 * \brief evaluates the x-derivative of the 2D tensor product basis at a
 *         point a in [0,1]^2
 **/
template <int P, typename Scalar>
inline void phiphi_dx_(Eigen::Matrix<Scalar, Eigen::Dynamic, 1> *c, Scalar w,
                       Eigen::Vector2d a) {
  constexpr int I = P + 1;
  double dX[I], Y[I];
  PSpecificShapeFunctionHandler<I - 1>::evalDerBasis(I - 1, dX, a(0));
  PSpecificShapeFunctionHandler<I - 1>::evalBasis(I - 1, Y, a(1));

  for (int iy = 0; iy < I; iy++)
    for (int ix = 0; ix < I; ix++) (*c)(iy * I + ix) += w * dX[ix] * Y[iy];

  return;
}

/**
 * \ingroup Spline
 * \brief evaluates the y-derivative of the 2D tensor product basis at a
 *         point a in [0,1]^2
 **/
template <int P, typename Scalar>
inline void phiphi_dy_(Eigen::Matrix<Scalar, Eigen::Dynamic, 1> *c, Scalar w,
                       Eigen::Vector2d a) {
  constexpr int I = P + 1;
  double X[I], dY[I];
  PSpecificShapeFunctionHandler<I - 1>::evalBasis(I - 1, X, a(0));
  PSpecificShapeFunctionHandler<I - 1>::evalDerBasis(I - 1, dY, a(1));

  for (int iy = 0; iy < I; iy++)
    for (int ix = 0; ix < I; ix++) (*c)(iy * I + ix) += w * X[ix] * dY[iy];
  return;
}

/**
 * \ingroup Spline
 * \brief evaluates the interaction of two phiphis, one at xi and one at
 *        eta. Used for e.g. gram matrices.
 **/
template <int P, typename Scalar>
void Phi_times_Phi_(Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> *c,
                    Scalar w, Eigen::Vector2d xi, Eigen::Vector2d eta) {
  constexpr int I = P + 1;
  Scalar a[I * I];
  double b[I * I], X[I], Y[I];

  PSpecificShapeFunctionHandler<I - 1>::evalBasis(I - 1, X, xi(0));
  PSpecificShapeFunctionHandler<I - 1>::evalBasis(I - 1, Y, xi(1));

  for (int iy = 0; iy < I; iy++)
    for (int ix = 0; ix < I; ix++) a[iy * I + ix] = w * X[ix] * Y[iy];

  PSpecificShapeFunctionHandler<I - 1>::evalBasis(I - 1, X, eta(0));
  PSpecificShapeFunctionHandler<I - 1>::evalBasis(I - 1, Y, eta(1));

  for (int iy = 0; iy < I; iy++)
    for (int ix = 0; ix < I; ix++) b[iy * I + ix] = X[ix] * Y[iy];

  for (int i = 0; i < (I * I); i++)
    for (int j = 0; j < (I * I); j++) (*c)(i, j) += a[i] * b[j];

  return;
}

/**
 * \ingroup Spline
 * \brief same as above, just using the divergence
 **/
template <int P, typename Scalar>
void Div_Phi_times_Div_Phi_(
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> *c, Scalar weight,
    Eigen::Vector2d xi, Eigen::Vector2d eta) {
  constexpr int I = P + 1;
  constexpr int I2 = I * I;
  Eigen::Matrix<Scalar, Eigen::Dynamic, 1> a_dx(I2);
  a_dx.setZero();
  Eigen::Matrix<Scalar, Eigen::Dynamic, 1> a_dy(I2);
  a_dy.setZero();
  Eigen::VectorXd b_dx(I2);
  b_dx.setZero();
  Eigen::VectorXd b_dy(I2);
  b_dy.setZero();

  phiphi_dx_<P>(&a_dx, weight, xi);
  phiphi_dy_<P>(&a_dy, weight, xi);
  phiphi_dx_<P>(&b_dx, 1., eta);
  phiphi_dy_<P>(&b_dy, 1., eta);

  for (int i = 0; i < I2; ++i)
    for (int j = 0; j < I2; ++j) (*c)(i, j) += a_dx[i] * b_dx[j];
  for (int i = 0; i < I2; ++i)
    for (int j = 0; j < I2; ++j) (*c)(i, j + I2) += a_dx[i] * b_dy[j];
  for (int i = 0; i < I2; ++i)
    for (int j = 0; j < I2; ++j) (*c)(i + I2, j) += a_dy[i] * b_dx[j];
  for (int i = 0; i < I2; ++i)
    for (int j = 0; j < I2; ++j) (*c)(i + I2, j + I2) += a_dy[i] * b_dy[j];

  return;
}

/**
 * \ingroup Spline
 * \brief The functions above have a fixed compile time polynomial degree. The
 *        PSpecificBasis handler is used to convert this to a runtime p through
 *        template recursion.
 **/
template <int P, typename Scalar>
class PSpecificBasisHandler : public PSpecificShapeFunctionHandler<P> {
 public:
  // These methods are for calling the functions
  static inline void phi(int p, Eigen::Matrix<Scalar, Eigen::Dynamic, 1> *c,
                         Scalar w, double x) {
    return P == p ? phi_<P, Scalar>(c, w, x)
                  : PSpecificBasisHandler<P - 1, Scalar>::phi(p, c, w, x);
  }
  static inline void phiDx(int p, Eigen::Matrix<Scalar, Eigen::Dynamic, 1> *c,
                           Scalar w, double x) {
    return P == p ? phi_dx_<P, Scalar>(c, w, x)
                  : PSpecificBasisHandler<P - 1, Scalar>::phiDx(p, c, w, x);
  }
  static inline void phiPhi(int p, Eigen::Matrix<Scalar, Eigen::Dynamic, 1> *c,
                            Scalar w, Eigen::Vector2d a) {
    return P == p ? phiphi_<P, Scalar>(c, w, a)
                  : PSpecificBasisHandler<P - 1, Scalar>::phiPhi(p, c, w, a);
  }
  static inline void phiPhiDx(int p,
                              Eigen::Matrix<Scalar, Eigen::Dynamic, 1> *c,
                              Scalar w, Eigen::Vector2d a) {
    return P == p ? phiphi_dx_<P, Scalar>(c, w, a)
                  : PSpecificBasisHandler<P - 1, Scalar>::phiPhiDx(p, c, w, a);
  }
  static inline void phiPhiDy(int p,
                              Eigen::Matrix<Scalar, Eigen::Dynamic, 1> *c,
                              Scalar w, Eigen::Vector2d a) {
    return P == p ? phiphi_dy_<P, Scalar>(c, w, a)
                  : PSpecificBasisHandler<P - 1, Scalar>::phiPhiDy(p, c, w, a);
  }
  static inline void phiTimesPhi(
      int p, Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> *c, Scalar w,
      Eigen::Vector2d xi, Eigen::Vector2d eta) {
    return P == p ? Phi_times_Phi_<P, Scalar>(c, w, xi, eta)
                  : PSpecificBasisHandler<P - 1, Scalar>::phiTimesPhi(p, c, w,
                                                                      xi, eta);
  }
  static inline void divPhiTimesDivPhi(
      int p, Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> *c,
      Scalar weight, Eigen::Vector2d xi, Eigen::Vector2d eta) {
    return P == p ? Div_Phi_times_Div_Phi_<P>(c, weight, xi, eta)
                  : PSpecificBasisHandler<P - 1, Scalar>::divPhiTimesDivPhi(
                        p, c, weight, xi, eta);
  }
  // These methods are for storing the functions for a given p, i.e., return
  // the function pointers
  static constexpr funptr_phi<Scalar> funPtrPhi(int p) {
    return P == p ? &phi_<P, Scalar>
                  : PSpecificBasisHandler<P - 1, Scalar>::funPtrPhi(p);
  }
  static constexpr funptr_phidx<Scalar> funPtrPhiDx(int p) {
    return P == p ? &phi_dx_<P, Scalar>
                  : PSpecificBasisHandler<P - 1, Scalar>::funPtrPhiDx(p);
  }
  static constexpr funptr_phiphi<Scalar> funPtrPhiPhi(int p) {
    return P == p ? &phiphi_<P, Scalar>
                  : PSpecificBasisHandler<P - 1, Scalar>::funPtrPhiPhi(p);
  }
  static constexpr funptr_phiphidx<Scalar> funPtrPhiPhiDx(int p) {
    return P == p ? &phiphi_dx_<P, Scalar>
                  : PSpecificBasisHandler<P - 1, Scalar>::funPtrPhiPhiDx(p);
  }
  static constexpr funptr_phiphidy<Scalar> funPtrPhiPhiDy(int p) {
    return P == p ? &phiphi_dy_<P, Scalar>
                  : PSpecificBasisHandler<P - 1, Scalar>::funPtrPhiPhiDy(p);
  }
  static constexpr funptr_phitimesphi<Scalar> funPtrPhiTimesPhi(int p) {
    return P == p ? &Phi_times_Phi_<P, Scalar>
                  : PSpecificBasisHandler<P - 1, Scalar>::funPtrPhiTimesPhi(p);
  }
  static constexpr funptr_divphitimesdivphi<Scalar> funPtrDivPhiTimesDivPhi(
      int p) {
    return P == p
               ? &Div_Phi_times_Div_Phi_<P, Scalar>
               : PSpecificBasisHandler<P - 1, Scalar>::funPtrDivPhiTimesDivPhi(
                     p);
  }
};

// Anchors of the recursions above
template <typename Scalar>
class PSpecificBasisHandler<0, Scalar>
    : public PSpecificShapeFunctionHandler<0> {
 public:
  static inline void phi(int p, Eigen::Matrix<Scalar, Eigen::Dynamic, 1> *c,
                         Scalar w, double x) {
    return phi_<0, Scalar>(c, w, x);
  }
  static inline void phiDx(int p, Eigen::Matrix<Scalar, Eigen::Dynamic, 1> *c,
                           Scalar w, double x) {
    return phi_dx_<0, Scalar>(c, w, x);
  }
  static inline void phiPhi(int p, Eigen::Matrix<Scalar, Eigen::Dynamic, 1> *c,
                            Scalar w, Eigen::Vector2d a) {
    return phiphi_<0, Scalar>(c, w, a);
  }
  static inline void phiPhiDx(int p,
                              Eigen::Matrix<Scalar, Eigen::Dynamic, 1> *c,
                              Scalar w, Eigen::Vector2d a) {
    return phiphi_dx_<0, Scalar>(c, w, a);
  }
  static inline void phiPhiDy(int p,
                              Eigen::Matrix<Scalar, Eigen::Dynamic, 1> *c,
                              Scalar w, Eigen::Vector2d a) {
    return phiphi_dy_<0, Scalar>(c, w, a);
  }
  static inline void phiTimesPhi(
      int p, Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> *c, Scalar w,
      Eigen::Vector2d xi, Eigen::Vector2d eta) {
    return Phi_times_Phi_<0, Scalar>(c, w, xi, eta);
  }
  static inline void divPhiTimesDivPhi(
      int p, Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> *c,
      Scalar weight, Eigen::Vector2d xi, Eigen::Vector2d eta) {
    return Div_Phi_times_Div_Phi_<0, Scalar>(c, weight, xi, eta);
  }

  static constexpr funptr_phi<Scalar> funPtrPhi(int p) {
    return &phi_<0, Scalar>;
  }
  static constexpr funptr_phidx<Scalar> funPtrPhiDx(int p) {
    return &phi_dx_<0, Scalar>;
  }
  static constexpr funptr_phiphi<Scalar> funPtrPhiPhi(int p) {
    return &phiphi_<0, Scalar>;
  }
  static constexpr funptr_phiphidx<Scalar> funPtrPhiPhiDx(int p) {
    return &phiphi_dx_<0, Scalar>;
  }
  static constexpr funptr_phiphidy<Scalar> funPtrPhiPhiDy(int p) {
    return &phiphi_dy_<0, Scalar>;
  }
  static constexpr funptr_phitimesphi<Scalar> funPtrPhiTimesPhi(int p) {
    return &Phi_times_Phi_<0, Scalar>;
  }
  static constexpr funptr_divphitimesdivphi<Scalar> funPtrDivPhiTimesDivPhi(
      int p) {
    return &Div_Phi_times_Div_Phi_<0, Scalar>;
  }
};

/// This instantiates the basishandler for a maximal p
template <typename Scalar>
using BasisHandler = PSpecificBasisHandler<Constants::MaxP, Scalar>;

}  // namespace Basis
}  // namespace Bembel

#endif
