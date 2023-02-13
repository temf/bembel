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
#ifndef BEMBEL_SRC_GEOMETRY_PATCH_HPP_
#define BEMBEL_SRC_GEOMETRY_PATCH_HPP_
namespace Bembel {

/**
 *  \ingroup Geometry
 *  \class Patch
 *  \brief handles a single patch
 **/
class Patch {
 public:
  //////////////////////////////////////////////////////////////////////////////
  /// constructors
  //////////////////////////////////////////////////////////////////////////////
  Patch() {}
  Patch(const std::vector<Eigen::Matrix<double, -1, -1>> &control_points,
        const std::vector<double> &knots_x,
        const std::vector<double> &knots_y) {
    init_Patch(control_points, knots_x, knots_y);
  }
  //////////////////////////////////////////////////////////////////////////////
  /// init
  //////////////////////////////////////////////////////////////////////////////
  inline void init_Patch(const std::vector<Eigen::Matrix<double, -1, -1>> &xyzw,
                         const std::vector<double> &x_knots,
                         const std::vector<double> &y_knots) {
    assert(xyzw.size() == 4);
    const int xyzw_cols = xyzw[0].cols();
    const int xyzw_rows = xyzw[0].rows();
    unique_knots_x_ = Spl::ExtractUniqueKnotVector(x_knots);
    unique_knots_y_ = Spl::ExtractUniqueKnotVector(y_knots);
    const int xnumpatch = unique_knots_x_.size() - 1;
    const int ynumpatch = unique_knots_y_.size() - 1;
    polynomial_degree_x_ = x_knots.size() - xyzw_cols;
    polynomial_degree_y_ = y_knots.size() - xyzw_rows;
    data_.resize(4 * (polynomial_degree_x_ * xnumpatch * polynomial_degree_y_ *
                      ynumpatch));
    {
      // Since its only for initialization, I do not care about speed.
      // Here I look weather the data given is already in bezier form.
      if (unique_knots_x_.size() == 2 && unique_knots_y_.size() == 2) {
        for (int i = 0; i < 4; i++) {
          Eigen::Matrix<double, -1, 1> tmp = Spl::Unroll(xyzw[i]);
          for (int j = 0; j < tmp.rows(); j++) data_[j * 4 + i] = (tmp[j]);
        }
      } else {
        // If not, I construct the dynamic projection (i.e. solve
        // systems for
        // the coeffs) and project to the superspace.

        Eigen::SparseMatrix<double> phi = Spl::MakeProjection(
            x_knots, y_knots, unique_knots_x_, unique_knots_y_,
            polynomial_degree_x_, polynomial_degree_y_);

        for (int i = 0; i < 4; i++) {
          Eigen::Matrix<double, -1, 1> tmp =
              Spl::Unroll(xyzw[i]).transpose() * phi.transpose();

          for (int j = 0; j < tmp.rows(); j++) data_[j * 4 + i] = (tmp[j]);
        }
      }
    }

    return;
  }
  /* eval() evaluates the geometry. I look up the position in the knot vector,
   * scale the input arguments, evaluate the 1D basis functions and sum over
   * them with the controll points from data. */

  Eigen::Vector3d eval(const Eigen::Vector2d &reference_point) const {
    const int x_location =
        Spl::FindLocationInKnotVector(reference_point(0), unique_knots_x_);
    const int y_location =
        Spl::FindLocationInKnotVector(reference_point(1), unique_knots_y_);
    const int numy = (unique_knots_y_.size() - 1) * polynomial_degree_y_;
    const double scaledx =
        Spl::Rescale(reference_point(0), unique_knots_x_[x_location],
                     unique_knots_x_[x_location + 1]);
    const double scaledy =
        Spl::Rescale(reference_point(1), unique_knots_y_[y_location],
                     unique_knots_y_[y_location + 1]);

    Eigen::Vector4d tmp = Eigen::Vector4d::Zero();
    Eigen::VectorXd xbasis = Bembel::Basis::ShapeFunctionHandler::evalBasis(
        polynomial_degree_x_ - 1, scaledx);
    Eigen::VectorXd ybasis = Bembel::Basis::ShapeFunctionHandler::evalBasis(
        polynomial_degree_y_ - 1, scaledy);

    for (int i = 0; i < polynomial_degree_x_; i++) {
      for (int j = 0; j < polynomial_degree_y_; j++) {
        const double tpbasisval = xbasis(i) * ybasis(j);
        const int accs = 4 * (numy * (polynomial_degree_x_ * x_location + i) +
                              polynomial_degree_y_ * y_location + j);
        for (int k = 0; k < 4; k++) tmp(k) += data_[accs + k] * tpbasisval;
      }
    }

    return tmp.head<3>() / tmp[3];
  }

  Eigen::Matrix<double, 3, 2> evalJacobian(
      const Eigen::Vector2d &reference_point) const {
    const int x_location =
        Spl::FindLocationInKnotVector(reference_point(0), unique_knots_x_);
    const int y_location =
        Spl::FindLocationInKnotVector(reference_point(1), unique_knots_y_);
    const int numy = (unique_knots_y_.size() - 1) * polynomial_degree_y_;
    const double scaledx =
        Spl::Rescale(reference_point(0), unique_knots_x_[x_location],
                     unique_knots_x_[x_location + 1]);
    const double scaledy =
        Spl::Rescale(reference_point(1), unique_knots_y_[y_location],
                     unique_knots_y_[y_location + 1]);

    Eigen::Vector4d tmp = Eigen::Vector4d::Zero();
    Eigen::Vector4d tmpDx = Eigen::Vector4d::Zero();
    Eigen::Vector4d tmpDy = Eigen::Vector4d::Zero();
    Eigen::VectorXd xbasis = Bembel::Basis::ShapeFunctionHandler::evalBasis(
        polynomial_degree_x_ - 1, scaledx);
    Eigen::VectorXd ybasis = Bembel::Basis::ShapeFunctionHandler::evalBasis(
        polynomial_degree_y_ - 1, scaledy);
    Eigen::VectorXd xbasisD = Bembel::Basis::ShapeFunctionHandler::evalDerBasis(
        polynomial_degree_x_ - 1, scaledx);
    Eigen::VectorXd ybasisD = Bembel::Basis::ShapeFunctionHandler::evalDerBasis(
        polynomial_degree_y_ - 1, scaledy);

    for (int i = 0; i < polynomial_degree_x_; i++) {
      for (int j = 0; j < polynomial_degree_y_; j++) {
        const double tpbasisval = xbasis(i) * ybasis(j);
        const double tpbasisvalDx = xbasisD(i) * ybasis(j);
        const double tpbasisvalDy = xbasis(i) * ybasisD(j);
        const int accs = 4 * (numy * (polynomial_degree_x_ * x_location + i) +
                              polynomial_degree_y_ * y_location + j);

        // Here I add up the values of the basis functions in the dc
        // basis
        for (int k = 0; k < 4; k++) {
          tmp(k) += data_[accs + k] * tpbasisval;
          tmpDx(k) += data_[accs + k] * tpbasisvalDx;
          tmpDy(k) += data_[accs + k] * tpbasisvalDy;
        }
      }
    }

    double botsqr = 1. / (tmp(3) * tmp(3));
    return botsqr * (Eigen::Matrix<double, 3, 2>()
                         << tmpDx.head<3>() * tmp(3) - tmp.head<3>() * tmpDx(3),
                     tmpDy.head<3>() * tmp(3) - tmp.head<3>() * tmpDy(3))
                        .finished();
  }

  inline Eigen::Matrix<double, 3, 1> evalNormal(
      const Eigen::Vector2d &reference_point) const {
    Eigen::Matrix<double, 3, 2> jac = evalJacobian(reference_point);
    return jac.col(0).cross(jac.col(1));
  }

  // Wrapper for legacy code
  inline Eigen::Vector3d eval(double x, double y) const {
    return eval(Eigen::Vector2d(x, y));
  }
  inline Eigen::Matrix<double, 3, 2> evalJacobian(double x, double y) const {
    return evalJacobian(Eigen::Vector2d(x, y));
  }
  inline Eigen::Matrix<double, 3, 1> evalNormal(double x, double y) const {
    return evalNormal(Eigen::Vector2d(x, y));
  }

  // This is a combination of eval und evalJacobian, to avoid duplication of
  // work. See SurfacePoint.hpp
  void updateSurfacePoint(SurfacePoint *srf_pt, const Eigen::Vector2d &ref_pt,
                          double w, const Eigen::Vector2d &xi) const {
    const int x_location =
        Spl::FindLocationInKnotVector(ref_pt(0), unique_knots_x_);
    const int y_location =
        Spl::FindLocationInKnotVector(ref_pt(1), unique_knots_y_);
    const int numy = (unique_knots_y_.size() - 1) * polynomial_degree_y_;
    const double scaledx = Spl::Rescale(ref_pt(0), unique_knots_x_[x_location],
                                        unique_knots_x_[x_location + 1]);
    const double scaledy = Spl::Rescale(ref_pt(1), unique_knots_y_[y_location],
                                        unique_knots_y_[y_location + 1]);

    // resize operations do nothing if the size does not change
    srf_pt->get_buffer().resize(6);
    srf_pt->get_buffer()[0].resize(4, 1);
    srf_pt->get_buffer()[1].resize(4, 2);
    srf_pt->get_buffer()[2].resize(polynomial_degree_x_, 1);
    srf_pt->get_buffer()[3].resize(polynomial_degree_y_, 1);
    srf_pt->get_buffer()[4].resize(polynomial_degree_x_, 1);
    srf_pt->get_buffer()[5].resize(polynomial_degree_y_, 1);

    // improve readability
    Eigen::MatrixXd &tmp = srf_pt->get_buffer()[0];
    Eigen::MatrixXd &tmpD = srf_pt->get_buffer()[1];
    Eigen::MatrixXd &xbasis = srf_pt->get_buffer()[2];
    Eigen::MatrixXd &ybasis = srf_pt->get_buffer()[3];
    Eigen::MatrixXd &xbasisD = srf_pt->get_buffer()[4];
    Eigen::MatrixXd &ybasisD = srf_pt->get_buffer()[5];

    tmp.setZero();
    tmpD.setZero();

    xbasis.col(0) = Bembel::Basis::ShapeFunctionHandler::evalBasis(
        polynomial_degree_x_ - 1, scaledx);
    ybasis.col(0) = Bembel::Basis::ShapeFunctionHandler::evalBasis(
        polynomial_degree_y_ - 1, scaledy);
    xbasisD.col(0) = Bembel::Basis::ShapeFunctionHandler::evalDerBasis(
        polynomial_degree_x_ - 1, scaledx);
    ybasisD.col(0) = Bembel::Basis::ShapeFunctionHandler::evalDerBasis(
        polynomial_degree_y_ - 1, scaledy);

    for (int i = 0; i < polynomial_degree_x_; ++i) {
      for (int j = 0; j < polynomial_degree_y_; ++j) {
        const double tpbasisval = xbasis(i) * ybasis(j);
        const double tpbasisvalDx = xbasisD(i) * ybasis(j);
        const double tpbasisvalDy = xbasis(i) * ybasisD(j);
        const int accs = 4 * (numy * (polynomial_degree_x_ * x_location + i) +
                              polynomial_degree_y_ * y_location + j);

        // Here I add up the values of the basis functions in the dc
        // basis
        for (int k = 0; k < 4; ++k) {
          tmp(k, 0) += data_[accs + k] * tpbasisval;
          tmpD(k, 0) += data_[accs + k] * tpbasisvalDx;
          tmpD(k, 1) += data_[accs + k] * tpbasisvalDy;
        }
      }
    }

    const double bot = 1. / tmp(3, 0);
    const double botsqr = bot * bot;

    srf_pt->set_xi(xi);
    srf_pt->set_w(w);
    srf_pt->set_f(bot * tmp.block<3, 1>(0, 0));
    // srf_pt->set_jacobian(
    //     botsqr * (Eigen::Matrix<double, 3, 2>()
    //                   << tmpDx.head<3>() * tmp(3) - tmp.head<3>() * tmpDx(3),
    //               tmpDy.head<3>() * tmp(3) - tmp.head<3>() * tmpDy(3))
    //                  .finished());
    srf_pt->set_jacobian(botsqr *
                         (tmpD.block<3, 2>(0, 0) * tmp(3, 0) -
                          tmp.block<3, 1>(0, 0) * tmpD.block<1, 2>(3, 0)));
    return;
  }

  //////////////////////////////////////////////////////////////////////////////
  /// getter
  //////////////////////////////////////////////////////////////////////////////

  std::vector<double> data_;  // Controllpoints in Bezier-Extracted Format.
  int polynomial_degree_x_;   // Degree in x
  int polynomial_degree_y_;   // Degree in y
  std::vector<double>
      unique_knots_x_;  // The knot vectors, where each knot is unique
  std::vector<double>
      unique_knots_y_;  // The knot vectors, where each knot is unique
};

inline std::vector<Patch> PatchShredder(const Patch &patch) noexcept {
  // Already a Bezier patch
  if (patch.unique_knots_y_.size() == 2 && patch.unique_knots_x_.size() == 2) {
    return {patch};
  }

  // number of subpatches in x and y directions
  const int xchips = patch.unique_knots_x_.size() - 1;
  const int ychips = patch.unique_knots_y_.size() - 1;

  const int xp = patch.polynomial_degree_x_;
  const int yp = patch.polynomial_degree_y_;
  const int numy = ychips * yp;

  std::vector<Patch> out(xchips * ychips);

  for (int ix = 0; ix < xchips; ix++) {
    for (int iy = 0; iy < ychips; iy++) {
      const int index = ix * ychips + iy;

      out[index].unique_knots_x_ = {0, 1};
      out[index].unique_knots_y_ = {0, 1};
      out[index].polynomial_degree_x_ = xp;
      out[index].polynomial_degree_y_ = yp;
      out[index].data_.reserve(xp * yp * 4);
    }
  }

  for (int ix = 0; ix < xchips; ix++) {
    for (int iy = 0; iy < ychips; iy++) {
      const int index = ix * ychips + iy;
      for (int jx = 0; jx < xp; jx++) {
        for (int jy = 0; jy < yp; jy++) {
          const int accs = 4 * (numy * (xp * ix + jx) + yp * iy + jy);
          for (int k = 0; k < 4; k++) {
            out[index].data_.push_back(patch.data_[accs + k]);
          }
        }
      }
    }
  }

  return out;
}

// Shredds a whole vector of Patches
inline std::vector<Patch> PatchShredder(
    const std::vector<Patch> &patches) noexcept {
  std::vector<Patch> out;
  const int input_size = patches.size();

  for (int i = 0; i < input_size; i++) {
    std::vector<Patch> tmp = PatchShredder(patches[i]);

    const int tmp_size = tmp.size();

    for (int j = 0; j < tmp_size; j++) out.push_back(tmp[j]);
  }

  out.shrink_to_fit();

  return out;
}

}  // namespace Bembel

#endif  // BEMBEL_SRC_GEOMETRY_PATCH_HPP_
