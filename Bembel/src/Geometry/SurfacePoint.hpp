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
#ifndef BEMBEL_SRC_GEOMETRY_SURFACEPOINT_HPP_
#define BEMBEL_SRC_GEOMETRY_SURFACEPOINT_HPP_
// #include <Eigen/StdVector>
/**
 * \ingroup Geometry
 * \brief This class stores quadrature and geometry information on quadrature
 *points and provides various methods to compute differential geometric
 *entities.
 **/
class SurfacePoint {
 public:
  // quadrature point
  Eigen::Vector2d get_xi() { return xi_; }
  const Eigen::Vector2d get_xi() const { return xi_; }

  // quadrature weight
  double get_w() { return w_; }
  const double get_w() const { return w_; }

  // geometry parametrization
  Eigen::Vector3d get_f() { return f_; }
  const Eigen::Vector3d get_f() const { return f_; }

  // geometry parametrization first derivatives
  Eigen::Vector3d get_f_dx() { return jacobian_.col(0); }
  const Eigen::Vector3d get_f_dx() const { return jacobian_.col(0); }
  Eigen::Vector3d get_f_dy() { return jacobian_.col(1); }
  const Eigen::Vector3d get_f_dy() const { return jacobian_.col(1); }
  Eigen::Matrix<double, 3, 2> get_jacobian() { return jacobian_; }
  const Eigen::Matrix<double, 3, 2> get_jacobian() const { return jacobian_; }

  // normal vectors
  Eigen::Vector3d get_normal() {
    return jacobian_.col(0).cross(jacobian_.col(1));
  }
  const Eigen::Vector3d get_normal() const {
    return jacobian_.col(0).cross(jacobian_.col(1));
  }
  Eigen::Vector3d get_unit_normal() { return get_normal().normalized(); }
  const Eigen::Vector3d get_unit_normal() const {
    return get_normal().normalized();
  }

  // surface measure
  double get_surface_measure() { return get_normal().norm(); }
  const double get_surface_measure() const { return get_normal().norm(); }

  // setter
  void set_xi(const Eigen::Vector2d &xi) { xi_ = xi; }
  void set_w(const double w) { w_ = w; }
  void set_f(const Eigen::Vector3d &f) { f_ = f; }
  void set_jacobian(const Eigen::Matrix<double, 3, 2> &jacobian) {
    jacobian_ = jacobian;
  }

 private:
  Eigen::Vector2d xi_;
  double w_;
  Eigen::Vector3d f_;
  Eigen::Matrix<double, 3, 2> jacobian_;
};

typedef std::vector<SurfacePoint, Eigen::aligned_allocator<SurfacePoint>>
    ElementSurfacePoints;

#endif  // BEMBEL_SRC_GEOMETRY_SURFACEPOINT_HPP_
