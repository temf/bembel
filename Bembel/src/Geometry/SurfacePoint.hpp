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
 * \brief Wrapper class for legacy SurfacePoint represented as an Eigen::Vector
 *
 * This typedef is essential for any evaluation of a bilinear form. It provides
 * all required geometry information, stored as follows:
 * (0) x-coordinate of the evaluation point in the parameter domain [0,1]^2
 *     of the current element, i.e we map [0,1]^2->element->surface
 * (1) y-coordinate of the evaluation point in the parameter domain [0,1]^2
 *     of the current element, i.e we map [0,1]^2->element->surface
 * (2) a quadrature weight. Can be left empty if not used as part of a
 *     quadrature.
 * (3) x-coordinate of patch eval in space
 * (4) y-coordinate of patch eval in space
 * (5) z-coordinate of patch eval in space
 * (6) x-component of derivative in x-dir
 * (7) y-component of derivative in x-dir
 * (8) z-component of derivative in x-dir
 * (9) x-component of derivative in y-dir
 * (10) y-component of derivative in y-dir
 * (11) z-component of derivative in y-dir
 * For application of the pull-back to the reference domain, one requires the
 * jacobian of any point on the surface. Calling eval and evalJacobian of the
 * Patch class introduces work that needs to be done twice. The
 * updateSurdacePoint method is specialized and should be used, since it avoids
 * redundant work.
 **/
class SurfacePoint {
 public:
  /**
   * \brief Supports legacy code.
   */
  double& operator()(const int j) { return (data_(j)); }
  /**
   * \brief Supports legacy code.
   */
  const double& operator()(const int j) const { return (data_(j)); }
  /**
   * \brief Supports legacy code.
   */
  template <int N>
  Eigen::Matrix<double, N, 1> segment(const int l) {
    return data_.segment<N>(l);
  }
  /**
   * \brief Supports legacy code.
   */
  template <int N>
  const Eigen::Matrix<double, N, 1> segment(const int l) const {
    return data_.segment<N>(l);
  }
  /**
   * \brief Supports legacy code.
   */
  Eigen::Matrix<double, Eigen::Dynamic, 1> segment(const int l, const int n) {
    return data_.segment(l, n);
  }
  /**
   * \brief Supports legacy code.
   */
  const Eigen::Matrix<double, Eigen::Dynamic, 1> segment(const int l,
                                                         const int n) const {
    return data_.segment(l, n);
  }
  /**
   * \brief Supports legacy code.
   */
  Eigen::Matrix<double, Eigen::Dynamic, 1> head(const int n) {
    return data_.head(n);
  }
  /**
   * \brief Supports legacy code.
   */
  const Eigen::Matrix<double, Eigen::Dynamic, 1> head(const int n) const {
    return data_.head(n);
  }
  /**
   * \brief Supports legacy code.
   */
  Eigen::Matrix<double, Eigen::Dynamic, 1> tail(const int n) {
    return data_.tail(n);
  }
  /**
   * \brief Supports legacy code.
   */
  const Eigen::Matrix<double, Eigen::Dynamic, 1> tail(const int n) const {
    return data_.tail(n);
  }
  /**
   * \brief Supports legacy code.
   */
  Eigen::Matrix<double, Eigen::Dynamic, 1> get_data() { return data_; }
  /**
   * \brief Supports legacy code.
   */
  const Eigen::Matrix<double, Eigen::Dynamic, 1> get_data() const {
    return data_;
  }

  // EIGEN_MAKE_ALIGNED_OPERATOR_NEW

 private:
  Eigen::Matrix<double, 12, 1> data_;
};

typedef std::vector<SurfacePoint, Eigen::aligned_allocator<SurfacePoint>>
    ElementSurfacePoints;

#endif  // BEMBEL_SRC_GEOMETRY_SURFACEPOINT_HPP_