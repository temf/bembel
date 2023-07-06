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
#ifndef BEMBEL_SRC_ANSATZSPACE_SUPERSPACE_HPP_
#define BEMBEL_SRC_ANSATZSPACE_SUPERSPACE_HPP_
namespace Bembel {
/**
 *  \ingroup AnsatzSpace
 *  \brief The superspace manages local polynomial bases on each element of the
 * mesh and provides an itnerface to evaluate them.
 */
template <typename Derived>
struct SuperSpace {
  typedef typename LinearOperatorTraits<Derived>::Scalar Scalar;
  //////////////////////////////////////////////////////////////////////////////
  //    constructors
  //////////////////////////////////////////////////////////////////////////////
  SuperSpace() {}
  SuperSpace(Geometry& geom, int M, int P) { init_SuperSpace(geom, M, P); }
  SuperSpace(const SuperSpace& other) {
    mesh_ = other.mesh_;
    phi = other.phi;
    phiDx = other.phiDx;
    phiPhi = other.phiPhi;
    phiPhiDx = other.phiPhiDx;
    phiPhiDy = other.phiPhiDy;
    phiTimesPhi = other.phiTimesPhi;
    // vPhiScalVPhi = other.vPhiScalVPhi;
    divPhiTimesDivPhi = other.divPhiTimesDivPhi;
    polynomial_degree = other.polynomial_degree;
    polynomial_degree_plus_one_squared =
        other.polynomial_degree_plus_one_squared;
  }
  SuperSpace(SuperSpace&& other) {
    mesh_ = other.mesh_;
    phi = other.phi;
    phiDx = other.phiDx;
    phiPhi = other.phiPhi;
    phiPhiDx = other.phiPhiDx;
    phiPhiDy = other.phiPhiDy;
    phiTimesPhi = other.phiTimesPhi;
    // vPhiScalVPhi = other.vPhiScalVPhi;
    divPhiTimesDivPhi = other.divPhiTimesDivPhi;
    polynomial_degree = other.polynomial_degree;
    polynomial_degree_plus_one_squared =
        other.polynomial_degree_plus_one_squared;
  }
  SuperSpace& operator=(SuperSpace other) {
    mesh_ = other.mesh_;
    phi = other.phi;
    phiDx = other.phiDx;
    phiPhi = other.phiPhi;
    phiPhiDx = other.phiPhiDx;
    phiPhiDy = other.phiPhiDy;
    phiTimesPhi = other.phiTimesPhi;
    // vPhiScalVPhi = other.vPhiScalVPhi;
    divPhiTimesDivPhi = other.divPhiTimesDivPhi;
    polynomial_degree = other.polynomial_degree;
    polynomial_degree_plus_one_squared =
        other.polynomial_degree_plus_one_squared;
    return *this;
  }
  //////////////////////////////////////////////////////////////////////////////
  //    getters
  //////////////////////////////////////////////////////////////////////////////
  int get_polynomial_degree() const { return polynomial_degree; }
  int get_polynomial_degree_plus_one_squared() const {
    return polynomial_degree_plus_one_squared;
  }
  int get_refinement_level() const { return mesh_->get_max_level(); }
  int get_number_of_elements() const { return mesh_->get_number_of_elements(); }
  int get_number_of_patches() const { return mesh_->get_geometry().size(); }
  const PatchVector& get_geometry() const { return mesh_->get_geometry(); }
  const ClusterTree& get_mesh() const { return *mesh_; }
  //////////////////////////////////////////////////////////////////////////////
  //    init_SuperSpace
  //////////////////////////////////////////////////////////////////////////////
  void init_SuperSpace(const Geometry& geom, int M, int P) {
    polynomial_degree = P;
    polynomial_degree_plus_one_squared =
        (polynomial_degree + 1) * (polynomial_degree + 1);
    phi = (Basis::BasisHandler<Scalar>::funPtrPhi(P));
    phiDx = (Basis::BasisHandler<Scalar>::funPtrPhiDx(P));
    phiPhi = (Basis::BasisHandler<Scalar>::funPtrPhiPhi(P));
    phiPhiDx = (Basis::BasisHandler<Scalar>::funPtrPhiPhiDx(P));
    phiPhiDy = (Basis::BasisHandler<Scalar>::funPtrPhiPhiDy(P));
    phiTimesPhi = (Basis::BasisHandler<Scalar>::funPtrPhiTimesPhi(P));
    // vPhiScalVPhi = (Basis::BasisHandler<typename
    // LinearOperatorTraits<Derived>::Scalar>::funPtrVPhiScalVPhi(P));
    divPhiTimesDivPhi =
        (Basis::BasisHandler<Scalar>::funPtrDivPhiTimesDivPhi(P));
    mesh_ = std::make_shared<ClusterTree>();
    mesh_->init_ClusterTree(geom, M);
    mesh_->checkOrientation();
    return;
  }
  //////////////////////////////////////////////////////////////////////////////
  //    map2surface
  //////////////////////////////////////////////////////////////////////////////
  void map2surface(const ElementTreeNode& e, const Eigen::Vector2d& xi,
                   double w, SurfacePoint* surf_pt) const {
    Eigen::Vector2d st = e.llc_ + e.get_h() * xi;
    mesh_->get_geometry()[e.patch_].updateSurfacePoint(surf_pt, st, w, xi);
    return;
  }
  //////////////////////////////////////////////////////////////////////////////
  //    Methods
  //////////////////////////////////////////////////////////////////////////////
  /**
   * \brief Compute all products of local shape functions on the unit square at
   * coordinates s,t, scale by w and add to intval.
   */
  void addScaledBasisInteraction(
      Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>* intval,
      typename LinearOperatorTraits<Derived>::Scalar w,
      const Eigen::Vector2d& s, const Eigen::Vector2d& t) const {
    phiTimesPhi(intval, w, s, t);
  }
  /**
   * \brief Compute all products of local shape functions on the unit square at
   * coordinates s,t.
   */
  Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> basisInteraction(
      const Eigen::Vector2d& s, const Eigen::Vector2d& t) const {
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> intval(
        polynomial_degree_plus_one_squared, polynomial_degree_plus_one_squared);
    intval.setZero();
    phiTimesPhi(&intval, 1., s, t);
    return intval;
  }

  /**
   * \brief Compute all products of surface curls of local shape functions
   * on the unit square at coordinates s,t.
   */
  void addScaledSurfaceCurlInteraction(
      Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>* intval, Scalar w,
      const SurfacePoint& p1, const SurfacePoint& p2) const {
    // surface measures
    double kappa1 = p1.segment<3>(6).cross(p1.segment<3>(9)).norm();
    double kappa2 = p2.segment<3>(6).cross(p2.segment<3>(9)).norm();
    // compute basis functions's surface curl. Each column of s_curl is a basis
    // function's surface curl at point s.
    Eigen::MatrixXd s_curl(3, polynomial_degree_plus_one_squared);
    s_curl = (1.0 / kappa1) *
             (-p1.segment<3>(6) * basisDy(p1.segment<2>(0)).transpose() +
              p1.segment<3>(9) * basisDx(p1.segment<2>(0)).transpose());
    Eigen::MatrixXd t_curl(3, polynomial_degree_plus_one_squared);
    t_curl = (1.0 / kappa2) *
             (-p2.segment<3>(6) * basisDy(p2.segment<2>(0)).transpose() +
              p2.segment<3>(9) * basisDx(p2.segment<2>(0)).transpose());
    // inner product of surface curls of any two basis functions
    for (int j = 0; j < polynomial_degree_plus_one_squared; ++j)
      for (int i = 0; i < polynomial_degree_plus_one_squared; ++i)
        (*intval)(j * polynomial_degree_plus_one_squared + i) +=
            w * s_curl.col(i).dot(t_curl.col(j));
  }

  /**
   * \brief Compute all products of surface gradients of local shape functions
   * on the unit square at coordinates s,t.
   */
  void addScaledSurfaceGradientInteraction(
      Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>* intval, Scalar w,
      const SurfacePoint& p1, const SurfacePoint& p2) const {
    // inner product of surface gradients of any two basis functions equals to
    // inner product of surface curls of any two basis functions
    addScaledSurfaceCurlInteraction(intval, w, p1, p2);
  }

  /**
   * \brief Compute all scalar products of vector valued local shape functions
   * on the surface points with reference coordinates s,t, scale by w and add to
   * intval.
   */
  void addScaledVectorBasisInteraction(
      Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>* intval, Scalar w,
      const Eigen::Vector2d& s, const Eigen::Vector2d& t,
      const Eigen::Vector3d x_f_dx, const Eigen::Vector3d x_f_dy,
      const Eigen::Vector3d y_f_dx, const Eigen::Vector3d y_f_dy) const {
    auto basis_interaction = basisInteraction(s, t);
    intval->block(0, 0, polynomial_degree_plus_one_squared,
                  polynomial_degree_plus_one_squared) +=
        w * x_f_dx.dot(y_f_dx) * basis_interaction;
    intval->block(0, polynomial_degree_plus_one_squared,
                  polynomial_degree_plus_one_squared,
                  polynomial_degree_plus_one_squared) +=
        w * x_f_dx.dot(y_f_dy) * basis_interaction;
    intval->block(polynomial_degree_plus_one_squared, 0,
                  polynomial_degree_plus_one_squared,
                  polynomial_degree_plus_one_squared) +=
        w * x_f_dy.dot(y_f_dx) * basis_interaction;
    intval->block(polynomial_degree_plus_one_squared,
                  polynomial_degree_plus_one_squared,
                  polynomial_degree_plus_one_squared,
                  polynomial_degree_plus_one_squared) +=
        w * x_f_dy.dot(y_f_dy) * basis_interaction;
  }
  /**
   * \brief Compute all scalar products of vector valued local shape functions
   * on the surface points with reference coordinates s,t.
   */
  Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> vectorBasisInteraction(
      const Eigen::Vector2d& s, const Eigen::Vector2d& t,
      const Eigen::Vector3d x_f_dx, const Eigen::Vector3d x_f_dy,
      const Eigen::Vector3d y_f_dx, const Eigen::Vector3d y_f_dy) const {
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> intval(
        2 * polynomial_degree_plus_one_squared,
        2 * polynomial_degree_plus_one_squared);
    intval.setZero();
    addScaledVectorBasisInteraction(&intval, 1., s, t, x_f_dx, x_f_dy, y_f_dx,
                                    y_f_dy);
    return intval;
  }
  /**
   * \brief Compute all products of divergences of local shape functions on the
   * unit square at coordinates s,t, scale by w and add to intval.
   */
  void addScaledVectorBasisDivergenceInteraction(
      Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>* intval, Scalar w,
      const Eigen::Vector2d& s, const Eigen::Vector2d& t) const {
    divPhiTimesDivPhi(intval, w, s, t);
  }
  /**
   * \brief Compute all products of divergences of local shape functions on the
   * unit square at coordinates s,t.
   */
  Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>
  vectorBasisDivergenceInteraction(const Eigen::Vector2d& s,
                                   const Eigen::Vector2d& t) const {
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> intval(
        2 * polynomial_degree_plus_one_squared,
        2 * polynomial_degree_plus_one_squared);
    intval.setZero();
    divPhiTimesDivPhi(&intval, 1., s, t);
    return intval;
  }
  /**
   * \brief Evaluate local shape functions on the unit square at coordinate s,
   * scale by w and add to intval.
   */
  void addScaledBasis(Eigen::Matrix<Scalar, Eigen::Dynamic, 1>* intval,
                      Scalar w, const Eigen::Vector2d& s) const {
    phiPhi(intval, w, s);
  }
  /**
   * \brief Evaluate local shape functions on the unit square at coordinate s.
   */
  Eigen::Matrix<Scalar, Eigen::Dynamic, 1> basis(
      const Eigen::Vector2d& s) const {
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> intval(
        polynomial_degree_plus_one_squared);
    intval.setZero();
    phiPhi(&intval, 1., s);
    return intval;
  }
  /**
   * \brief Evaluate derivatives in x direction of local shape functions on the
   * unit square at coordinate s, scale by w and add to intval.
   */
  void addScaledBasisDx(Eigen::Matrix<Scalar, Eigen::Dynamic, 1>* intval,
                        typename LinearOperatorTraits<Derived>::Scalar w,
                        const Eigen::Vector2d& s) const {
    phiPhiDx(intval, w, s);
  }
  /**
   * \brief Evaluate derivatives in x direction of local shape functions on the
   * unit square at coordinate s.
   */
  Eigen::Matrix<Scalar, Eigen::Dynamic, 1> basisDx(
      const Eigen::Vector2d& s) const {
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> intval(
        polynomial_degree_plus_one_squared);
    intval.setZero();
    phiPhiDx(&intval, 1., s);
    return intval;
  }
  /**
   * \brief Evaluate derivatives in y direction of local shape functions on the
   * unit square at coordinate s, scale by w and add to intval.
   */
  void addScaledBasisDy(Eigen::Matrix<Scalar, Eigen::Dynamic, 1>* intval,
                        typename LinearOperatorTraits<Derived>::Scalar w,
                        const Eigen::Vector2d& s) const {
    phiPhiDy(intval, w, s);
  }
  /**
   * \brief Evaluate derivatives in y direction of local shape functions on the
   * unit square at coordinate s.
   */
  Eigen::Matrix<Scalar, Eigen::Dynamic, 1> basisDy(
      const Eigen::Vector2d& s) const {
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> intval(
        polynomial_degree_plus_one_squared);
    intval.setZero();
    phiPhiDy(&intval, 1., s);
    return intval;
  }
  /**
   * \brief Evaluate local shape functions on the unit interval at coordinate s,
   * scale by w and add to intval.
   */
  void addScaledBasis1D(Eigen::Matrix<Scalar, Eigen::Dynamic, 1>* intval,
                        Scalar w, double s) const {
    phi(intval, w, s);
  }
  /**
   * \brief Evaluate local shape functions on the unit interval at coordinate s.
   */
  Eigen::Matrix<Scalar, Eigen::Dynamic, 1> basis1D(double s) const {
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> intval(polynomial_degree + 1);
    intval.setZero();
    phi(&intval, 1., s);
    return intval;
  }
  /**
   * \brief Evaluate derivatives of local shape functions on the unit interval
   * at coordinate s, scale by w and add to intval.
   */
  void addScaledBasis1DDx(Eigen::Matrix<Scalar, Eigen::Dynamic, 1>* intval,
                          Scalar w, double s) const {
    phiDx(intval, w, s);
  }
  /**
   * \brief Evaluate derivatives of local shape functions on the unit interval
   * at coordinate s.
   */
  Eigen::Matrix<Scalar, Eigen::Dynamic, 1> basis1DDx(double s) const {
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> intval(polynomial_degree + 1);
    intval.setZero();
    phiDx(&intval, 1., s);
    return intval;
  }
  //////////////////////////////////////////////////////////////////////////////
  //    member variables
  //////////////////////////////////////////////////////////////////////////////
 private:
  std::shared_ptr<ClusterTree> mesh_;
  Basis::funptr_phi<Scalar> phi;
  Basis::funptr_phidx<Scalar> phiDx;
  Basis::funptr_phiphi<Scalar> phiPhi;
  Basis::funptr_phiphidx<Scalar> phiPhiDx;
  Basis::funptr_phiphidy<Scalar> phiPhiDy;
  Basis::funptr_phitimesphi<Scalar> phiTimesPhi;
  // Basis::funptr_vphiscalvphi<typename LinearOperatorTraits<Derived>::Scalar>
  // vPhiScalVPhi;
  Basis::funptr_divphitimesdivphi<Scalar> divPhiTimesDivPhi;
  int polynomial_degree;
  int polynomial_degree_plus_one_squared;
};  // namespace Bembel
}  // namespace Bembel
#endif  // BEMBEL_SRC_ANSATZSPACE_SUPERSPACE_HPP_
