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
#ifndef BEMBEL_SRC_CLUSTERTREE_CLUSTERTREE_HPP_
#define BEMBEL_SRC_CLUSTERTREE_CLUSTERTREE_HPP_

namespace Bembel {

/**
 * \ingroup ClusterTree
 * \brief The ClusterTree class introduces an element structure on a Geometry
 * object. Note that we do not introduce a mesh in the classical sense, but only
 * introduce a system of local coordinates via an ElementTree.
 */
class ClusterTree {
 public:
  /**
   * \brief Copy constructor (deleted).
   *
   * \param other The ClusterTree instance to copy from.
   */
  ClusterTree(const ClusterTree& other) = delete;

  /**
   * \brief Move constructor (deleted).
   *
   * \param other The ClusterTree instance to move from.
   */
  ClusterTree(ClusterTree&& other) = delete;

  /**
   * \brief Copy assignment operator (deleted).
   *
   * \param other The ClusterTree instance to copy from.
   * \return Reference to this ClusterTree instance.
   */
  ClusterTree& operator=(const ClusterTree& other) = delete;

  /**
   * \brief Move assignment operator (deleted).
   *
   * \param other The ClusterTree instance to move from.
   * \return Reference to this ClusterTree instance.
   */
  ClusterTree& operator=(ClusterTree&& other) = delete;
  //////////////////////////////////////////////////////////////////////////////
  /// constructors
  //////////////////////////////////////////////////////////////////////////////
  /**
   * \brief Default constructor.
   */
  ClusterTree() {}
  /**
   * \brief Constructs a ClusterTree object for a given refinement level.
   *
   * Initializes a new instance of the ClusterTree class with the given geometry
   * and refinement level.
   *
   * \param geom The geometry object to construct a cluster tree on.
   * \param M refinement level of the ElementTree.
   */
  ClusterTree(const Geometry& geom, int M) { init_ClusterTree(geom, M); }
  //////////////////////////////////////////////////////////////////////////////
  /// init
  //////////////////////////////////////////////////////////////////////////////
  /**
   * \brief Initializes a ClusterTree object for a given refinement level.
   *
   * This function initializes the ElementTree for an uniform refinement level.
   * After that the point list is filled.
   *
   * \param geom The geometry object to construct a cluster tree on.
   * \param M refinement level of the ElementTree.
   */
  void init_ClusterTree(const Geometry& geom, int M) {
    element_tree_.init_ElementTree(geom, M);
    points_ = element_tree_.computeElementEnclosings();
    return;
  }
  //////////////////////////////////////////////////////////////////////////////
  /// getter
  //////////////////////////////////////////////////////////////////////////////
  /**
   * \brief Return reference to the ElementTree.
   *
   * \return Reference to the ElementTree.
   */
  ElementTree& get_element_tree() { return element_tree_; }
  /**
   * \brief Return const reference to the ElementTree.
   *
   * \return Const Reference to the ElementTree.
   */
  const ElementTree& get_element_tree() const { return element_tree_; }
  /**
   * \brief Return point list with coordinates.
   *
   * \return 3xN matrix with N number of points in the ElementTree.
   */
  const Eigen::MatrixXd& get_points() const { return points_; }
  /**
   * \brief Return const reference to the Geometry.
   *
   * \return Const Reference to the PatchVector.
   */
  const PatchVector& get_geometry() const {
    return element_tree_.get_geometry();
  }
  /**
   * \brief Return maximum level of refinement.
   *
   * \return Level of refinement.
   */
  int get_max_level() const { return element_tree_.get_max_level(); }
  /**
   * \brief Return number of elements in the ElementTree.
   *
   * \return Number of elements.
   */
  int get_number_of_elements() const {
    return element_tree_.get_number_of_elements();
  }
  //////////////////////////////////////////////////////////////////////////////
  /// member functions
  //////////////////////////////////////////////////////////////////////////////
  /**
   * \brief Check all patches at the mid points of the edges if the
   * parametrization and normal vector match.
   */
  void checkOrientation() {
    std::vector<std::array<int, 4>> edges = element_tree_.patchTopologyInfo();
    std::vector<Eigen::Vector2d> edge_midpoints = {
        Eigen::Vector2d(Constants::edgemps[0][0], Constants::edgemps[1][0]),  //
        Eigen::Vector2d(Constants::edgemps[0][1], Constants::edgemps[1][1]),  //
        Eigen::Vector2d(Constants::edgemps[0][2], Constants::edgemps[1][2]),  //
        Eigen::Vector2d(Constants::edgemps[0][3], Constants::edgemps[1][3])   //
    };
    const PatchVector& geo = element_tree_.get_geometry();
    for (auto edge : edges) {
      if (edge[2] > 0 && edge[3] > 0) {
        Eigen::Vector3d a = geo[edge[0]].eval(edge_midpoints[edge[2]]);
        Eigen::Vector3d b = geo[edge[1]].eval(edge_midpoints[edge[3]]);
        Eigen::Vector3d na = geo[edge[0]].evalNormal(edge_midpoints[edge[2]]);
        Eigen::Vector3d nb = geo[edge[1]].evalNormal(edge_midpoints[edge[3]]);
        assert((a - b).norm() < Constants::pt_comp_tolerance &&
               "These points should coincide according to the element tree");
        assert(a.dot(b) > 0 &&
               "Normals across patches are oriented the same way");
      }
    }
    return;
  }
  //////////////////////////////////////////////////////////////////////////////
  /// private members
  //////////////////////////////////////////////////////////////////////////////
 private:
  ElementTree element_tree_;
  Eigen::MatrixXd points_;
  //////////////////////////////////////////////////////////////////////////////
};
}  // namespace Bembel
#endif  // BEMBEL_SRC_CLUSTERTREE_CLUSTERTREE_HPP_
