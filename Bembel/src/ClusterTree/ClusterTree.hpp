// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef BEMBEL_CLUSTERTREE_CLUSTERTREE_H_
#define BEMBEL_CLUSTERTREE_CLUSTERTREE_H_

namespace Bembel {

/**
 *  \ingroup ClusterTree
 *  \brief The ClusterTree class introduces an element structure on a Geometry object.
 * Note that we do not introduce a mesh in the classical sense, but only
 * introduce a system of local coordinates via an ElementTree.
 */
class ClusterTree {
 public:
  //////////////////////////////////////////////////////////////////////////////
  /// constructors
  //////////////////////////////////////////////////////////////////////////////
  ClusterTree() {}
  ClusterTree(const Geometry& geom, int M) { init_ClusterTree(geom, M); }
  //////////////////////////////////////////////////////////////////////////////
  /// init
  //////////////////////////////////////////////////////////////////////////////
  void init_ClusterTree(const Geometry& geom, int M) {
    max_level_ = M;
    geom_ = geom.get_geometry_ptr();
    element_tree_.init_ElementTree(geom, max_level_);
    points_ = element_tree_.computeElementEnclosings();
    return;
  }
  //////////////////////////////////////////////////////////////////////////////
  /// getter
  //////////////////////////////////////////////////////////////////////////////
  const ElementTree& get_element_tree() const { return element_tree_; }
  const Eigen::MatrixXd& get_points() const { return points_; }
  const PatchVector& get_geometry() const { return *geom_; }
  int get_max_level() const { return max_level_; }
  int get_number_of_elements() const {
    return std::distance(element_tree_.cpbegin(), element_tree_.cpend());
  }
  //////////////////////////////////////////////////////////////////////////////
  /// member functions
  //////////////////////////////////////////////////////////////////////////////
  void checkOrientation() {
    std::vector<std::array<int, 4>> edges = element_tree_.patchTopologyInfo();
    std::vector<Eigen::Vector2d> edge_midpoints = {
        Eigen::Vector2d(Constants::edgemps[0][0], Constants::edgemps[1][0]),  //
        Eigen::Vector2d(Constants::edgemps[0][1], Constants::edgemps[1][1]),  //
        Eigen::Vector2d(Constants::edgemps[0][2], Constants::edgemps[1][2]),  //
        Eigen::Vector2d(Constants::edgemps[0][3], Constants::edgemps[1][3])   //
    };
    for (auto edge : edges) {
      if (edge[2] > 0 && edge[3] > 0) {
        Eigen::Vector3d a = (*geom_)[edge[0]].eval(edge_midpoints[edge[2]]);
        Eigen::Vector3d b = (*geom_)[edge[1]].eval(edge_midpoints[edge[3]]);
        Eigen::Vector3d na =
            (*geom_)[edge[0]].evalNormal(edge_midpoints[edge[2]]);
        Eigen::Vector3d nb =
            (*geom_)[edge[1]].evalNormal(edge_midpoints[edge[3]]);
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
  // we declare functionality which has not been implemented (yet)
  // to be private
  ElementTree element_tree_;
  ClusterTree(const ClusterTree& other);
  ClusterTree(ClusterTree&& other);
  ClusterTree& operator=(const ClusterTree& other);
  ClusterTree& operator=(ClusterTree&& other);
  Eigen::MatrixXd points_;
  int max_level_;
  std::shared_ptr<PatchVector> geom_;
  //////////////////////////////////////////////////////////////////////////////
};
}  // namespace Bembel
#endif
