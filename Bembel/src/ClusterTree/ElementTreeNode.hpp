// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.

#ifndef BEMBEL_CLUSTERTREE_ELEMENTTREENODE_H_
#define BEMBEL_CLUSTERTREE_ELEMENTTREENODE_H_

namespace Bembel {

// forward declaration of memory is necessary here
struct ElementTreeMemory;

/**
 *  \ingroup ClusterTree
 *  \brief The ElementTreeNode correposnds to an element in the element tree.
 */
class ElementTreeNode {
 public:
  //////////////////////////////////////////////////////////////////////////////
  /// constructors
  //////////////////////////////////////////////////////////////////////////////
  ElementTreeNode()
      : id_(-1),
        level_(-1),
        patch_(-1),
        radius_(std::numeric_limits<double>::infinity()) {
    midpoint_ << 0., 0., 0.;
    llc_ << 0., 0.;
  }
  ElementTreeNode(ElementTreeNode &&other) noexcept {
    midpoint_ = std::move(other.midpoint_);
    llc_ = std::move(other.llc_);
    sons_ = std::move(other.sons_);
    vertices_ = std::move(other.vertices_);
    adjcents_ = std::move(other.adjcents_);
    id_ = other.id_;
    level_ = other.level_;
    patch_ = other.patch_;
    radius_ = other.radius_;
    memory_ = other.memory_;
  }
  ElementTreeNode(const ElementTreeNode &other) noexcept {
    midpoint_ = other.midpoint_;
    llc_ = other.llc_;
    sons_ = other.sons_;
    vertices_ = other.vertices_;
    adjcents_ = other.adjcents_;
    id_ = other.id_;
    level_ = other.level_;
    patch_ = other.patch_;
    radius_ = other.radius_;
    memory_ = other.memory_;
  }
  //////////////////////////////////////////////////////////////////////////////
  /// operators
  //////////////////////////////////////////////////////////////////////////////
  ElementTreeNode &operator=(ElementTreeNode other) noexcept {
    midpoint_.swap(other.midpoint_);
    llc_.swap(other.llc_);
    sons_.swap(other.sons_);
    vertices_.swap(other.vertices_);
    adjcents_.swap(other.adjcents_);
    id_ = other.id_;
    level_ = other.level_;
    patch_ = other.patch_;
    radius_ = other.radius_;
    memory_ = other.memory_;
    return *this;
  }
  //////////////////////////////////////////////////////////////////////////////
  /// methods
  //////////////////////////////////////////////////////////////////////////////
  void print() const {
    std::cout << "{" << std::endl;
    std::cout << "midpoint:   " << midpoint_.transpose() << std::endl;
    std::cout << "llc:        " << llc_.transpose() << std::endl;
    std::cout << "children:   ";
    for (auto i = 0; i < sons_.size(); ++i) std::cout << sons_[i] << " ";
    std::cout << std::endl;
    std::cout << "neighbours: ";
    for (auto i = 0; i < adjcents_.size(); ++i)
      std::cout << adjcents_[i] << " ";
    std::cout << std::endl;
    std::cout << "vertices:   ";
    for (auto i = 0; i < vertices_.size(); ++i)
      std::cout << vertices_[i] << " ";
    std::cout << std::endl;
    std::cout << "radius:     " << radius_ << std::endl;
    std::cout << "id:         " << id_ << std::endl;
    std::cout << "level:      " << level_ << std::endl;
    std::cout << "patch:      " << patch_ << std::endl;
    std::cout << "}" << std::endl;
    return;
  }
  Eigen::Vector2d mapToReferenceElement(const Eigen::Vector2d &in) const {
    Eigen::Vector2d out = (in - llc_) / get_h();
    assert(out(0) >= 0. && out(0) <= 1. && out(1) >= 0. && out(1) <= 1.);
    return out;
  }
  Eigen::Vector2d referenceMidpoint() const {
    return llc_+Eigen::Vector2d(0.5,0.5)*get_h();
  }
  //////////////////////////////////////////////////////////////////////////////
  /// setter
  //////////////////////////////////////////////////////////////////////////////
  void set_memory(std::shared_ptr<ElementTreeMemory> memory) {
    memory_ = memory;
    return;
  }
  //////////////////////////////////////////////////////////////////////////////
  /// getter
  //////////////////////////////////////////////////////////////////////////////
  constexpr double get_h() const { return 1. / double(1 << level_); }
  constexpr int get_level() const { return level_; }
  const std::shared_ptr<ElementTreeMemory> get_memory() const { return memory_; }
  //////////////////////////////////////////////////////////////////////////////
  /// member variables
  //////////////////////////////////////////////////////////////////////////////
  Eigen::Vector3d midpoint_;   /// midpoint of the element
  Eigen::Vector2d llc_;        /// lower left corner on [0,1]^2
  std::vector<int> vertices_;  /// indices of the vertices
  double radius_;              /// radius of the element
  int id_;                     /// element id with respect to the level
  int level_;                  /// level of the element
  int patch_;                  /// patch of the element
  std::vector<int> sons_;      /// children
  std::vector<int> adjcents_;  /// neighbouring elements indices
  std::shared_ptr<ElementTreeMemory> memory_;
};
}  // namespace Bembel
#endif
