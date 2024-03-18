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

#ifndef BEMBEL_SRC_CLUSTERTREE_ELEMENTTREENODE_HPP_
#define BEMBEL_SRC_CLUSTERTREE_ELEMENTTREENODE_HPP_

namespace Bembel {

/**
 *  \ingroup ClusterTree
 *  \brief The ElementTreeNode corresponds to an element in the element tree.
 */
class ElementTreeNode {
 public:
  /**
   * \brief iterator struct for element tree nodes. They may be used to iterator
   * over the elements in a cluster. To do so, however, the cluster must be set
   * up by ElementTree beforehand.
   */
  struct const_iterator {
    using iterator_category = std::forward_iterator_tag;
    using difference_type = std::ptrdiff_t;
    using value_type = ElementTreeNode;
    using pointer = value_type *;
    using reference = value_type &;

    /**
     * \brief Constructs a new iterator.
     */
    explicit const_iterator(pointer ptr) : m_ptr(ptr) {}

    /**
     * \brief Accesses the pointed-to element.
     */
    reference operator*() const { return *m_ptr; }
    /**
     * \brief Accesses the pointed-to element.
     */
    const pointer operator->() const { return m_ptr; }

    /**
     * \brief Prefix increment.
     */
    const_iterator &operator++() {
      m_ptr = m_ptr->next_;
      return *this;
    }

    /**
     * \brief Postfix increment.
     */
    const_iterator operator++(int) {
      const_iterator tmp = *this;
      ++(*this);
      return tmp;
    }

    /**
     * \brief Compares the underlying iterators.
     */
    friend bool operator==(const const_iterator &a, const const_iterator &b) {
      return a.m_ptr == b.m_ptr;
    }
    /**
     * \brief Compares the underlying iterators.
     */
    friend bool operator!=(const const_iterator &a, const const_iterator &b) {
      return a.m_ptr != b.m_ptr;
    }

   private:
    pointer m_ptr;
  };
  //////////////////////////////////////////////////////////////////////////////
  /// constructors
  //////////////////////////////////////////////////////////////////////////////
  /**
   * \brief Default constructor.
   */
  ElementTreeNode() noexcept
      : prev_(nullptr),
        next_(nullptr),
        id_(-1),
        level_(-1),
        patch_(-1),
        radius_(std::numeric_limits<double>::infinity()) {
    midpoint_ << 0., 0., 0.;
    llc_ << 0., 0.;
  }
  /**
   * \brief Move constructor.
   */
  ElementTreeNode(ElementTreeNode &&other) noexcept {
    midpoint_.swap(other.midpoint_);
    llc_.swap(other.llc_);
    vertices_ = std::move(other.vertices_);
    radius_ = other.radius_;
    id_ = other.id_;
    patch_ = other.patch_;
    sons_ = std::move(other.sons_);
    adjcents_ = std::move(other.adjcents_);
  }

  /**
   * \brief Copy constructor (deleted).
   *
   * \param other The ElementTreeNode instance to copy from.
   */
  ElementTreeNode(const ElementTreeNode &other) = delete;
  /**
   * \brief Copy assignment operator (deleted).
   *
   * \param other The ElementTreeNode instance to copy from.
   * \return Reference to this ElementTreeNode instance.
   */
  ElementTreeNode &operator=(const ElementTreeNode &other) = delete;
  /**
   * \brief Move assignment operator (deleted).
   *
   * \param other The ElementTreeNode instance to move from.
   * \return Reference to this ElementTreeNode instance.
   */
  ElementTreeNode &operator=(ElementTreeNode &&other) = delete;
  //////////////////////////////////////////////////////////////////////////////
  /// methods
  //////////////////////////////////////////////////////////////////////////////
  /**
   * \brief Prints the member variables of the element.
   */
  void print() const {
    std::cout << "{" << std::endl;
    std::cout << "midpoint:   " << midpoint_.transpose() << std::endl;
    std::cout << "llc:        " << llc_.transpose() << std::endl;
    std::cout << "p s n:      " << prev_ << " " << this << " " << next_
              << std::endl;
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
  //////////////////////////////////////////////////////////////////////////////
  /**
   * \brief Maps a point in the reference domain of a patch to the reference
   * element of this element.
   *
   * This function asserts that the correct element is chosen. So the return
   * value must be with in [0,1]^2.
   *
   * \param in Point in reference domain of the patch.
   * \return Point in reference domain of the element.
   */
  Eigen::Vector2d mapToReferenceElement(const Eigen::Vector2d &in) const {
    Eigen::Vector2d out = (in - llc_) / get_h();
    assert(out(0) >= 0. && out(0) <= 1. && out(1) >= 0. && out(1) <= 1.);
    return out;
  }
  //////////////////////////////////////////////////////////////////////////////
  /**
   * \brief Get the midpoint of this element with respect to the patch.
   *
   * \return Midpoint with respect to the patch containing the element.
   */
  Eigen::Vector2d referenceMidpoint() const {
    return llc_ + Eigen::Vector2d(0.5, 0.5) * get_h();
  }
  //////////////////////////////////////////////////////////////////////////////
  /// getter
  //////////////////////////////////////////////////////////////////////////////
  /**
   * \brief Get with of the reference element with respect to the reference
   * domain of the patch.
   *
   * \return Width of the reference element.
   */
  double get_h() const { return double(1) / double(1 << level_); }
  //////////////////////////////////////////////////////////////////////////////
  /**
   * \brief Get the level of refinement of the element.
   *
   * \return Level of refinement.
   */
  int get_level() const { return level_; }
  //////////////////////////////////////////////////////////////////////////////
  /**
   * \brief Returns a const reference to the first son if any or itself.
   *
   * \return Const reference to the first son or this element itself.
   */
  const ElementTreeNode &front() const {
    const ElementTreeNode *ptr = this;
    while (ptr->sons_.size()) ptr = std::addressof(ptr->sons_.front());
    return *ptr;
  }
  //////////////////////////////////////////////////////////////////////////////
  /**
   * \brief Returns a const reference to the last son if any or itself.
   *
   * \return Const reference to the last son or this element itself.
   */
  const ElementTreeNode &back() const {
    const ElementTreeNode *ptr = this;
    while (ptr->sons_.size()) ptr = std::addressof(ptr->sons_.back());
    return *ptr;
  }
  //////////////////////////////////////////////////////////////////////////////
  /**
   * \brief Returns an iterator pointing to the element in the sequence before
   * this ElementTreeNodes.
   *
   * \return Returns a ElementTreeNode::const_iterator object.
   */

  const_iterator cbegin() const {
    const ElementTreeNode &el = this->front();
    return const_iterator(const_cast<ElementTreeNode *>(std::addressof(el)));
  }
  //////////////////////////////////////////////////////////////////////////////
  /**
   * \brief Returns an iterator pointing to the element in the sequence after
   * this ElementTreeNode.
   *
   * \return Returns a ElementTreeNode::const_iterator object.
   */
  const_iterator cend() const {
    const ElementTreeNode &el = this->back();
    return const_iterator(const_cast<ElementTreeNode *>(el.next_));
  }
  //////////////////////////////////////////////////////////////////////////////
  /**
   * \brief Returns an iterator pointing to the element in the sequence before
   * this ElementTreeNodes.
   *
   * \return Returns a ElementTreeNode::const_iterator object.
   */
  const_iterator begin() const { return cbegin(); }
  //////////////////////////////////////////////////////////////////////////////
  /**
   * \brief Returns an iterator pointing to the element in the sequence after
   * this ElementTreeNode.
   *
   * \return Returns a ElementTreeNode::const_iterator object.
   */
  const_iterator end() const {
    return cend();
  }
  //////////////////////////////////////////////////////////////////////////////
  /// member variables
  //////////////////////////////////////////////////////////////////////////////
  std::vector<ElementTreeNode> sons_;        /// children
  std::vector<ElementTreeNode *> adjcents_;  /// neighbouring elements indices
  std::vector<int> vertices_;                /// indices of the vertices
  Eigen::Vector3d midpoint_;                 /// midpoint of the element
  Eigen::Vector2d llc_;                      /// lower left corner on [0,1]^2
  ElementTreeNode *prev_;
  ElementTreeNode *next_;
  double radius_;  /// radius of the element
  int id_;         /// element id with respect to the level
  int level_;      /// level of the element
  int patch_;      /// patch of the element
};
}  // namespace Bembel
#endif  // BEMBEL_SRC_CLUSTERTREE_ELEMENTTREENODE_HPP_
