// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.

#ifndef BEMBEL_CLUSTERTREE_ELEMENTTREEMEMORY_H_
#define BEMBEL_CLUSTERTREE_ELEMENTTREEMEMORY_H_

namespace Bembel {

/**
 *  \ingroup ClusterTree
 *  \brief this nice struct keeps track of the entire memory management of the
 *         element tree. In fact, due to this struct, the tree is now easily
 *         copyable.
 */
struct ElementTreeMemory {
  //////////////////////////////////////////////////////////////////////////////
  std::vector<int>::size_type nsons(ElementTreeNode *etn) {
    return etn->sons_.size();
  }

  ElementTreeNode &son(ElementTreeNode &etn, std::vector<int>::size_type id) {
    return (*memory_)[etn.sons_[id]];
  }

  const ElementTreeNode &son(const ElementTreeNode &etn,
                             std::vector<int>::size_type id) const {
    return (*memory_)[etn.sons_[id]];
  }

  ElementTreeNode &adjcent(ElementTreeNode &etn,
                           std::vector<int>::size_type id) {
    return (*memory_)[etn.adjcents_[id]];
  }

  const ElementTreeNode &adjcent(const ElementTreeNode &etn,
                                 std::vector<int>::size_type id) const {
    return (*memory_)[etn.adjcents_[id]];
  }

  int cumNumElements(int l) const {
    return l == -1 ? 1 : number_of_patches_ * ((1 << (2 * l + 2)) - 1) / 3 + 1;
  }

  ElementTreeNode &get_root() { return (*memory_)[0]; }

  const ElementTreeNode &get_root() const { return (*memory_)[0]; }

  ElementTreeNode &get_element(std::vector<ElementTreeNode>::size_type id) {
    return (*memory_)[id];
  }

  const ElementTreeNode &get_element(
      std::vector<ElementTreeNode>::size_type id) const {
    return (*memory_)[id];
  }

  //////////////////////////////////////////////////////////////////////////////
  /// iterators
  //////////////////////////////////////////////////////////////////////////////
  std::vector<ElementTreeNode>::const_iterator cluster_begin(
      const ElementTreeNode &etn) const {
    const ElementTreeNode *left = std::addressof(etn);
    while (left->sons_.size()) left = std::addressof(son(*left, 0));
    assert(left->level_ == max_level_ && "panels on different levels");
    size_t inc = cumNumElements(max_level_ - 1) + left->id_;
    return (*memory_).begin() + inc;
  }

  std::vector<ElementTreeNode>::const_iterator cluster_end(
      const ElementTreeNode &etn) const {
    const ElementTreeNode *right = std::addressof(etn);
    while (right->sons_.size())
      right = std::addressof(son(*right, right->sons_.size() - 1));
    assert(right->level_ == max_level_ && "panels on different levels");
    size_t inc = cumNumElements(max_level_ - 1) + right->id_ + 1;
    return (*memory_).begin() + inc;
  }

  std::vector<ElementTreeNode>::const_iterator cpbegin() const {
    size_t inc = cumNumElements(max_level_ - 1);
    return (*memory_).cbegin() + inc;
  }

  std::vector<ElementTreeNode>::const_iterator cpend() const {
    return (*memory_).cend();
  }

  std::vector<ElementTreeNode>::iterator pbegin() {
    size_t inc = cumNumElements(max_level_ - 1);
    return (*memory_).begin() + inc;
  }

  std::vector<ElementTreeNode>::const_iterator pend() {
    return (*memory_).end();
  }

  std::vector<ElementTreeNode>::iterator lbegin(unsigned int level) {
    size_t inc = cumNumElements(int(level) - 1);
    return (*memory_).begin() + inc;
  }

  std::vector<ElementTreeNode>::const_iterator lend(unsigned int level) {
    size_t inc = cumNumElements(int(level));
    return (*memory_).begin() + inc;
  }

  std::vector<ElementTreeNode>::const_iterator clbegin(
      unsigned int level) const {
    size_t inc = cumNumElements(int(level) - 1);
    return (*memory_).cbegin() + inc;
  }

  std::vector<ElementTreeNode>::const_iterator clend(unsigned int level) const {
    size_t inc = cumNumElements(int(level));
    return (*memory_).cbegin() + inc;
  }
  //////////////////////////////////////////////////////////////////////////////
  /// member variables
  //////////////////////////////////////////////////////////////////////////////
  std::shared_ptr<std::vector<ElementTreeNode>> memory_;
  int number_of_patches_;
  int max_level_;
};
}  // namespace Bembel
#endif
