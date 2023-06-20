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

#ifndef BEMBEL_SRC_CLUSTERTREE_ELEMENTTREE_HPP_
#define BEMBEL_SRC_CLUSTERTREE_ELEMENTTREE_HPP_

namespace Bembel {
/**
 *  \ingroup ClusterTree
 *  \brief This class organizes an element structure on a Geometry object and
 * handles refinement.
 *
 * \todo Describe the ElementTree
 */
class ElementTree {
 public:
  //////////////////////////////////////////////////////////////////////////////
  ElementTree(const ElementTree &) = delete;
  ElementTree(ElementTree &&) = delete;
  ElementTree &operator=(const ElementTree &) = delete;
  ElementTree &operator=(ElementTree &&) = delete;
  //////////////////////////////////////////////////////////////////////////////
  /// constructors
  //////////////////////////////////////////////////////////////////////////////
  ElementTree() {}
  explicit ElementTree(const Geometry &g, unsigned int max_level = 0) {
    init_ElementTree(g, max_level);
  }
  //////////////////////////////////////////////////////////////////////////////
  /// init
  //////////////////////////////////////////////////////////////////////////////
  void init_ElementTree(const Geometry &g, unsigned int max_level) {
    // initialise the data fields
    geometry_ = g.get_geometry_ptr();
    max_level_ = max_level;
    number_of_patches_ = geometry_->size();
    // create the patches and set up the topology
    {
      std::vector<Eigen::Vector3d> uniquePts;
      std::vector<ElementTreeNode *> patches;
      Eigen::Vector3d v;
      number_of_points_ = 0;
      number_of_elements_ = number_of_patches_;
      ElementTreeNode &root = root_;
      root.sons_.resize(number_of_patches_);
      for (auto i = 0; i < number_of_patches_; ++i) {
        root.sons_[i].adjcents_ = std::vector<ElementTreeNode *>(4, nullptr);
        root.sons_[i].vertices_.resize(4);
        root.sons_[i].id_ = i;
        root.sons_[i].level_ = 0;
        root.sons_[i].patch_ = i;
        // add linked list structure to the panels
        if (i == 0) {
          root.sons_[0].prev_ = nullptr;
          pfirst_ = std::addressof(root.sons_[0]);
        } else {
          root.sons_[i].prev_ = std::addressof(root.sons_[i - 1]);
        }
        if (i == number_of_patches_ - 1) {
          root.sons_[i].next_ = nullptr;
          plast_ = std::addressof(root.sons_[i]);
        } else {
          root.sons_[i].next_ = std::addressof(root.sons_[i + 1]);
        }
        // get images of the four corners of the unit square under the
        // diffeomorphism geo[i]
        for (auto j = 0; j < 4; ++j) {
          v = (*geometry_)[i].eval(Constants::corners[0][j],
                                   Constants::corners[1][j]);
          unsigned int index = 0;
          for (; index < uniquePts.size(); ++index)
            if ((uniquePts[index] - v).norm() < Constants::pt_comp_tolerance)
              break;
          if (index != uniquePts.size()) {
            root.sons_[i].vertices_[j] = index;
          } else {
            uniquePts.push_back(v);
            root.sons_[i].vertices_[j] = number_of_points_;
            ++number_of_points_;
          }
        }
        patches.push_back(std::addressof(root.sons_[i]));
      }
      updateTopology(patches);
    }
    for (auto i = 0; i < max_level; ++i) refineUniformly();
    return;
  }
  //////////////////////////////////////////////////////////////////////////////
  void refineUniformly_recursion(ElementTreeNode &el) {
    if (el.sons_.size()) {
      for (auto i = 0; i < el.sons_.size(); ++i)
        refineUniformly_recursion(el.sons_[i]);
    } else {
      refineLeaf(el);
    }
    return;
  }
  //////////////////////////////////////////////////////////////////////////////
  void refineUniformly() {
    refineUniformly_recursion(root_);
    return;
  }
  //////////////////////////////////////////////////////////////////////////////
  void refinePatch(int patch) {
    refineUniformly_recursion(root_.sons_[patch]);
    return;
  }
  //////////////////////////////////////////////////////////////////////////////
  Eigen::MatrixXd generatePointList(Eigen::VectorXi *idct = nullptr) const {
    Eigen::MatrixXd pts(3, number_of_points_);
    if (idct != nullptr) {
      idct->resize(number_of_points_);
      idct->setZero();
    }
    for (auto it = pbegin(); it != pend(); ++it) {
      double h = double(1) / double(1 << it->level_);
      pts.col(it->vertices_[0]) =
          (*geometry_)[it->patch_].eval(it->llc_(0), it->llc_(1));
      pts.col(it->vertices_[1]) =
          (*geometry_)[it->patch_].eval(it->llc_(0) + h, it->llc_(1));
      pts.col(it->vertices_[2]) =
          (*geometry_)[it->patch_].eval(it->llc_(0) + h, it->llc_(1) + h);
      pts.col(it->vertices_[3]) =
          (*geometry_)[it->patch_].eval(it->llc_(0), it->llc_(1) + h);
      if (idct != nullptr) {
        ++((*idct)(it->vertices_[0]));
        ++((*idct)(it->vertices_[1]));
        ++((*idct)(it->vertices_[2]));
        ++((*idct)(it->vertices_[3]));
      }
    }
    return pts;
  }
  //////////////////////////////////////////////////////////////////////////////
  Eigen::MatrixXi generateElementList() const {
    unsigned int elID = 0;
    Eigen::MatrixXi retval(4, number_of_elements_);
    for (auto it = pbegin(); it != pend(); ++it) {
      retval.col(elID) << it->vertices_[0], it->vertices_[1], it->vertices_[2],
          it->vertices_[3];
      ++elID;
    }
    return retval;
  }
  //////////////////////////////////////////////////////////////////////////////
  /// getter
  //////////////////////////////////////////////////////////////////////////////
  int get_number_of_points() const { return number_of_points_; }
  int get_number_of_elements() const { return number_of_elements_; }
  int get_max_level() const { return max_level_; }
  ElementTreeNode &root() { return root_; }
  const ElementTreeNode &root() const { return root_; }
  const PatchVector &get_geometry() const { return *geometry_; }
  //////////////////////////////////////////////////////////////////////////////
  /// iterators
  //////////////////////////////////////////////////////////////////////////////
  ElementTreeNode::const_iterator pbegin() const {
    return ElementTreeNode::const_iterator(pfirst_);
  }
  ElementTreeNode::const_iterator pend() const {
    return ElementTreeNode::const_iterator(plast_->next_);
  }
  ElementTreeNode::const_iterator cpbegin() const { return pbegin(); }
  ElementTreeNode::const_iterator cpend() const { return pend(); }
  ElementTreeNode::const_iterator cluster_begin(
      const ElementTreeNode &cl) const {
    return cl.cbegin();
  }
  ElementTreeNode::const_iterator cluster_end(const ElementTreeNode &cl) const {
    return cl.cend();
  }
  //////////////////////////////////////////////////////////////////////////////
  Eigen::MatrixXd computeElementEnclosings() {
    // compute point list
    Eigen::MatrixXd P = generatePointList();
    for (auto i = 0; i < root_.sons_.size(); ++i)
      computeElementEnclosings_recursion(root_.sons_[i], P);
    return P;
  }
  //////////////////////////////////////////////////////////////////////////////
  void computeElementEnclosings_recursion(ElementTreeNode &el,
                                          const Eigen::MatrixXd &P) {
    Eigen::Vector3d mp1, mp2;
    double r1, r2;
    // compute point list
    if (!el.sons_.size()) {
      // assign enclosing balls to leafs
      computeEnclosingBall(&mp1, &r1, P.col(el.vertices_[0]), 0,
                           P.col(el.vertices_[2]), 0);
      computeEnclosingBall(&mp2, &r2, P.col(el.vertices_[1]), 0,
                           P.col(el.vertices_[3]), 0);
      computeEnclosingBall(&(el.midpoint_), &(el.radius_), mp1, r1, mp2, r2);
    } else {
      // handle the four(!!!) children
      for (auto i = 0; i < 4; ++i)
        computeElementEnclosings_recursion(el.sons_[i], P);
      // assign enclosing balls to fathers bottom up
      computeEnclosingBall(&mp1, &r1, el.sons_[0].midpoint_,
                           el.sons_[0].radius_, el.sons_[2].midpoint_,
                           el.sons_[2].radius_);
      computeEnclosingBall(&mp2, &r2, el.sons_[1].midpoint_,
                           el.sons_[1].radius_, el.sons_[3].midpoint_,
                           el.sons_[3].radius_);
      computeEnclosingBall(&(el.midpoint_), &(el.radius_), mp1, r1, mp2, r2);
    }
    return;
  }
  //////////////////////////////////////////////////////////////////////////////
  void printPanels() const {
    auto i = 0;
    for (auto it = pbegin(); it != pend(); ++it) {
      std::cout << "element[" << i++ << "] = ";
      it->print();
    }
    return;
  }
  //////////////////////////////////////////////////////////////////////////////
  /// other Stuff
  //////////////////////////////////////////////////////////////////////////////
  Eigen::MatrixXd generateMidpointList() const {
    Eigen::MatrixXd retval(3, number_of_elements_);
    unsigned int i = 0;
    for (auto it = pbegin(); it != pend(); ++it) {
      retval.col(i) = it->midpoint_;
      ++i;
    }
    return retval;
  }
  //////////////////////////////////////////////////////////////////////////////
  Eigen::MatrixXd generateRadiusList() const {
    Eigen::VectorXd retval(number_of_elements_);
    unsigned int i = 0;
    for (auto it = pbegin(); it != pend(); ++it) {
      retval(i) = it->radius_;
      ++i;
    }
    return retval;
  }
  //////////////////////////////////////////////////////////////////////////////
  Eigen::VectorXi generateElementLabels() const {
    Eigen::VectorXi retval(number_of_elements_);
    unsigned int i = 0;
    for (auto it = pbegin(); it != pend(); ++it) {
      retval(i) = compute_global_id(*it);
      ++i;
    }
    return retval;
  }
  //////////////////////////////////////////////////////////////////////////////
  Eigen::VectorXi generatePatchBoundaryLabels() const {
    Eigen::VectorXi retval(number_of_elements_);
    retval.setZero();
    unsigned int i = 0;
    for (auto it = pbegin(); it != pend(); ++it) {
      for (auto j = 0; j < 4; ++j)
        if (it->adjcents_[j] == nullptr) {
          retval[i] = -1;
          break;
        } else if (it->adjcents_[j]->patch_ != it->patch_) {
          ++(retval[i]);
        }
      ++i;
    }
    return retval;
  }
  //////////////////////////////////////////////////////////////////////////////
  Eigen::VectorXi identifyPatch(unsigned int pn) const {
    Eigen::VectorXi retval(number_of_elements_);
    retval.setZero();
    assert(pn < number_of_patches_);
    unsigned int i = 0;
    for (auto it = pbegin(); it != pend(); ++it) {
      retval[i] = it->patch_ == pn ? 1 : 0;
      ++i;
    }
    return retval;
  }
  //////////////////////////////////////////////////////////////////////////////
  /// static members
  //////////////////////////////////////////////////////////////////////////////
  /**
   *  \brief computes a ball enclosing the union of B_r1(mp1) and B_r2(mp2),
   * i.e B(mp,r)\supset B_r1(mp1) \cup B_r2(mp2)
   */
  static void computeEnclosingBall(Eigen::Vector3d *mp, double *r,
                                   const Eigen::Vector3d &mp1, double r1,
                                   const Eigen::Vector3d &mp2, double r2) {
    // compute distance vector of the two spheres
    auto z = mp1 - mp2;
    auto norm = (mp1 - mp2).norm();
    /// B(d2,r2) \subset B(d1,r1)
    if (norm + r2 <= r1) {
      *mp = mp1;
      *r = r1;
      /// B(d1,r1) \subset B(d2,r2)
    } else if (norm + r1 <= r2) {
      *mp = mp2;
      *r = r2;
      /// the union is not a ball
    } else {
      *mp = 0.5 * (mp1 + mp2 + (r1 - r2) / norm * z);
      *r = 0.5 * (r1 + r2 + norm);
      *r = 0.5 * (r1 + r2 + norm);
    }
    return;
  }
  /**
   *  \brief
   */
  //////////////////////////////////////////////////////////////////////////////
  std::vector<std::array<int, 4>> patchTopologyInfo() const {
    std::vector<std::array<int, 4>> retval;
    for (auto it = root_.sons_.begin(); it != root_.sons_.end(); ++it) {
      for (auto j = 0; j < 4; ++j) {
        // do we have a neighbour?
        if (it->adjcents_[j] != nullptr) {
          const ElementTreeNode &cur_neighbour = *(it->adjcents_[j]);
          // add the edge only if it->id_ < neighbour->id_
          // otherwise the edge has already been added
          if (it->id_ < cur_neighbour.id_) {
            int k = 0;
            for (; k < 4; ++k)
              if (cur_neighbour.adjcents_[k] != nullptr &&
                  cur_neighbour.adjcents_[k]->id_ == it->id_)
                break;
            retval.push_back({it->id_, cur_neighbour.id_, j, k});
          }
        } else {
          retval.push_back({it->id_, -1, j, -1});
        }
      }
    }
    return retval;
  }
  // The ordering of elements in the element tree does not correspond to the
  // element order underlying the coefficient vector. This reordering can be
  // computed for look ups by this function.
  // TODO(Max) This function assumes that everything is refined uniformly!
  //////////////////////////////////////////////////////////////////////////////
  std::vector<int> computeReorderingVector() const {
    std::vector<int> out(number_of_elements_);
    for (auto it = pbegin(); it != pend(); ++it) {
      const double h = it->get_h();
      const Eigen::Vector2d ref_mid =
          it->llc_ + Eigen::Vector2d(.5 * h, .5 * h);
      const int x_idx = std::floor(ref_mid(0) / h);
      const int y_idx = std::floor(ref_mid(1) / h);
      const int num_in_one_dir = 1 << max_level_;
      const int num_per_patch = num_in_one_dir * num_in_one_dir;
      const int tp_idx =
          it->patch_ * num_per_patch + y_idx * num_in_one_dir + x_idx;
      out[tp_idx] = it->id_;
    }
    return out;
  }
  //////////////////////////////////////////////////////////////////////////////
  size_t compute_global_id(const ElementTreeNode &el) const {
    return number_of_patches_ *
               (((1 << el.level_) * (1 << el.level_) - 1) / 3) +
           el.id_;
  }
  //////////////////////////////////////////////////////////////////////////////
  /// private members
  //////////////////////////////////////////////////////////////////////////////
 private:
  //////////////////////////////////////////////////////////////////////////////
  /// methods
  //////////////////////////////////////////////////////////////////////////////
  template <typename T>
  struct isEqual {
    bool operator()(const T &v1, const T &v2) const {
      return ((v1 - v2).norm() < Constants::pt_comp_tolerance);
    }
  };
  /**
   * \brief function to set up the local topology, i.e. the adjacents_
   *        of a refined element
   **/
  //////////////////////////////////////////////////////////////////////////////
  void updateTopology(const std::vector<ElementTreeNode *> &elements) {
    std::map<std::array<int, 2>, ElementTreeNode *> edges;
    std::array<int, 2> e1, e2;
    // generate edge list for all elements in question
    for (auto i = 0; i < elements.size(); ++i)
      for (auto j = 0; j < 4; ++j) {
        // compute a unique id for each edge
        const int v1 = elements[i]->vertices_[j];
        const int v2 = elements[i]->vertices_[(j + 1) % 4];
        e1 = {v1 < v2 ? v1 : v2, v1 < v2 ? v2 : v1};
        // perferm a look up if the edge is already existing.
        auto it = edges.find(e1);
        // if so, identify the two neighbours
        if (it != edges.end()) {
          elements[i]->adjcents_[j] = it->second;
          // now find the edge also in the patch that added it to the set
          for (auto k = 0; k < 4; ++k) {
            const int v3 = it->second->vertices_[k];
            const int v4 = it->second->vertices_[(k + 1) % 4];
            e2 = {v3 < v4 ? v3 : v4, v3 < v4 ? v4 : v3};
            if (e1 == e2) {
              it->second->adjcents_[k] = elements[i];
              break;
            }
          }
          // otherwise add the edge to the list
        } else {
          edges.insert(std::make_pair(e1, elements[i]));
        }
      }
    return;
  }
  //////////////////////////////////////////////////////////////////////////////
  void refineLeaf(ElementTreeNode &cur_el) {
    // check if we have actually a panel
    if (cur_el.sons_.size()) return;
    std::vector<int> ptIds(5);
    std::vector<int> refNeighbours(4);
    std::vector<ElementTreeNode *> elements;
    // determine position of p with respect to its neighbours
    for (auto i = 0; i < 4; ++i) {
      // is there a neighbour?
      if (cur_el.adjcents_[i] != nullptr) {
        for (auto j = 0; j < 4; ++j)
          if (cur_el.adjcents_[i]->adjcents_[j] == std::addressof(cur_el)) {
            refNeighbours[i] = j;
            break;
          }
      } else {
        refNeighbours[i] = -1;
      }
    }
    //  determine new points
    for (auto i = 0; i < 4; ++i) {
      // is there a neighbour?
      if (cur_el.adjcents_[i] != nullptr) {
        ElementTreeNode &ref_cur_neighbour = *(cur_el.adjcents_[i]);
        // is the neighbour already refined?
        if (ref_cur_neighbour.sons_.size()) {
          // this is the midpoint of the shared edge
          ptIds[i] = ref_cur_neighbour.sons_[refNeighbours[i]]
                         .vertices_[(refNeighbours[i] + 1) % 4];
          // these are the two elements adjacent to the edge
          elements.push_back(
              std::addressof(ref_cur_neighbour.sons_[refNeighbours[i]]));
          elements.push_back(std::addressof(
              ref_cur_neighbour.sons_[(refNeighbours[i] + 1) % 4]));
        } else {
          // otherwise add the point id to the tree
          ptIds[i] = number_of_points_;
          ++number_of_points_;
        }
      } else {
        // otherwise add the point id to the tree
        ptIds[i] = number_of_points_;
        ++number_of_points_;
      }
    }
    // add midpoint of the current element, which is always a new point
    ptIds[4] = number_of_points_;
    ++number_of_points_;
    // set up new elements
    cur_el.sons_.resize(4);
    number_of_elements_ += 3;
    max_level_ =
        max_level_ < cur_el.level_ + 1 ? cur_el.level_ + 1 : max_level_;
    //  auto it = leafs_.find(compute_global_id(cur_el));
    //  if (it != leafs_.end()) leafs_.erase(it);
    for (auto i = 0; i < 4; ++i) {
      // add linked list structure to the panels
      //////////////////////////////////////////////////////////////////////////
      if (i == 0) {
        cur_el.sons_[i].prev_ = cur_el.prev_;
        if (cur_el.prev_ != nullptr)
          (cur_el.prev_)->next_ = std::addressof(cur_el.sons_[i]);
        cur_el.prev_ = nullptr;
        if (std::addressof(cur_el) == pfirst_)
          pfirst_ = std::addressof(cur_el.sons_[i]);

      } else {
        cur_el.sons_[i].prev_ = std::addressof(cur_el.sons_[i - 1]);
      }
      if (i == 3) {
        cur_el.sons_[i].next_ = cur_el.next_;
        if (cur_el.next_ != nullptr)
          (cur_el.next_)->prev_ = std::addressof(cur_el.sons_[i]);
        cur_el.next_ = nullptr;
        if (std::addressof(cur_el) == plast_)
          plast_ = std::addressof(cur_el.sons_[i]);
      } else {
        cur_el.sons_[i].next_ = std::addressof(cur_el.sons_[i + 1]);
      }
      //////////////////////////////////////////////////////////////////////////
      cur_el.sons_[i].patch_ = cur_el.patch_;
      cur_el.sons_[i].level_ = cur_el.level_ + 1;
      cur_el.sons_[i].id_ = 4 * cur_el.id_ + i;
      cur_el.sons_[i].adjcents_ = std::vector<ElementTreeNode *>(4, nullptr);
      cur_el.sons_[i].llc_(0) =
          cur_el.llc_(0) + Constants::llcs[0][i] / double(1 << cur_el.level_);
      cur_el.sons_[i].llc_(1) =
          cur_el.llc_(1) + Constants::llcs[1][i] / double(1 << cur_el.level_);
      elements.push_back(std::addressof(cur_el.sons_[i]));
    }
    // set vertices
    // first element
    cur_el.sons_[0].vertices_.push_back(cur_el.vertices_[0]);
    cur_el.sons_[0].vertices_.push_back(ptIds[0]);
    cur_el.sons_[0].vertices_.push_back(ptIds[4]);
    cur_el.sons_[0].vertices_.push_back(ptIds[3]);
    // second element
    cur_el.sons_[1].vertices_.push_back(ptIds[0]);
    cur_el.sons_[1].vertices_.push_back(cur_el.vertices_[1]);
    cur_el.sons_[1].vertices_.push_back(ptIds[1]);
    cur_el.sons_[1].vertices_.push_back(ptIds[4]);
    // third element
    cur_el.sons_[2].vertices_.push_back(ptIds[4]);
    cur_el.sons_[2].vertices_.push_back(ptIds[1]);
    cur_el.sons_[2].vertices_.push_back(cur_el.vertices_[2]);
    cur_el.sons_[2].vertices_.push_back(ptIds[2]);
    // fourth element
    cur_el.sons_[3].vertices_.push_back(ptIds[3]);
    cur_el.sons_[3].vertices_.push_back(ptIds[4]);
    cur_el.sons_[3].vertices_.push_back(ptIds[2]);
    cur_el.sons_[3].vertices_.push_back(cur_el.vertices_[3]);
    // fix adjecency relations
    updateTopology(elements);

    return;
  }
  //////////////////////////////////////////////////////////////////////////////
  /// member variables
  //////////////////////////////////////////////////////////////////////////////
  std::shared_ptr<PatchVector> geometry_;
  ElementTreeNode root_;
  ElementTreeNode *pfirst_;
  ElementTreeNode *plast_;
  int number_of_patches_;
  int max_level_;
  int number_of_points_;
  int number_of_elements_;
  //////////////////////////////////////////////////////////////////////////////
};
}  // namespace Bembel
#endif  // BEMBEL_SRC_CLUSTERTREE_ELEMENTTREE_HPP_
