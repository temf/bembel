// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.

#ifndef BEMBEL_CLUSTERTREE_ELEMENTTREE_H_
#define BEMBEL_CLUSTERTREE_ELEMENTTREE_H_

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
  ElementTree(const Geometry &g, unsigned int max_level = 0) {
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
        // get images of the four corners of the unit square under the
        // diffeomorphism geo[i]
        for (auto j = 0; j < 4; ++j) {
          v = (*geometry_)[i].eval(Constants::corners[0][j],
                                   Constants::corners[1][j]);
          unsigned int index = 0;
          for (; index < uniquePts.size(); ++index)
            if ((uniquePts[index] - v).norm() < Constants::pt_comp_tolerance)
              break;
          if (index != uniquePts.size())
            root.sons_[i].vertices_[j] = index;
          else {
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
  void refineUniformly_recursion(ElementTreeNode &el) {
    if (el.sons_.size())
      for (auto i = 0; i < el.sons_.size(); ++i)
        refineUniformly_recursion(el.sons_[i]);
    else
      refineLeaf(el);
    return;
  }
  void refineUniformly() {
    refineUniformly_recursion(root_);
    return;
  }
  //////////////////////////////////////////////////////////////////////////////
  void generatePointList_recursion(Eigen::MatrixXd *pts,
                                   const ElementTreeNode &el) const {
    if (el.sons_.size())
      for (auto i = 0; i < el.sons_.size(); ++i)
        generatePointList_recursion(pts, el.sons_[i]);
    else {
      double h = double(1) / double(1 << el.level_);
      pts->col(el.vertices_[0]) =
          (*geometry_)[el.patch_].eval(el.llc_(0), el.llc_(1));
      pts->col(el.vertices_[1]) =
          (*geometry_)[el.patch_].eval(el.llc_(0) + h, el.llc_(1));
      pts->col(el.vertices_[2]) =
          (*geometry_)[el.patch_].eval(el.llc_(0) + h, el.llc_(1) + h);
      pts->col(el.vertices_[3]) =
          (*geometry_)[el.patch_].eval(el.llc_(0), el.llc_(1) + h);
    }
    return;
  }
  //////////////////////////////////////////////////////////////////////////////
  Eigen::MatrixXd generatePointList() const {
    Eigen::MatrixXd retval(3, number_of_points_);
    generatePointList_recursion(&retval, root_);
    return retval;
  }
  //////////////////////////////////////////////////////////////////////////////
  void generateElementList_recursion(Eigen::MatrixXi *elements,
                                     const ElementTreeNode &el,
                                     unsigned int &elID) const {
    if (el.sons_.size())
      for (auto i = 0; i < el.sons_.size(); ++i)
        generateElementList_recursion(elements, el.sons_[i], elID);
    else {
      elements->col(elID) << el.vertices_[0], el.vertices_[1], el.vertices_[2],
          el.vertices_[3];
      ++elID;
    }
    return;
  }
  //////////////////////////////////////////////////////////////////////////////
  Eigen::MatrixXi generateElementList() const {
    unsigned int elID = 0;
    Eigen::MatrixXi retval(4, number_of_elements_);
    generateElementList_recursion(&retval, root_, elID);
    return retval;
  }
#if 0
  //////////////////////////////////////////////////////////////////////////////
  Eigen::MatrixXd computeElementEnclosings() {
    // compute point list
    Eigen::MatrixXd P = generatePointList();
    // assign enclosing balls to leafs
    for (auto it = mem_->pbegin(); it != mem_->pend(); ++it) {
      Eigen::Vector3d mp1, mp2;
      double r1, r2;
      computeEnclosingBall(&mp1, &r1, P.col(it->vertices_[0]), 0,
                           P.col(it->vertices_[2]), 0);
      computeEnclosingBall(&mp2, &r2, P.col(it->vertices_[1]), 0,
                           P.col(it->vertices_[3]), 0);
      computeEnclosingBall(&(it->midpoint_), &(it->radius_), mp1, r1, mp2, r2);
    }
    // assign enclosing balls to fathers bottom up
    for (auto level = max_level_ - 1; level >= 0; --level) {
      for (auto it = mem_->lbegin(level); it != mem_->lend(level); ++it) {
        Eigen::Vector3d mp1, mp2;
        double r1, r2;
        computeEnclosingBall(
            &mp1, &r1, mem_->son(*it, 0).midpoint_, mem_->son(*it, 0).radius_,
            mem_->son(*it, 2).midpoint_, mem_->son(*it, 2).radius_);
        computeEnclosingBall(
            &mp2, &r2, mem_->son(*it, 1).midpoint_, mem_->son(*it, 1).radius_,
            mem_->son(*it, 3).midpoint_, mem_->son(*it, 3).radius_);
        computeEnclosingBall(&(it->midpoint_), &(it->radius_), mp1, r1, mp2,
                             r2);
      }
    }
    return P;
  }
  //////////////////////////////////////////////////////////////////////////////
  void print() const {
    auto i = 0;
    for (auto it = mem_->begin(); it != mem_->end(); ++it) {
      std::cout << "element[" << i++ << "] = ";
      it->print();
    }
    return;
  }
  //////////////////////////////////////////////////////////////////////////////
  /// getter
  //////////////////////////////////////////////////////////////////////////////
  int get_number_of_points() const { return number_of_points_; }
  int get_number_of_elements() const { return mem_->cpend() - mem_->cpbegin(); }
  ElementTreeNode &root() { return mem_->get_root(); }
  const ElementTreeNode &root() const { return mem_->get_root(); }
  //////////////////////////////////////////////////////////////////////////////
  /// iterators
  //////////////////////////////////////////////////////////////////////////////
  std::vector<ElementTreeNode>::const_iterator cpbegin() const {
    return mem_->cpbegin();
  }
  std::vector<ElementTreeNode>::const_iterator cpend() const {
    return mem_->cpend();
  }
  std::vector<ElementTreeNode>::iterator pbegin() { return mem_->pbegin(); }
  std::vector<ElementTreeNode>::const_iterator pend() { return mem_->pend(); }
  std::vector<ElementTreeNode>::iterator lbegin(unsigned int level) {
    return mem_->lbegin(level);
  }
  std::vector<ElementTreeNode>::const_iterator lend(unsigned int level) {
    return mem_->lend(level);
  }
  std::vector<ElementTreeNode>::const_iterator clbegin(unsigned int level) {
    return mem_->clbegin(level);
  }
  std::vector<ElementTreeNode>::const_iterator clend(unsigned int level) const {
    return mem_->clend(level);
  }
  //////////////////////////////////////////////////////////////////////////////
  /// other Stuff
  //////////////////////////////////////////////////////////////////////////////
  Eigen::MatrixXd generateMidpointList() const {
    Eigen::MatrixXd retval(3, mem_->cumNumElements(mem_->get_max_level()));
    unsigned int i = 0;
    for (auto it = mem_->begin(); it != mem_->end(); ++it)
      retval.col(i++) = it->midpoint_;
    return retval;
  }
  //////////////////////////////////////////////////////////////////////////////
  Eigen::MatrixXd generateRadiusList() const {
    Eigen::VectorXd retval(mem_->cumNumElements(mem_->get_max_level()));
    unsigned int i = 0;
    for (auto it = mem_->begin(); it != mem_->end(); ++it)
      retval(i++) = it->radius_;
    return retval;
  }

  //////////////////////////////////////////////////////////////////////////////
  Eigen::MatrixXi generateElementList() const {
    Eigen::MatrixXi retval(4, number_of_patches_ * (1 << 2 * max_level_));
    unsigned int i = 0;
    for (auto it = mem_->cpbegin(); it != mem_->cpend(); ++it) {
      retval.col(i++) << it->vertices_[0], it->vertices_[1], it->vertices_[2],
          it->vertices_[3];
    }
    return retval;
  }
  //////////////////////////////////////////////////////////////////////////////
  Eigen::VectorXi generateElementLabels() const {
    Eigen::VectorXi retval(number_of_patches_ * (1 << 2 * max_level_));
    unsigned int i = 0;
    for (auto it = mem_->cpbegin(); it != mem_->cpend(); ++it) {
      retval(i++) = it->id_;
    }
    return retval;
  }
  //////////////////////////////////////////////////////////////////////////////
  Eigen::VectorXi generatePatchBoundaryLabels() const {
    Eigen::VectorXi retval(number_of_patches_ * (1 << 2 * max_level_));
    retval.setZero();
    unsigned int i = 0;
    for (auto it = mem_->cpbegin(); it != mem_->cpend(); ++it) {
      for (auto j = 0; j < 4; ++j)
        if (it->adjcents_[j] == -1) {
          retval[i] = -1;
          break;
        } else if (mem_->adjcent(*it, j).patch_ != it->patch_)
          ++(retval[i]);
      ++i;
    }
    return retval;
  }
  //////////////////////////////////////////////////////////////////////////////
  Eigen::VectorXi identifyPatch(unsigned int pn) const {
    Eigen::VectorXi retval(number_of_patches_ * (1 << 2 * max_level_));
    retval.setZero();
    assert(pn < number_of_patches_);
    ElementTreeNode &patch = *(mem_->begin() + 1 + pn);
    auto begin = mem_->cluster_begin(patch);
    auto end = mem_->cluster_end(patch);
    for (auto it = begin; it != end; ++it) {
      retval[it->id_] = 1;
    }
    return retval;
  }
  //////////////////////////////////////////////////////////////////////////////
  /// static members
  //////////////////////////////////////////////////////////////////////////////
  /**
   *  \brief computes a ball enclosing the union of B_r1(mp1) and B_r2(mp2), i.e
   *         B(mp,r)\supset B_r1(mp1) \cup B_r2(mp2)
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
    for (auto it = mem_->clbegin(0); it != mem_->clend(0); ++it) {
      for (auto j = 0; j < 4; ++j) {
        // do we have a neighbour?
        if (it->adjcents_[j] != -1) {
          const ElementTreeNode &cur_neighbour = mem_->adjcent(*it, j);
          // add the edge only if it->id_ < neighbour->id_
          // otherwise the edge has already been added
          if (it->id_ < cur_neighbour.id_) {
            int k = 0;
            for (; k < 4; ++k)
              if (cur_neighbour.adjcents_[k] != -1 &&
                  mem_->adjcent(cur_neighbour, k).id_ == it->id_)
                break;
            retval.push_back({it->id_, cur_neighbour.id_, j, k});
          }
        } else
          retval.push_back({it->id_, -1, j, -1});
      }
    }
    return retval;
  }
  // The ordering of elements in the element tree does not correspond to the
  // element order underlying the coefficient vector. This reordering can be
  // computet for look ups by this function.
  //////////////////////////////////////////////////////////////////////////////
  std::vector<int> computeReorderingVector() const {
    std::vector<int> out(get_number_of_elements());
    for (auto it = mem_->cpbegin(); it != mem_->cpend(); ++it) {
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
#endif
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
        } else
          edges.insert(std::make_pair(e1, elements[i]));
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
      } else
        refNeighbours[i] = -1;
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
        }
        // otherwise add the point id to the tree
        else {
          ptIds[i] = number_of_points_;
          ++number_of_points_;
        }
      }
      // otherwise add the point id to the tree
      else {
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
    for (auto i = 0; i < 4; ++i) {
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
  int number_of_patches_;
  int max_level_;
  int number_of_points_;
  int number_of_elements_;
  //////////////////////////////////////////////////////////////////////////////
};
}  // namespace Bembel
#endif
