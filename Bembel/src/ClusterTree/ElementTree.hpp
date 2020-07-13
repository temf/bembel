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
  /// constructors
  //////////////////////////////////////////////////////////////////////////////
  ElementTree() {}
  ElementTree(const Geometry &g, unsigned int max_level) {
    init_ElementTree(g, max_level);
  }
  //////////////////////////////////////////////////////////////////////////////
  /// init
  //////////////////////////////////////////////////////////////////////////////
  void init_ElementTree(const Geometry &g, unsigned int max_level) {
    // initialise the data fields
    geometry_ = g.get_geometry_ptr();
    mem_ = std::make_shared<ElementTreeMemory>();
    mem_->number_of_patches_ = geometry_->size();
    mem_->max_level_ = max_level;
    max_level_ = max_level;
    number_of_patches_ = mem_->number_of_patches_;
    // balanced quadTree > the number of elements per branch is
    // \sum_{l=0}^L 4^l = (4^{L+1}-1) / 3;
    mem_->memory_ = std::make_shared<std::vector<ElementTreeNode>>();
    mem_->memory_->resize(mem_->cumNumElements(max_level));
    // create the patches and set up the topology
    {
      std::vector<Eigen::Vector3d> uniquePts;
      std::vector<ElementTreeNode *> patches;
      Eigen::Vector3d v;
      number_of_points_ = 0;
      // set hold all element
      ElementTreeNode &root = mem_->get_root();
      root.sons_.resize(mem_->number_of_patches_);
      root.set_memory(mem_);

      for (auto i = 0; i < number_of_patches_; ++i) {
        root.sons_[i] = i + 1;
        mem_->son(root, i).adjcents_ = std::vector<int>(4, -1);
        mem_->son(root, i).vertices_.resize(4);
        mem_->son(root, i).id_ = i;
        mem_->son(root, i).level_ = 0;
        mem_->son(root, i).patch_ = i;
        mem_->son(root, i).set_memory(mem_);
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
            mem_->son(root, i).vertices_[j] = index;
          else {
            uniquePts.push_back(v);
            mem_->son(root, i).vertices_[j] = number_of_points_;
            ++number_of_points_;
          }
        }
        patches.push_back(std::addressof(mem_->son(root, i)));
      }
      updateTopology(patches);
    }
    // refine mesh
    {
      for (auto level = 0; level < max_level_; ++level) {
        auto offset1 = mem_->cumNumElements(level - 1);
        auto offset2 = mem_->cumNumElements(level);
        for (auto j = offset1; j < offset2; ++j) refineLeaf(j);
      }
    }
    return;
  }

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
  void print() const {
    auto i = 0;
    for (auto it = mem_->memory_->begin(); it != mem_->memory_->end(); ++it) {
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
    Eigen::MatrixXd retval(3, mem_->cumNumElements(mem_->max_level_));
    unsigned int i = 0;
    for (auto it = mem_->memory_->begin(); it != mem_->memory_->end(); ++it)
      retval.col(i++) = it->midpoint_;
    return retval;
  }
  Eigen::MatrixXd generateRadiusList() const {
    Eigen::VectorXd retval(mem_->cumNumElements(mem_->max_level_));
    unsigned int i = 0;
    for (auto it = mem_->memory_->begin(); it != mem_->memory_->end(); ++it)
      retval(i++) = it->radius_;
    return retval;
  }

  Eigen::MatrixXd generatePointList() const {
    Eigen::MatrixXd retval(3, number_of_points_);
    double h = mem_->cpbegin()->get_h();
    for (auto it = mem_->cpbegin(); it != mem_->cpend(); ++it) {
      retval.col(it->vertices_[0]) =
          (*geometry_)[it->patch_].eval(it->llc_(0), it->llc_(1));
      retval.col(it->vertices_[1]) =
          (*geometry_)[it->patch_].eval(it->llc_(0) + h, it->llc_(1));
      retval.col(it->vertices_[2]) =
          (*geometry_)[it->patch_].eval(it->llc_(0) + h, it->llc_(1) + h);
      retval.col(it->vertices_[3]) =
          (*geometry_)[it->patch_].eval(it->llc_(0), it->llc_(1) + h);
    }
    return retval;
  }
  Eigen::MatrixXi generateElementList() const {
    Eigen::MatrixXi retval(4, number_of_patches_ * (1 << 2 * max_level_));
    unsigned int i = 0;
    for (auto it = mem_->cpbegin(); it != mem_->cpend(); ++it) {
      retval.col(i++) << it->vertices_[0], it->vertices_[1], it->vertices_[2],
          it->vertices_[3];
    }
    return retval;
  }
  Eigen::VectorXi generateElementLabels() const {
    Eigen::VectorXi retval(number_of_patches_ * (1 << 2 * max_level_));
    unsigned int i = 0;
    for (auto it = mem_->cpbegin(); it != mem_->cpend(); ++it) {
      retval(i++) = it->id_;
    }
    return retval;
  }
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
  Eigen::VectorXi identifyPatch(unsigned int pn) const {
    Eigen::VectorXi retval(number_of_patches_ * (1 << 2 * max_level_));
    retval.setZero();
    assert(pn < number_of_patches_);
    ElementTreeNode &patch = (*(mem_->memory_))[1 + pn];
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
  //////////////////////////////////////////////////////////////////////////////
  /// private members
  //////////////////////////////////////////////////////////////////////////////
 private:
  //////////////////////////////////////////////////////////////////////////////
  /// member variables
  //////////////////////////////////////////////////////////////////////////////
  std::shared_ptr<PatchVector> geometry_;
  std::shared_ptr<ElementTreeMemory> mem_;
  int number_of_patches_;
  int max_level_;
  int number_of_points_;
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
   *  \brief comparison functor for edges
   */
  template <typename T>
  struct edgeComp {
    bool operator()(const T &v1, const T &v2) const {
      return std::lexicographical_compare(v1.data(), v1.data() + 2, v2.data(),
                                          v2.data() + 2);
    }
  };
  /**
   * \brief function to set up the local topolog, i.e. the adjacents_
   *        y of a refined element
   **/
  void updateTopology(const std::vector<ElementTreeNode *> &elements) {
    std::set<std::array<int, 3>, edgeComp<std::array<int, 3>>> edges;
    std::array<int, 3> e1, e2;
    // generate edge list for all elements in question
    for (auto i = 0; i < elements.size(); ++i) {
      // iterate over the four edges of each element
      for (auto j = 0; j < 4; ++j) {
        // compute a unique id for each edge
        if (elements[i]->vertices_[j] < elements[i]->vertices_[(j + 1) % 4]) {
          e1[0] = elements[i]->vertices_[j];
          e1[1] = elements[i]->vertices_[(j + 1) % 4];
        } else {
          e1[1] = elements[i]->vertices_[j];
          e1[0] = elements[i]->vertices_[(j + 1) % 4];
        }
        e1[2] = i;
        // perferm a look up if the edge is already existing.
        auto it = edges.find(e1);
        // if so, identify the two neighbours
        if (it != edges.end()) {
          elements[i]->adjcents_[j] =
              mem_->cumNumElements(elements[(*it)[2]]->level_ - 1) +
              elements[(*it)[2]]->id_;
          // now find the edge also in the patch that added it to the set
          for (auto k = 0; k < 4; ++k) {
            if (elements[(*it)[2]]->vertices_[k] <
                elements[(*it)[2]]->vertices_[(k + 1) % 4]) {
              e2[0] = elements[(*it)[2]]->vertices_[k];
              e2[1] = elements[(*it)[2]]->vertices_[(k + 1) % 4];
            } else {
              e2[1] = elements[(*it)[2]]->vertices_[k];
              e2[0] = elements[(*it)[2]]->vertices_[(k + 1) % 4];
            }
            if (e1[0] == e2[0] && e1[1] == e2[1])
              elements[(*it)[2]]->adjcents_[k] =
                  mem_->cumNumElements(elements[i]->level_ - 1) +
                  elements[i]->id_;
          }
          // otherwise add the edge to the list
        } else
          edges.insert(e1);
      }
    }
    return;
  }
  void refineLeaf(int global_id) {
    // check if we have actually a panel
    ElementTreeNode &cur_el = mem_->get_element(global_id);
    if (cur_el.sons_.size()) return;
    std::vector<int> ptIds(5);
    std::vector<int> refNeighbours(4);
    std::vector<ElementTreeNode *> elements;
    // determine position of p with respect to its neighbours
    for (auto i = 0; i < 4; ++i) {
      // is there a neighbour?
      if (cur_el.adjcents_[i] != -1) {
        for (auto j = 0; j < 4; ++j)
          if (mem_->adjcent(cur_el, i).adjcents_[j] == global_id) {
            refNeighbours[i] = j;
            break;
          }
      } else
        refNeighbours[i] = -1;
    }
    //  determine new points
    for (auto i = 0; i < 4; ++i) {
      auto cur_neighbour = cur_el.adjcents_[i];
      if (cur_neighbour != -1) {
        ElementTreeNode &ref_cur_neighbour = mem_->adjcent(cur_el, i);
        // is the neighbour already refined?
        if (ref_cur_neighbour.sons_.size()) {
          // this is the midpoint of the shared edge
          ptIds[i] = mem_->son(ref_cur_neighbour, refNeighbours[i])
                         .vertices_[(refNeighbours[i] + 1) % 4];
          // these are the two elements adjacent to the edge
          elements.push_back(
              std::addressof(mem_->son(ref_cur_neighbour, refNeighbours[i])));
          elements.push_back(std::addressof(
              mem_->son(ref_cur_neighbour, (refNeighbours[i] + 1) % 4)));

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
    for (auto i = 0; i < 4; ++i) {
      auto child_ind = 4 * cur_el.id_ + i + mem_->cumNumElements(cur_el.level_);
      cur_el.sons_[i] = child_ind;
      mem_->son(cur_el, i).patch_ = cur_el.patch_;
      mem_->son(cur_el, i).level_ = cur_el.level_ + 1;
      mem_->son(cur_el, i).id_ = 4 * cur_el.id_ + i;
      mem_->son(cur_el, i).adjcents_ = std::vector<int>(4, -1);
      mem_->son(cur_el, i).llc_(0) =
          cur_el.llc_(0) +
          double(Constants::llcs[0][i]) / double(1 << cur_el.level_);
      mem_->son(cur_el, i).llc_(1) =
          cur_el.llc_(1) +
          double(Constants::llcs[1][i]) / double(1 << cur_el.level_);
      mem_->son(cur_el, i).set_memory(mem_);
      elements.push_back(std::addressof(mem_->son(cur_el, i)));
    }
    // set vertices
    mem_->son(cur_el, 0).vertices_.push_back(cur_el.vertices_[0]);
    mem_->son(cur_el, 0).vertices_.push_back(ptIds[0]);
    mem_->son(cur_el, 0).vertices_.push_back(ptIds[4]);
    mem_->son(cur_el, 0).vertices_.push_back(ptIds[3]);

    mem_->son(cur_el, 1).vertices_.push_back(ptIds[0]);
    mem_->son(cur_el, 1).vertices_.push_back(cur_el.vertices_[1]);
    mem_->son(cur_el, 1).vertices_.push_back(ptIds[1]);
    mem_->son(cur_el, 1).vertices_.push_back(ptIds[4]);

    mem_->son(cur_el, 2).vertices_.push_back(ptIds[4]);
    mem_->son(cur_el, 2).vertices_.push_back(ptIds[1]);
    mem_->son(cur_el, 2).vertices_.push_back(cur_el.vertices_[2]);
    mem_->son(cur_el, 2).vertices_.push_back(ptIds[2]);

    mem_->son(cur_el, 3).vertices_.push_back(ptIds[3]);
    mem_->son(cur_el, 3).vertices_.push_back(ptIds[4]);
    mem_->son(cur_el, 3).vertices_.push_back(ptIds[2]);
    mem_->son(cur_el, 3).vertices_.push_back(cur_el.vertices_[3]);

    updateTopology(elements);

    return;
  }
  //////////////////////////////////////////////////////////////////////////////
};
}  // namespace Bembel
#endif
