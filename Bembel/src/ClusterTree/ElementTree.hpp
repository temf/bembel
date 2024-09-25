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
 * \ingroup ClusterTree
 * \brief This class organizes an element structure on a Geometry object and
 * handles refinement.
 *
 * This class implements a tree which realizes the refinement of the geometry.
 * On the first level of refinement the leafs/elements are all patches of the
 * geometry. Currently, the construction of the tree allows uniform refinement.
 * In the refinement process each element get four sons and the neighborhood
 * relations get updated.
 *
 * The heart of the ElementTree is the leaf iterator which implements a Morton
 * Z-curve. If an element get refined it gets replaced by it four sons element.
 * This core routine can handle adaptive refinement but with some limitations in
 * the resolution of neighborhood relations.
 */
class ElementTree {
 public:
  //////////////////////////////////////////////////////////////////////////////
  /**
   * \brief Deleted copy constructor for the ElementTree class.
   *
   * This copy constructor is explicitly deleted to prevent copying of
   * ElementTree objects.
   */
  ElementTree(const ElementTree &) = delete;

  /**
   * \brief Deleted move constructor for the ElementTree class.
   *
   * This move constructor is explicitly deleted to prevent moving of
   * ElementTree objects.
   */
  ElementTree(ElementTree &&) = delete;

  /**
   * \brief Deleted copy assignment operator for the ElementTree class.
   *
   * This copy assignment operator is explicitly deleted to prevent copying of
   * ElementTree objects.
   *
   * \return A reference to the updated ElementTree object.
   */
  ElementTree &operator=(const ElementTree &) = delete;

  /**
   * \brief Deleted move assignment operator for the ElementTree class.
   *
   * This move assignment operator is explicitly deleted to prevent moving of
   * ElementTree objects.
   *
   * \return A reference to the updated ElementTree object.
   */
  ElementTree &operator=(ElementTree &&) = delete;
  //////////////////////////////////////////////////////////////////////////////
  /// constructors
  //////////////////////////////////////////////////////////////////////////////
  /**
   * \brief Default constructor for the ElementTree class.
   */
  ElementTree() {}
  /**
   * \brief Explicit constructor for the ElementTree class.
   *
   * This constructor initializes an ElementTree object for the provided
   * geometry and maximum level. Currently this constructor is implemented for
   * uniform refinement.
   *
   * \param g The geometry object defining the patches.
   * \param max_level (optional) The maximum level of the ElementTree. Default
   * is 0.
   */
  explicit ElementTree(const Geometry &g, unsigned int max_level = 0) {
    init_ElementTree(g, max_level);
  }
  //////////////////////////////////////////////////////////////////////////////
  /// init
  //////////////////////////////////////////////////////////////////////////////
  /**
   * \brief This function initializes the ElementTree
   *
   * First Each Patch becomes an ElementTreeNode. After resolving the patch
   * neighborhood relations each patch gets refined to the common maximum level
   * of refinement.
   *
   * \param g The geometry object defining the patches.
   * \param max_level (optional) The maximum level of the ElementTree. Default
   * is 0.
   */
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
  /**
   * \brief Refine the given element or it's sons.
   *
   * This function refines an element recursively. If the given element is
   * already refined than iterate all sons and refine them and so on.
   */
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
  /**
   * \brief Refine all patches uniformly
   */
  void refineUniformly() {
    refineUniformly_recursion(root_);
    return;
  }
  //////////////////////////////////////////////////////////////////////////////
  /**
   * \brief Refine a given patch.
   */
  void refinePatch(int patch) {
    refineUniformly_recursion(root_.sons_[patch]);
    return;
  }
  //////////////////////////////////////////////////////////////////////////////
  /**
   * \brief Return the coordinates of all points of the elements.
   *
   * This function iterates all elements and returns the coordinates of the
   * vertices.
   *
   * \param idct Pointer to an Eigen::VectorXi to count how often a vertex is
   *        part of all elements.
   *
   * \return A 3xN Matrix where N is the number of all vertices in the geometry.
   */
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
  /**
   * \brief Return list of elements containing the indices of the vertices.
   *
   * This function stores all indices of an element in the column of the
   * returned matrix.
   *
   * \return A 4xN integer Matrix where N is the number elements.
   */
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
  /**
   * \brief Return number of points in the ElementTree.
   *
   * \return Number of points.
   */
  int get_number_of_points() const { return number_of_points_; }
  /**
   * \brief Return number of elements in the ElementTree.
   *
   * \return Number of elements.
   */
  int get_number_of_elements() const { return number_of_elements_; }
  /**
   * \brief Return maximum level of refinement.
   *
   * \return Level of refinement.
   */
  int get_max_level() const { return max_level_; }
  /**
   * \brief Return reference to the root ElementTreeNode.
   *
   * \return Reference to the root ElementTreeNode.
   */
  ElementTreeNode &root() { return root_; }
  /**
   * \brief Return const reference to the root ElementTreeNode.
   *
   * \return Const Reference to the root ElementTreeNode.
   */
  const ElementTreeNode &root() const { return root_; }
  /**
   * \brief Return const reference to the Geometry.
   *
   * \return Const Reference to the PatchVector.
   */
  const PatchVector &get_geometry() const { return *geometry_; }
  //////////////////////////////////////////////////////////////////////////////
  /// iterators
  //////////////////////////////////////////////////////////////////////////////
  /**
   * \brief Returns an iterator to the beginning of the sequence represented by
   * the leafs as ElementTreeNodes of the ElementTree.
   *
   * \return Returns a ElementTreeNode::const_iterator object.
   */
  ElementTreeNode::const_iterator pbegin() const {
    return ElementTreeNode::const_iterator(pfirst_);
  }
  /**
   * \brief Returns an iterator one past the end of the sequence represented by
   * the leafs as ElementTreeNodes of the ElementTree.
   *
   * \return Returns a ElementTreeNode::const_iterator object.
   */
  ElementTreeNode::const_iterator pend() const {
    return ElementTreeNode::const_iterator(plast_->next_);
  }
  /**
   * \brief Returns an iterator to the beginning of the sequence represented by
   * the leafs as ElementTreeNodes of the ElementTree.
   *
   * \return Returns a ElementTreeNode::const_iterator object.
   */
  ElementTreeNode::const_iterator cpbegin() const { return pbegin(); }
  /**
   * \brief Returns an iterator one past the end of the sequence represented by
   * the leafs as ElementTreeNodes of the ElementTree.
   *
   * \return Returns a ElementTreeNode::const_iterator object.
   */
  ElementTreeNode::const_iterator cpend() const { return pend(); }
  /**
   * \brief Returns a cluster iterator to the beginning of the sequence
   * represented by the the given ElementTreeNode.
   *
   * \param cl Const reference to the ElementTreeNode to start the iterator.
   * \return Returns a ElementTreeNode::const_iterator object.
   */
  ElementTreeNode::const_iterator cluster_begin(
      const ElementTreeNode &cl) const {
    return cl.cbegin();
  }
  /**
   * \brief Returns a cluster iterator one past the end of the sequence
   * represented by the the given ElementTreeNode.
   *
   * \param cl Const reference to the ElementTreeNode to start the iterator.
   * \return Returns a ElementTreeNode::const_iterator object.
   */
  ElementTreeNode::const_iterator cluster_end(const ElementTreeNode &cl) const {
    return cl.cend();
  }
  //////////////////////////////////////////////////////////////////////////////
  /**
   * \brief Computes enclosing balls surrounding all elements.
   *
   * This functions sets the parameters midpoint_ radius_ of the ElementTreeNode
   * objects stored in the tree.
   *
   * \return The point list as 3xN matrix with N the number of points.
   */
  Eigen::MatrixXd computeElementEnclosings() {
    // compute point list
    Eigen::MatrixXd P = generatePointList();
    for (auto i = 0; i < root_.sons_.size(); ++i)
      computeElementEnclosings_recursion(root_.sons_[i], P);
    return P;
  }
  //////////////////////////////////////////////////////////////////////////////
  /**
   * \brief Computes enclosing balls of one element and all its sons.
   *
   * \param el ElementTreeNode to start the recursion.
   * \param P Point list of the vertices.
   */
  void computeElementEnclosings_recursion(ElementTreeNode &el,
                                          const Eigen::MatrixXd &P) {
    Eigen::Vector3d mp1, mp2;
    double r1, r2;
    // compute point list
    if (!el.sons_.size()) {
      // assign enclosing balls to leafs
      util::computeEnclosingBall(&mp1, &r1, P.col(el.vertices_[0]), 0,
                                 P.col(el.vertices_[2]), 0);
      util::computeEnclosingBall(&mp2, &r2, P.col(el.vertices_[1]), 0,
                                 P.col(el.vertices_[3]), 0);
      util::computeEnclosingBall(&(el.midpoint_), &(el.radius_), mp1, r1, mp2,
                                 r2);
    } else {
      // handle the four(!!!) children
      for (auto i = 0; i < 4; ++i)
        computeElementEnclosings_recursion(el.sons_[i], P);
      // assign enclosing balls to fathers bottom up
      util::computeEnclosingBall(&mp1, &r1, el.sons_[0].midpoint_,
                                 el.sons_[0].radius_, el.sons_[2].midpoint_,
                                 el.sons_[2].radius_);
      util::computeEnclosingBall(&mp2, &r2, el.sons_[1].midpoint_,
                                 el.sons_[1].radius_, el.sons_[3].midpoint_,
                                 el.sons_[3].radius_);
      util::computeEnclosingBall(&(el.midpoint_), &(el.radius_), mp1, r1, mp2,
                                 r2);
    }
    return;
  }
  //////////////////////////////////////////////////////////////////////////////
  /**
   * \brief Prints all Elements of the Tree.
   */
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
  /**
   * \brief Return a matrix with all midpoints of the elements.
   *
   * \return The midpoint list as 3xN matrix with N the number of points.
   */
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
  /**
   * \brief Return a matrix with all radii of the element enclosing.
   *
   * \return A vector containing the radii of the element enclosing.
   */
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
  /**
   * \brief Return a vector with computed global indices of the elements.
   *
   * \return A vector with global indices of the elements.
   */
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
  /**
   * \brief Generate list of labels if elements are on the patch boundary or at
   * the boundary of the geometry.
   *
   * The value at the index of the element is 1 if the element is at the patch
   * boundary and there is a neighbor patch. If there is no neighbor patch then
   * the value is -1.
   *
   * \return A vector of integers of size N with N the number of elements.
   */
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
  /**
   * \brief Generate list of labels if elements contained within a given patch.
   *
   * The value at the index of the element is 1 if the element is contained
   * within the given patch. Otherwise its zero.
   *
   * \return A vector of integers of size N with N the number of elements.
   */
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
  /**
   * \brief Resolves neighborhood relations of the patches.
   *
   * \return A vector where each entry defines a patch interface or boundary.
   * The entries correspond to [patchIndex1, edgeCase1, patchIndex2, edgeCase2].
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
  /**
   * \brief The ordering of elements in the element tree does not correspond to
   * the element order underlying the coefficient vector. This reordering can be
   * computed for look ups by this function.
   *
   * Limitation to the uniform case!
   *
   * \return Vector with the tensor product index of the elements.
   */
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
   * \brief function to set up the local topology, i.e. the adjacents
   *        of a refined element
   */
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
  /**
   * \brief This function refines the given ElementTreeNode.
   *
   * The function introduces 4 new elements and takes care of the newly
   * introduces vertices. Furthermore, all the neighborhood relation of all
   * surrounding elements are resolved.
   */
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
