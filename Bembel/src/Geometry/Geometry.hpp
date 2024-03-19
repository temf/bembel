// This file is part of Bembel, the higher order C++ boundary element library.
//
// Copyright (C) 2024 see <http://www.bembel.eu>
//
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef BEMBEL_SRC_GEOMETRY_GEOMETRY_HPP_
#define BEMBEL_SRC_GEOMETRY_GEOMETRY_HPP_

namespace Bembel {
/**
 *  \ingroup Geometry
 *  \brief this class wraps a GeometryVector and provides some basic
 *         functionality, like reading Geometry files
 */
class Geometry {
 public:
  //////////////////////////////////////////////////////////////////////////////
  //    Constructors
  //////////////////////////////////////////////////////////////////////////////
  /**
   * \brief Default constructor.
   */
  Geometry() {}
  /**
   * \brief Constructor.
   *
   * This constructor initializes a geometry specified in a given file. We can
   * handle .dat and .igs at the moment.
   * 
   * \param filename File name of the .dat or .igs file.
   */
  explicit Geometry(const std::string &filename) { init_Geometry(filename); }
  /**
   * \brief Move constructor.
   */
  Geometry(Geometry &&other) { geometry_ = std::move(other.geometry_); }
  /**
   * \brief Copy constructor.
   *
   * Though we are using a shared pointer, we are creating an actual copy here.
   * might be useful if we want to modify the geometry object.
   */
  Geometry(const Geometry &other) {
    geometry_ = std::make_shared<PatchVector>();
    *geometry_ = *(other.geometry_);
  }
  /**
   * \brief Copy constructor.
   *
   * Though we are using a shared pointer, we are creating an actual copy here.
   * might be useful if we want to modify the geometry object.
   *
   * \param in Geometry provided as PatchVector.
   */
  explicit Geometry(const PatchVector &in) {
    geometry_ = std::make_shared<PatchVector>();
    *geometry_ = in;
  }
  /**
   * \brief Copy assignment operator for the Geometry class.
   *
   * This copy assignment operator is explicitly deleted to prevent copying of
   * Geometry objects.
   *
   * \return A reference to the updated Geometry object.
   */
  Geometry &operator=(Geometry other) {
    std::swap(geometry_, other.geometry_);
    return *this;
  }
  //////////////////////////////////////////////////////////////////////////////
  //    init_Geometry
  //////////////////////////////////////////////////////////////////////////////
  /**
   * \brief Initialize the geometry from a geometry given by a file.
   *
   * Note that the Shredder is required. The order of ansatz functions allows to
   * be chosen higher than the smoothness of the NÙRBS mappings. Thus, we need
   * to shredder the geometry mappings to have Bezier patches. You can achieve
   * the higher regularity by changing coefficients in the projector.
   *
   * \param filename Filename of a .dat or .igs file.
   */
  inline void init_Geometry(const std::string &filename) {
    std::string file_suffix = filename.substr(filename.find('.') + 1);
    PatchVector tmp;

    if (file_suffix.compare("dat") == 0)
      tmp = PatchShredder(LoadGeometryFileDAT(filename));
    else if (file_suffix.compare("igs") == 0)
      tmp = PatchShredder(LoadGeometryFileIGS(filename));
    else
      assert(!"File type unknown!");

    geometry_ = std::make_shared<PatchVector>();
    *geometry_ = tmp;
  }
  //////////////////////////////////////////////////////////////////////////////
  //    getters
  //////////////////////////////////////////////////////////////////////////////
  /**
   * \brief Return const reference to the geometry.
   * 
   * \return Const reference to the PatchVector.
   */
  const PatchVector &get_geometry() const { return *geometry_; }
  /**
   * \brief Return reference to the geometry.
   * 
   * \return Reference to the PatchVector.
   */
  PatchVector &get_geometry() { return *geometry_; }
  /**
   * \brief Return const pointer to the geometry.
   * 
   * \return Const shared pointer to the PatchVector.
   */
  const std::shared_ptr<PatchVector> get_geometry_ptr() const {
    return geometry_;
  }
  /**
   * \brief Return pointer to the geometry.
   * 
   * \return Shared pointer to the PatchVector.
   */
  std::shared_ptr<PatchVector> get_geometry_ptr() { return geometry_; }

  /**
   * \brief Get number of patches.
   * 
   * \return Number of patches.
   */
  int get_number_of_patches() const { return geometry_->size(); }
  //////////////////////////////////////////////////////////////////////////////
  //    private member variables
  //////////////////////////////////////////////////////////////////////////////
 private:
  std::shared_ptr<PatchVector> geometry_;
};
}  // namespace Bembel
#endif  // BEMBEL_SRC_GEOMETRY_GEOMETRY_HPP_
