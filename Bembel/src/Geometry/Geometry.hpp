// This file is part of Bembel, the higher order C++ boundary element library.
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
  Geometry() {}
  Geometry(const std::string &filename) { init_Geometry(filename); }
  Geometry(Geometry &&other) { geometry_ = std::move(other.geometry_); }
  // though we are using a shared pointer, we are creating an actual
  // copy here. might be useful if we want to modify the geometry object
  Geometry(const Geometry &other) {
    geometry_ = std::make_shared<PatchVector>();
    *geometry_ = *(other.geometry_);
  }
  Geometry(const PatchVector &in) {
    geometry_ = std::make_shared<PatchVector>();
    *geometry_ = in;
  }
  Geometry &operator=(Geometry other) {
    std::swap(geometry_, other.geometry_);
    return *this;
  }
  //////////////////////////////////////////////////////////////////////////////
  //    init_Geometry
  //////////////////////////////////////////////////////////////////////////////
  inline void init_Geometry(const std::string &filename) {
    // Note that the Shredder is required. The order of ansatz functions allows
    // to be chosen higher than the smoothness of the NÙRBS mappings. Thus, we
    // need to shredder the geometry mappings to have Bezier patches. You can
    // achieve the higher regularity by changing coefficients in the projector.
    auto tmp = Bembel::PatchShredder(Bembel::LoadGeometryFile(filename));
    geometry_ = std::make_shared<PatchVector>();
    *geometry_ = tmp;
  }
  //////////////////////////////////////////////////////////////////////////////
  //    getters
  //////////////////////////////////////////////////////////////////////////////
  const PatchVector &get_geometry() const { return *geometry_; }
  PatchVector &get_geometry() { return *geometry_; }
  const std::shared_ptr<PatchVector> get_geometry_ptr() const {
    return geometry_;
  }
  std::shared_ptr<PatchVector> get_geometry_ptr() { return geometry_; }
  int get_number_of_patches() const { return geometry_->size(); };
  //////////////////////////////////////////////////////////////////////////////
  //    private member variables
  //////////////////////////////////////////////////////////////////////////////
 private:
  std::shared_ptr<PatchVector> geometry_;
};
}  // namespace Bembel
#endif  // BEMBEL_SRC_GEOMETRY_GEOMETRY_HPP_
