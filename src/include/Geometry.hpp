// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef __BEMBEL_GEOMETRY__
#define __BEMBEL_GEOMETRY__

#include "Spline.hpp"

namespace Bembel {

typedef std::vector<Spl::Patch> geometry;

class Geometry {
 public:
  // Constructors
  Geometry() {}
  Geometry(std::string filename) {
    // Note that the Shredder is required. The order of ansatz functions allows
    // to be chosen higher than the smoothness of the NÙRBS mappings. Thus, we
    // need to shredder the geometry mappings to have Bezier patches. You can
    // achieve the higher regularity by changing coefficients in the projector.
    auto tmp = Spl::PatchShredder(Spl::loadGeometryFile(filename));
    assert(Spl::check_geometry(tmp) &&
           "Normals are not all directed the same way");
    _geom = tmp;
  }
  // Geometry(const GeomName &geomname) { init_geometry(geomname); }
  Geometry(Geometry &&other) { _geom = std::move(other._geom); }
  Geometry(const Geometry &other) { _geom = other._geom; }
  Geometry(std::vector<Spl::Patch> in) {
    assert(Spl::check_geometry(in) &&
           "Normals are not all directed the same way");
    _geom = in;
  }
  // Destructors
  //
  // Operators
  Geometry &operator=(Geometry other) {
    std::swap(_geom, other._geom);
    return *this;
  }
  // Methods
  void init_geometry(std::string filename) {
    // Same as constructor
    auto tmp = Spl::PatchShredder(Spl::loadGeometryFile(filename));
    assert(Spl::check_geometry(tmp) &&
           "Normals are not all directed the same way");
    _geom = tmp;
  }
  const geometry &get_geometry() const { return _geom; }
  geometry &get_geometry() { return _geom; }

  void deform(Eigen::Matrix<double, 3, 3> trafo);
  void rotate(double x, double y, double z);
  void stretch(double x, double y, double z);

 private:
  geometry _geom;
};
}  // namespace Bembel
#endif
