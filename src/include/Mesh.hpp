// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef __BEMBEL_MESH__
#define __BEMBEL_MESH__

#include "meshdata.h"
namespace Bembel {
class Mesh {
 public:
  // Constructors
  Mesh() {}
  Mesh(const Geometry &geom, int M) { init_mesh(geom, M); }
  // Destructors
  ~Mesh() { free_meshdata(&_mesh); }
  // Methods
  void init_mesh(const Geometry &geom, int M) {
    _mesh = get_meshdata(geom.get_geometry(), M);
  }
  const meshdata &get_mesh() const { return _mesh; }
  meshdata &get_mesh() { return _mesh; }

 private:
  // we declare functionality which has not been implemented (yet)
  // to be private
  Mesh(const Mesh& other);
  Mesh(Mesh&& other);
  Mesh& operator=(const Mesh& other);
  Mesh& operator=(Mesh&& other);
  meshdata _mesh;
};
}  // namespace Bembel
#endif
