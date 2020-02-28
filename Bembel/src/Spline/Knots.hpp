// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef BEMBEL_SPLINE_KNOTS_H_
#define BEMBEL_SPLINE_KNOTS_H_

namespace Bembel {
namespace Spl {
/**
 *  \ingroup Spline
 *  \brief Here, routines for the creation and processing of knot vectors are
 *         defined.
 */
inline std::vector<double> MakeBezierKnotVector(
    int polynomial_degree) noexcept {
  std::vector<double> out;
  out.reserve(polynomial_degree * 2);
  for (int i = 0; i < polynomial_degree; i++) out.push_back(0);
  for (int i = 0; i < polynomial_degree; i++) out.push_back(1);
  return out;
}

inline int GetPolynomialDegreeFromKnotVector(
    const std::vector<double> &knot_vector) {
  constexpr double tol = .0000001;
  int i = 0;
  const int sz = knot_vector.size();
  while (++i < sz) {
    if (knot_vector[i] > tol) {
      return i;
    };
  }
  return 0;
}

inline std::vector<double> MakeUniformKnotVector(int polynomial_degree,
                                                 int interior,
                                                 int knotrepetition) noexcept {
  std::vector<double> out;
  const double h = 1. / (interior + 1);
  out.reserve(polynomial_degree * 2);
  for (int i = 0; i < polynomial_degree; i++) out.push_back(0);
  for (int i = 1; i < interior + 1; i++)
    for (int k = 0; k < knotrepetition; k++) out.push_back(i * h);
  for (int i = 0; i < polynomial_degree; i++) out.push_back(1);
  return out;
}

inline std::vector<double> MakeUniformKnotVector(int p, int lvl) {
  return MakeUniformKnotVector(p, lvl, 1);
}

// Chunk knots eats knotvectors and returns knotvectors in which each knot is
// unique, up to tolerance tol.
inline std::vector<double> ExtractUniqueKnotVector(
    const std::vector<double> &in) {
  constexpr double tol = Bembel::Constants::generic_tolerance;
  std::vector<double> out;
  const int size = in.size();
  // std::cout<< "chunk " << in[0]  <<"\n";
  out.push_back(in[0]);
  for (int i = 1; i < size; i++) {
    if (in[i] > (in[i - 1] + tol)) out.push_back(in[i]);
  }
  return out;
}

// Up tp tol, fins knot element index for a given x.
inline int FindLocationInKnotVector(const double &x,
                                    const std::vector<double> &v) {
  constexpr double tol = Bembel::Constants::generic_tolerance;
  int size = v.size();
  if (((x - tol) < 0.) && ((x + tol) > 0.)) return 0;
  for (int i = 0; i < size - 1; i++) {
    if (v[i] <= x && v[i + 1] > x) return i;
  }
  return size - 2;
}
}  // namespace Spl
}  // namespace Bembel
#endif
